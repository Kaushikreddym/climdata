"""Read CMIP6 data from a local ESGF-style mirror.

This is a *reader*, not a provider: it fetches nothing over the network, it
indexes and opens NetCDF files that are already on disk. That is why it lives
under ``climdata.utils`` rather than ``climdata.datasets`` — it has no
:meth:`download` to offer and does not belong in the dataset registry.

Two entry points:

* :func:`index_cmip6` walks a mirror and returns a pandas catalogue — one row
  per NetCDF file, with the DRS facets (variable, table, model, experiment,
  member, grid, version) and the time span each file covers.
* :func:`open_cmip6` picks one model/experiment out of that catalogue and hands
  back a single :class:`xarray.Dataset`, longitudes normalised to -180..180 and
  latitudes ascending, optionally cropped in space and time.

:func:`cmip6_coverage` summarises a catalogue as a variable × experiment table,
which is the quickest way to see which models are complete enough to use.

The mirrors are messy in ways a strict DRS parser trips over, so the facets are
read from the *filename* (which is canonical) rather than from fixed path
depths:

* macOS AppleDouble sidecars (``._foo.nc``) sit next to the real files and are
  not valid NetCDF — skipped.
* Version directories can carry the ESGF host they came from, split across
  extra directory levels (``v20190903esgf/ceda/ac/uk/``), so a file can be 3
  levels deeper than the DRS says.
* The same version is often mirrored from two hosts, and some models publish
  more than one grid (``gr1`` and ``gr2``), so a naive open of everything
  matching a variable yields duplicated or conflicting time steps.
  :func:`open_cmip6` resolves grid, version and member down to one set of files
  before opening anything.

Example:
    >>> from climdata.utils.cmip6_local import index_cmip6, open_cmip6
    >>> df = index_cmip6(table_id="day")                      # doctest: +SKIP
    >>> sorted(df.source_id.unique())                          # doctest: +SKIP
    ['ACCESS-ESM1-5', 'IPSL-CM6A-LR', 'MIROC6', 'MRI-ESM2-0']
    >>> ds = open_cmip6("MIROC6", "hist-GHG",                  # doctest: +SKIP
    ...                 ["pr", "tas"], index=df)
"""

import calendar
import re
import warnings
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Union

import pandas as pd
import xarray as xr

__all__ = ["index_cmip6", "open_cmip6", "cmip6_coverage", "DEFAULT_BASE_DIR"]

DEFAULT_BASE_DIR = "/data01/FDS/muduchuru/Atmos/cmip6_data"

# {variable}_{table}_{model}_{experiment}_{member}_{grid}[_{start}-{end}].nc
# The time range is absent for fixed fields (e.g. orography), hence optional.
_FILENAME_RE = re.compile(
    r"^(?P<variable_id>[^_]+)_(?P<table_id>[^_]+)_(?P<source_id>[^_]+)_"
    r"(?P<experiment_id>[^_]+)_(?P<member_id>[^_]+)_(?P<grid_label>[^_]+)"
    r"(?:_(?P<time_range>\d+-\d+))?\.nc$"
)

# A version directory starts with vYYYYMMDD; anything after that is the ESGF
# host glued on by whatever copied the mirror.
_VERSION_RE = re.compile(r"^(v\d+)")

# Native grid first, then regridded ones — used to break a grid_label tie.
_GRID_PREFERENCE = ("gn", "gr", "gr1", "gr2", "gr3")

# Dates xarray can hold as numpy datetime64[ns]; outside it, it uses cftime.
_NS_RANGE = (pd.Timestamp.min.ceil("D"), pd.Timestamp.max.floor("D"))

_INDEX_COLUMNS = [
    "path", "activity_id", "institution_id", "source_id", "experiment_id",
    "member_id", "table_id", "variable_id", "grid_label", "version",
    "start_date", "end_date", "years",
]


def _as_list(value) -> Optional[List[str]]:
    """Normalise a str-or-iterable filter argument to a list (None passes through)."""
    if value is None:
        return None
    if isinstance(value, str):
        return [value]
    return list(value)


def _stamp_to_timestamp(stamp: str, last: bool):
    """One DRS date stamp -> Timestamp, at the first or last day it covers.

    ``pd.Period`` is not usable here: its ``end_time`` is nanosecond-based and
    overflows on the SSP extension files that run to 2300, so the end of a
    month/year stamp is worked out directly instead.
    """
    year = int(stamp[:4])
    if len(stamp) == 8:
        return pd.Timestamp(year=year, month=int(stamp[4:6]), day=int(stamp[6:8]))
    if len(stamp) == 6:
        month = int(stamp[4:6])
        day = calendar.monthrange(year, month)[1] if last else 1
        return pd.Timestamp(year=year, month=month, day=day)
    if len(stamp) == 4:
        return pd.Timestamp(year=year, month=12, day=31) if last \
            else pd.Timestamp(year=year, month=1, day=1)
    raise ValueError(f"unrecognised DRS date stamp: {stamp!r}")


def _parse_time_range(time_range: Optional[str]):
    """DRS time range -> (start, inclusive end) timestamps.

    Handles the daily/monthly/yearly stampings (YYYYMMDD, YYYYMM, YYYY). The end
    is the *last day covered*, so a file stamped ``185001-185012`` ends
    1850-12-31 rather than 1850-12-01 — otherwise coverage totals systematically
    undercount.
    """
    if not time_range:
        return pd.NaT, pd.NaT
    try:
        start_str, end_str = time_range.split("-")
        return (_stamp_to_timestamp(start_str, last=False),
                _stamp_to_timestamp(end_str, last=True))
    except (ValueError, pd.errors.OutOfBoundsDatetime):
        return pd.NaT, pd.NaT


def _facets_from_path(path: Path, source_id: str, grid_label: str) -> Dict[str, Optional[str]]:
    """Recover activity_id, institution_id and version from the directory tree.

    The filename has no activity or institution, and the version directory can
    be followed by extra host levels, so both are located by anchoring on facets
    the filename *does* carry rather than by counting from the mirror root.
    """
    parts = path.parts
    activity_id = institution_id = version = None

    # .../<activity_id>/<institution_id>/<source_id>/<experiment_id>/...
    if source_id in parts:
        i = parts.index(source_id)
        if i >= 2:
            activity_id, institution_id = parts[i - 2], parts[i - 1]

    # .../<grid_label>/<version>[/<esgf host levels>]/<file>.nc — search from the
    # right, since a grid_label like 'gn' could also appear earlier by accident.
    if grid_label in parts:
        j = len(parts) - 1 - parts[::-1].index(grid_label)
        if j + 1 < len(parts) - 1:
            m = _VERSION_RE.match(parts[j + 1])
            version = m.group(1) if m else parts[j + 1]

    return {"activity_id": activity_id, "institution_id": institution_id,
            "version": version}


def index_cmip6(
    base_dir: Union[str, Path, Sequence[Union[str, Path]]] = DEFAULT_BASE_DIR,
    *,
    table_id: Optional[Union[str, Iterable[str]]] = None,
    variables: Optional[Union[str, Iterable[str]]] = None,
    experiments: Optional[Union[str, Iterable[str]]] = None,
    models: Optional[Union[str, Iterable[str]]] = None,
    members: Optional[Union[str, Iterable[str]]] = None,
    require_all_experiments: bool = False,
    min_years: Optional[float] = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """Index a local CMIP6 mirror into a catalogue of NetCDF files.

    Args:
        base_dir (str, Path or sequence): Mirror root(s) to walk. A sequence
            indexes several mirrors into one catalogue, which is useful when a
            model lives under one root and the rest under another.
        table_id (str or iterable, optional): Keep only these CMIP tables
            (e.g. 'day', 'Amon'). None keeps all.
        variables (str or iterable, optional): Keep only these variable ids.
        experiments (str or iterable, optional): Keep only these experiment ids.
        models (str or iterable, optional): Keep only these source ids.
        members (str or iterable, optional): Keep only these member ids.
        require_all_experiments (bool): Drop any model that does not have every
            experiment in ``experiments``. Requires ``experiments``.
        min_years (float, optional): Drop any (model, experiment, member,
            variable, grid) group whose files cover fewer than this many years
            in total.
        verbose (bool): Print a one-line summary and any skipped-file count.

    Returns:
        pd.DataFrame: One row per file with columns ``path``, ``activity_id``,
        ``institution_id``, ``source_id``, ``experiment_id``, ``member_id``,
        ``table_id``, ``variable_id``, ``grid_label``, ``version``,
        ``start_date``, ``end_date``, ``years``. Empty (but with those columns)
        when nothing matches.

    Raises:
        FileNotFoundError: If a ``base_dir`` does not exist.
        ValueError: If ``require_all_experiments`` is set without ``experiments``.
    """
    bases = [base_dir] if isinstance(base_dir, (str, Path)) else list(base_dir)
    bases = [Path(b) for b in bases]
    for b in bases:
        if not b.is_dir():
            raise FileNotFoundError(f"CMIP6 mirror root not found: {b}")

    table_id    = _as_list(table_id)
    variables   = _as_list(variables)
    experiments = _as_list(experiments)
    models      = _as_list(models)
    members     = _as_list(members)
    if require_all_experiments and not experiments:
        raise ValueError("require_all_experiments needs `experiments` to be given.")

    rows, n_skipped = [], 0
    for base in bases:
        for path in base.rglob("*.nc"):
            # AppleDouble sidecars from a macOS copy: not NetCDF, same name.
            if path.name.startswith("._"):
                continue
            match = _FILENAME_RE.match(path.name)
            if match is None:
                n_skipped += 1          # not DRS-named (e.g. a postproc/ file)
                continue
            f = match.groupdict()

            # Cheap facet filters before touching the path, since the filename
            # already carries everything they test.
            if table_id    and f["table_id"]      not in table_id:    continue
            if variables   and f["variable_id"]   not in variables:   continue
            if experiments and f["experiment_id"] not in experiments: continue
            if models      and f["source_id"]     not in models:      continue
            if members     and f["member_id"]     not in members:     continue

            start, end = _parse_time_range(f.pop("time_range"))
            rows.append({
                "path": str(path),
                **_facets_from_path(path, f["source_id"], f["grid_label"]),
                **f,
                "start_date": start,
                "end_date": end,
            })

    df = pd.DataFrame(rows)
    if df.empty:
        if verbose:
            print(f"No CMIP6 files matched under {', '.join(str(b) for b in bases)}"
                  f" ({n_skipped} non-DRS files skipped)")
        return pd.DataFrame(columns=_INDEX_COLUMNS)

    # Second-resolution datetimes, not the default nanosecond ones: SSP
    # extension files run past 2262 and would otherwise become NaT.
    for col in ("start_date", "end_date"):
        df[col] = pd.Series(df[col].tolist(), dtype="datetime64[s]")
    df["years"] = (df["end_date"] - df["start_date"]).dt.days / 365.25

    if require_all_experiments:
        have = df.groupby("source_id")["experiment_id"].agg(set)
        keep = have[have.apply(lambda e: set(experiments).issubset(e))].index
        dropped = sorted(set(df["source_id"]) - set(keep))
        df = df[df["source_id"].isin(keep)]
        if verbose and dropped:
            print(f"  dropped {len(dropped)} model(s) missing one of {experiments}: "
                  f"{', '.join(dropped)}")

    if min_years is not None and not df.empty:
        group = ["source_id", "experiment_id", "member_id", "variable_id", "grid_label"]
        total = df.groupby(group)["years"].transform("sum")
        df = df[total >= min_years]

    df = df[_INDEX_COLUMNS].sort_values(
        ["source_id", "experiment_id", "member_id", "variable_id", "start_date"]
    ).reset_index(drop=True)

    if verbose:
        print(f"Indexed {len(df)} CMIP6 files: "
              f"{df['source_id'].nunique()} model(s), "
              f"{df['experiment_id'].nunique()} experiment(s), "
              f"{df['variable_id'].nunique()} variable(s)"
              + (f" ({n_skipped} non-DRS files skipped)" if n_skipped else ""))
    return df


def cmip6_coverage(index: pd.DataFrame, value: str = "years") -> pd.DataFrame:
    """Summarise a catalogue as a (model, experiment) × variable table.

    Reports the coverage of the *best single* member/grid combination, not the
    sum over all of them: a mirror that holds the same version twice (two ESGF
    hosts) or two published grids would otherwise show double the years that are
    actually readable.

    Args:
        index (pd.DataFrame): Catalogue from :func:`index_cmip6`.
        value (str): 'years' for total years covered, 'files' for a file count.

    Returns:
        pd.DataFrame: Rows are (source_id, experiment_id), columns are
        variable ids, zeros where a model has nothing. A model row-block with no
        zeros and no short entries is one that can be used as-is.
    """
    if index.empty:
        return pd.DataFrame()
    if value not in ("years", "files"):
        raise ValueError(f"value must be 'years' or 'files', got {value!r}")

    # One row per distinct time slice: collapses the same file mirrored from two
    # hosts under different version directories.
    unique = index.drop_duplicates(
        subset=["source_id", "experiment_id", "member_id", "variable_id",
                "grid_label", "start_date", "end_date"])

    agg = ("path", "count") if value == "files" else ("years", "sum")
    per_choice = (unique.groupby(["source_id", "experiment_id", "member_id",
                                  "variable_id", "grid_label"])
                        .agg(**{value: agg})
                        .reset_index())
    return (per_choice.groupby(["source_id", "experiment_id", "variable_id"])[value]
                      .max()                       # best member/grid, not their sum
                      .reset_index()
                      .pivot(index=["source_id", "experiment_id"],
                             columns="variable_id", values=value)
                      .fillna(0))


def _select_files(index: pd.DataFrame, source_id: str, experiment_id: str,
                  variable_id: str, member_id: Optional[str],
                  grid_label: Optional[str]) -> pd.DataFrame:
    """Narrow a catalogue to the one member/grid/version to actually open.

    Mirrors routinely hold the same variable more than once — the same version
    pulled from two ESGF hosts, or two published grids (gr1/gr2). Opening all of
    them gives duplicated or conflicting time steps, so this collapses the
    choice before any file is read.
    """
    sub = index[(index["source_id"] == source_id) &
                (index["experiment_id"] == experiment_id) &
                (index["variable_id"] == variable_id)]
    if sub.empty:
        return sub

    if member_id is not None:
        sub = sub[sub["member_id"] == member_id]
    else:
        available = sorted(sub["member_id"].unique())
        pick = "r1i1p1f1" if "r1i1p1f1" in available else available[0]
        sub = sub[sub["member_id"] == pick]

    if grid_label is not None:
        sub = sub[sub["grid_label"] == grid_label]
    elif sub["grid_label"].nunique() > 1:
        available = sorted(sub["grid_label"].unique())
        pick = next((g for g in _GRID_PREFERENCE if g in available), available[0])
        sub = sub[sub["grid_label"] == pick]

    # Resolve versions *per time slice*, not globally: two versions often cover
    # disjoint periods (2015-2100 published separately from the 2101-2300
    # extension), so keeping only the newest version would silently drop years.
    # Within one time slice the newest version wins, which also collapses the
    # same file mirrored from two ESGF hosts.
    sub = sub.sort_values("version", na_position="first")
    return sub.drop_duplicates(subset=["start_date", "end_date"], keep="last")


def _normalise_grid(ds: xr.Dataset) -> xr.Dataset:
    """Longitudes to -180..180 and both axes ascending, so ``sel`` with a plain
    slice behaves the same for every model."""
    if "lon" in ds.coords and float(ds.lon.max()) > 180:
        ds = ds.assign_coords(lon=((ds.lon + 180) % 360) - 180).sortby("lon")
    for axis in ("lat", "lon"):
        if axis in ds.coords and ds[axis].size > 1 and ds[axis][0] > ds[axis][-1]:
            ds = ds.sortby(axis)
    return ds


def open_cmip6(
    source_id: str,
    experiment_id: str,
    variables: Union[str, Iterable[str]],
    *,
    base_dir: Union[str, Path, Sequence[Union[str, Path]]] = DEFAULT_BASE_DIR,
    index: Optional[pd.DataFrame] = None,
    table_id: str = "day",
    member_id: Optional[str] = None,
    grid_label: Optional[str] = None,
    bounds: Optional[Dict[str, float]] = None,
    time_range: Optional[Sequence[Optional[str]]] = None,
    chunks: Optional[Dict[str, int]] = None,
    use_cftime: Optional[bool] = None,
    strict: bool = False,
    verbose: bool = True,
) -> xr.Dataset:
    """Open one model/experiment from a local CMIP6 mirror as a Dataset.

    Args:
        source_id (str): Model, e.g. 'MIROC6'.
        experiment_id (str): Experiment, e.g. 'historical', 'hist-GHG'.
        variables (str or iterable): Variable id(s) to read.
        base_dir (str, Path or sequence): Mirror root(s). Ignored when ``index``
            is given.
        index (pd.DataFrame, optional): A catalogue from :func:`index_cmip6`,
            to avoid re-walking the mirror on every call.
        table_id (str): CMIP table to read (default 'day'). Ignored when
            ``index`` is given — filter the catalogue instead.
        member_id (str, optional): Ensemble member. Default: 'r1i1p1f1' when
            present, else the first sorted.
        grid_label (str, optional): Grid to read when the model publishes more
            than one. Default: native ('gn') if present, else the first sorted.
        bounds (dict, optional): Spatial crop with keys ``lat_min``, ``lat_max``,
            ``lon_min``, ``lon_max``. Applied after the -180..180 normalisation.
        time_range (sequence, optional): ``(start, end)`` passed to
            ``sel(time=slice(...))``; either may be None for an open end.
        chunks (dict, optional): Dask chunks (default ``{'time': 365}``).
        use_cftime (bool, optional): Force cftime (True) or numpy datetime64
            (False) time decoding. Default None picks cftime when the selected
            files reach outside the datetime64[ns] range (an SSP extension to
            2300, a piControl stamped from year 1), which is the case where
            mixed decoding would otherwise break the concatenation.
        strict (bool): Raise if any requested variable is missing, instead of
            warning and returning the rest.
        verbose (bool): Print what was chosen for each variable.

    Returns:
        xr.Dataset: The requested variables on the model's native grid, in
        native units, with ``source_id``/``experiment_id``/``member_id``/
        ``grid_label`` recorded in ``.attrs``.

    Raises:
        FileNotFoundError: If no requested variable is present, or if ``strict``
            and any single one is missing.
    """
    variables = _as_list(variables)
    if index is None:
        index = index_cmip6(base_dir, table_id=table_id, models=source_id,
                            experiments=experiment_id, variables=variables,
                            verbose=verbose)
    if index.empty:
        raise FileNotFoundError(
            f"No CMIP6 files for {source_id}/{experiment_id} under {base_dir}."
        )

    # Resolve every variable's files before opening any of them: whether to
    # decode times as cftime has to be decided once for the whole call.
    selected = {}
    for var in variables:
        sub = _select_files(index, source_id, experiment_id, var, member_id, grid_label)
        if not sub.empty:
            selected[var] = sub

    if use_cftime is None:
        # xarray falls back to cftime per *file* when its dates leave the
        # datetime64[ns] range, so an SSP whose extension runs to 2300 decodes
        # its 2015-2100 files as Timestamps and its 2251-2300 files as cftime —
        # and combine_by_coords then refuses the mixed time axis. Force cftime
        # for the whole call when any selected file reaches outside the range.
        #
        # Only ever forced to True, never to False: a model on a non-standard
        # calendar (GFDL-ESM4 is noleap, others are 360_day) has no datetime64
        # representation at all, and an explicit False makes decoding fail
        # outright. Leaving it None lets xarray pick cftime for those.
        spans = pd.concat([s[["start_date", "end_date"]] for s in selected.values()]) \
            if selected else pd.DataFrame(columns=["start_date", "end_date"])
        out_of_range = bool(len(spans)) and bool(
            (spans["end_date"] > _NS_RANGE[1]).any() or
            (spans["start_date"] < _NS_RANGE[0]).any()
        )
        use_cftime = True if out_of_range else None
        if out_of_range and verbose:
            print("  (dates outside the datetime64[ns] range — decoding times as cftime)")

    data_vars, missing, chosen = {}, [], {}
    for var in variables:
        sub = selected.get(var)
        if sub is None:
            missing.append(var)
            continue
        paths = sorted(sub["path"].tolist())
        ds = xr.open_mfdataset(paths, combine="by_coords",
                               chunks=chunks or {"time": 365},
                               use_cftime=use_cftime)
        if var not in ds:
            missing.append(var)
            ds.close()
            continue
        data_vars[var] = ds[var]
        chosen[var] = (sub["member_id"].iloc[0], sub["grid_label"].iloc[0],
                       sub["version"].iloc[0], len(paths))
        if verbose:
            m, g, v, n = chosen[var]
            print(f"  [{var}] {n} file(s)  member={m} grid={g} version={v}")

    if missing:
        msg = (f"{source_id}/{experiment_id}: no files for {missing} "
               f"(present: {sorted(index['variable_id'].unique())})")
        if strict or not data_vars:
            raise FileNotFoundError(msg)
        warnings.warn(msg)

    ds = _normalise_grid(xr.Dataset(data_vars))

    if bounds:
        ds = ds.sel(lat=slice(bounds["lat_min"], bounds["lat_max"]),
                    lon=slice(bounds["lon_min"], bounds["lon_max"]))
    if time_range:
        ds = ds.sel(time=slice(time_range[0], time_range[1]))

    first = next(iter(chosen.values()))
    ds.attrs.update({
        "source_id": source_id,
        "experiment_id": experiment_id,
        "member_id": first[0],
        "grid_label": first[1],
        "version": first[2] or "",
    })
    return ds
