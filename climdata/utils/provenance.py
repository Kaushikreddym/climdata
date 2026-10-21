"""Provenance records for climdata runs and the FAIR crate built from them.

:class:`~climdata.utils.wrapper_workflow.ClimateExtractor` logs one entry per
processing step (extract, calc_index, impute, reproject, to_nc, ...) into
``extractor.provenance``. :meth:`ClimateExtractor.to_nc` and
:meth:`ClimateExtractor.to_csv` write their file into a crate folder and turn
that log into three files that travel with the data:

.. code-block:: text

    <stem>_crate/
        <stem>.nc / <stem>.csv   # the data; the NetCDF's CF ``history`` lists every step
        meta.json                # full record: config, steps, software, checksums
        ro-crate-metadata.json   # RO-Crate 1.1 manifest, one CreateAction per step
        ro-crate-preview.html    # human-readable summary of all of the above

``meta.json`` is the source of truth; the RO-Crate manifest and the HTML preview
are both rendered from it, so the three cannot disagree.

Everything here is pure: functions take plain Python / xarray objects and
return dicts or strings, so they are testable without a configured extractor.
"""

import datetime
import hashlib
import html
import json
import platform
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
import xarray as xr

#: Version of the ``meta.json`` layout. Bump on incompatible changes.
META_SCHEMA_VERSION = "1.0"

#: Config groups composed in from ``conf/mappings`` — whole catalogues, not
#: choices made for this run, so they are trimmed from the config snapshot.
_CATALOGUE_KEYS = ("dsinfo", "varinfo", "extinfo", "actionsinfo", "imputeinfo", "bounds")

#: Packages whose versions decide what a run produces.
_TRACKED_PACKAGES = (
    "climdata", "xarray", "numpy", "pandas", "xclim", "dask", "netCDF4",
    "h5netcdf", "rioxarray", "rasterio", "pyproj", "xesmf", "hydra-core",
)


def now_iso() -> str:
    """Current UTC time as ``YYYY-MM-DDTHH:MM:SSZ``."""
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def jsonable(obj: Any) -> Any:
    """Convert ``obj`` into something :func:`json.dumps` accepts.

    Handles OmegaConf containers, numpy scalars and arrays, paths, timestamps
    and tuples. Anything else unknown is stringified rather than dropped, so a
    provenance record never fails to serialise over one odd parameter.
    """
    try:
        from omegaconf import DictConfig, ListConfig, OmegaConf
        if isinstance(obj, (DictConfig, ListConfig)):
            return jsonable(OmegaConf.to_container(obj, resolve=True))
    except ImportError:  # pragma: no cover - omegaconf is a hard dependency
        pass
    if obj is None or isinstance(obj, (bool, int, str)):
        return obj
    if isinstance(obj, float):
        return obj if np.isfinite(obj) else str(obj)
    if isinstance(obj, dict):
        return {str(k): jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set)):
        return [jsonable(v) for v in obj]
    if isinstance(obj, np.generic):
        return jsonable(obj.item())
    if isinstance(obj, np.ndarray):
        return jsonable(obj.tolist())
    if isinstance(obj, (pd.Timestamp, datetime.datetime, datetime.date, np.datetime64)):
        return str(pd.Timestamp(obj).isoformat())
    if isinstance(obj, Path):
        return str(obj)
    return str(obj)


def sha256(path) -> str:
    """Hex SHA-256 of a file, read in 1 MiB blocks."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def file_record(path, name: Optional[str] = None) -> Dict[str, Any]:
    """Describe a file (or a directory of files) on disk.

    Args:
        path: File or directory.
        name (str, optional): Name to record instead of the absolute path —
            used for files inside a crate, which are referenced relatively.

    Returns:
        dict: ``path``, ``size_bytes`` and ``sha256`` for a file; for a
        directory, ``path`` and a ``files`` list of the same for each file.
    """
    p = Path(path)
    if p.is_dir():
        return {
            "path": name or str(p.resolve()),
            "files": [file_record(f) for f in sorted(p.rglob("*")) if f.is_file()],
        }
    return {
        "path": name or str(p.resolve()),
        "size_bytes": p.stat().st_size,
        "sha256": sha256(p),
    }


def _coord(ds: xr.Dataset, names):
    for n in names:
        if n in ds.coords or n in ds.variables:
            return ds[n]
    return None


def temporal_coverage(ds: xr.Dataset) -> Optional[Dict[str, str]]:
    """First and last time step of ``ds``, or ``None`` without a time axis."""
    t = _coord(ds, ("time", "date"))
    if t is None or t.size == 0:
        return None
    vals = t.values.reshape(-1)
    try:
        start, end = pd.to_datetime(vals.min()), pd.to_datetime(vals.max())
        return {"start": start.strftime("%Y-%m-%d"), "end": end.strftime("%Y-%m-%d")}
    except (TypeError, ValueError):  # cftime calendars
        return {"start": str(vals.min())[:10], "end": str(vals.max())[:10]}


def spatial_coverage(ds: xr.Dataset) -> Optional[Dict[str, float]]:
    """Lat/lon bounding box of ``ds``, from 1-D or 2-D coordinates.

    Returns ``None`` when the dataset carries no geographic coordinates (a
    projected grid without auxiliary ``lat``/``lon``).
    """
    lat = _coord(ds, ("lat", "latitude"))
    lon = _coord(ds, ("lon", "longitude"))
    if lat is None or lon is None:
        return None
    lat_v = np.asarray(lat.values, dtype=float)
    lon_v = np.asarray(lon.values, dtype=float)
    return {
        "lat_min": round(float(np.nanmin(lat_v)), 6),
        "lat_max": round(float(np.nanmax(lat_v)), 6),
        "lon_min": round(float(np.nanmin(lon_v)), 6),
        "lon_max": round(float(np.nanmax(lon_v)), 6),
    }


def summarize_dataset(ds: xr.Dataset) -> Dict[str, Any]:
    """Structural summary of a dataset: variables, dimensions, extent, CRS.

    Reads only coordinates and attributes, never data values, so it is cheap
    on a lazy, Dask-backed dataset.
    """
    variables = []
    for name, da in ds.data_vars.items():
        variables.append({
            "name": name,
            "long_name": da.attrs.get("long_name", name),
            "standard_name": da.attrs.get("standard_name"),
            "units": da.attrs.get("units"),
            "dims": list(da.dims),
            "shape": list(da.shape),
            "dtype": str(da.dtype),
        })
    crs = None
    if "spatial_ref" in ds.variables:
        crs = ds["spatial_ref"].attrs.get("crs_wkt") or ds["spatial_ref"].attrs.get("spatial_ref")
    return {
        "type": "xarray.Dataset",
        "variables": variables,
        "dims": {str(k): int(v) for k, v in ds.sizes.items()},
        "temporal_coverage": temporal_coverage(ds),
        "spatial_coverage": spatial_coverage(ds),
        "crs": crs.splitlines()[0][:120] if isinstance(crs, str) else None,
    }


def summarize_output(obj: Any) -> Optional[Dict[str, Any]]:
    """Summarise whatever a processing step returned."""
    if obj is None:
        return None
    if isinstance(obj, xr.Dataset):
        return summarize_dataset(obj)
    if isinstance(obj, pd.DataFrame):
        return {"type": "pandas.DataFrame", "rows": int(len(obj)), "columns": list(map(str, obj.columns))}
    if isinstance(obj, (str, Path)) and Path(obj).exists():
        return {"type": "file", **file_record(obj)}
    return {"type": type(obj).__name__}


def software_environment() -> Dict[str, Any]:
    """Python, platform and versions of the packages that shape the output."""
    from importlib import metadata

    import climdata

    # The imported module, not the installed distribution: a source checkout on
    # PYTHONPATH can shadow an older installed climdata.
    packages = {"climdata": getattr(climdata, "__version__", "unknown")}
    for pkg in _TRACKED_PACKAGES[1:]:
        try:
            packages[pkg] = metadata.version(pkg)
        except metadata.PackageNotFoundError:
            continue
    return {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "packages": packages,
    }


def config_snapshot(cfg) -> Dict[str, Any]:
    """The run's configuration, minus the catalogues composed into it.

    ``cfg`` carries the full dataset, variable, index and region catalogues
    from ``conf/mappings``. Only the entries this run used are kept: the chosen
    dataset's description, the requested variables' metadata, the selected
    index definition and imputation method, and the selected region's bounds.
    """
    if cfg is None:
        return {}
    full = jsonable(cfg) or {}
    snap = {k: v for k, v in full.items() if k not in _CATALOGUE_KEYS}

    used: Dict[str, Any] = {}
    dataset = str(full.get("dataset") or "")
    dsinfo = (full.get("dsinfo") or {}).get(dataset.upper())
    if dsinfo:
        used["dataset"] = {k: v for k, v in dsinfo.items() if k not in ("variables", "params")}
    variables = full.get("variables") or []
    varinfo = full.get("varinfo") or {}
    if varinfo:
        used["variables"] = {v: varinfo[v] for v in variables if v in varinfo}
    index = full.get("index")
    indices = (full.get("extinfo") or {}).get("indices", {})
    if index and index in indices:
        used["index"] = {index: indices[index]}
    method = full.get("impute")
    if method and method in (full.get("imputeinfo") or {}):
        used["impute"] = {method: full["imputeinfo"][method]}
    region = full.get("region")
    if region and region in (full.get("bounds") or {}):
        used["region"] = {region: full["bounds"][region]}
    if used:
        snap["resolved"] = used
    return snap


def history_lines(steps: List[Dict[str, Any]], version: str) -> List[str]:
    """One CF-``history`` line per completed step, oldest first."""
    lines = []
    for s in steps:
        if s.get("status") != "completed":
            continue
        params = " ".join(
            f"{k}={json.dumps(v, separators=(',', ':')) if not isinstance(v, str) else v}"
            for k, v in (s.get("parameters") or {}).items()
            if not isinstance(v, (dict, list)) or len(json.dumps(v)) < 120
        )
        lines.append(f"{s.get('ended', '')} climdata {version} {s['action']} {params}".rstrip())
    return lines


# ---------------------------------------------------------------------------
# RO-Crate manifest
# ---------------------------------------------------------------------------
_ACTION_STATUS = {
    "completed": "http://schema.org/CompletedActionStatus",
    "failed": "http://schema.org/FailedActionStatus",
    "skipped": "http://schema.org/PotentialActionStatus",
}


def _is_url(s) -> bool:
    return isinstance(s, str) and s.startswith(("http://", "https://"))


def _license_ref(lic):
    """A URL licence is a linked entity; a licence name is plain text."""
    return {"@id": lic} if _is_url(lic) else lic


def _property_values(params: Dict[str, Any]) -> List[Dict[str, Any]]:
    out = []
    for k, v in (params or {}).items():
        value = v if isinstance(v, (str, int, float, bool)) or v is None else json.dumps(v)
        out.append({"@type": "PropertyValue", "name": k, "value": value})
    return out


def build_ro_crate(meta: Dict[str, Any]) -> Dict[str, Any]:
    """Render ``meta.json`` content as an RO-Crate 1.1 JSON-LD manifest.

    Each processing step becomes a ``CreateAction`` (Process Run Crate
    profile) whose ``instrument`` is climdata, with its parameters as
    ``PropertyValue`` entries. Input files of upload steps become ``object``;
    the files a step wrote into the crate become its ``result``.
    """
    data = meta["data"]
    data_files = [f for f in meta["files"] if f.get("role") == "data"]
    software = meta["software"]
    version = software["packages"].get("climdata", "unknown")
    tc, sc = data.get("temporal_coverage"), data.get("spatial_coverage")

    root: Dict[str, Any] = {
        "@id": "./",
        "@type": "Dataset",
        "identifier": meta["crate_id"],
        "name": meta["title"],
        "description": meta["description"],
        "datePublished": meta["created"],
        "license": _license_ref(meta["license"]),
        "keywords": [v["name"] for v in data["variables"]],
        "hasPart": [{"@id": f["path"]} for f in data_files] + [{"@id": "meta.json"}],
        "mentions": [{"@id": f"#step-{s['step']}"} for s in meta["processing"]],
    }
    if meta.get("creator"):
        root["creator"] = {"@id": "#creator"}
    if tc:
        root["temporalCoverage"] = f"{tc['start']}/{tc['end']}"
    if sc:
        root["spatialCoverage"] = {"@id": "#extent"}

    graph: List[Dict[str, Any]] = [
        {
            "@id": "ro-crate-metadata.json",
            "@type": "CreativeWork",
            "conformsTo": [
                {"@id": "https://w3id.org/ro/crate/1.1"},
                {"@id": "https://w3id.org/ro/wfrun/process/0.5"},
            ],
            "about": {"@id": "./"},
        },
        root,
        *[{
            "@id": f["path"],
            "@type": "File",
            "name": f["path"],
            "description": f.get("description"),
            "encodingFormat": f.get("encodingFormat"),
            "contentSize": str(f["size_bytes"]),
            "sha256": f["sha256"],
            "variableMeasured": [
                {"@type": "PropertyValue", "name": v["name"],
                 "description": v.get("long_name") or v["name"],
                 "unitText": v.get("units") or "1"}
                for v in data["variables"]
            ],
            **({"temporalCoverage": f"{tc['start']}/{tc['end']}"} if tc else {}),
        } for f in data_files],
        {
            "@id": "meta.json",
            "@type": "File",
            "name": "meta.json",
            "description": "Full processing record: configuration, every step "
                           "with parameters and timings, software environment, checksums.",
            "encodingFormat": "application/json",
        },
        {
            "@id": "#climdata",
            "@type": "SoftwareApplication",
            "name": "climdata",
            "version": version,
            "url": "https://github.com/Kaushikreddym/climdata",
        },
    ]
    if _is_url(meta["license"]):
        graph.append({"@id": meta["license"], "@type": "CreativeWork", "name": meta["license"]})
    if meta.get("creator"):
        graph.append({"@id": "#creator", "@type": "Person", "name": meta["creator"]})
    if sc:
        graph.append({
            "@id": "#extent",
            "@type": "Place",
            "geo": {"@type": "GeoShape",
                    "box": f"{sc['lat_min']} {sc['lon_min']} {sc['lat_max']} {sc['lon_max']}"},
        })

    for s in meta["processing"]:
        action: Dict[str, Any] = {
            "@id": f"#step-{s['step']}",
            "@type": "CreateAction",
            "name": f"{s['step']}. {s['action']}",
            "instrument": {"@id": "#climdata"},
            "startTime": s.get("started"),
            "endTime": s.get("ended"),
            "actionStatus": {"@id": _ACTION_STATUS.get(s.get("status"), _ACTION_STATUS["completed"])},
            "additionalProperty": _property_values(s.get("parameters")),
        }
        if s.get("error"):
            action["error"] = s["error"]
        inputs = [p for p in (s.get("inputs") or [])]
        if inputs:
            action["object"] = [{"@id": f["path"]} for f in inputs]
            for f in inputs:
                graph.append({"@id": f["path"], "@type": "File",
                              "contentSize": str(f.get("size_bytes", "")),
                              "sha256": f.get("sha256")})
        if s.get("results"):
            action["result"] = [{"@id": r} for r in s["results"]]
        graph.append(action)

    return {
        "@context": [
            "https://w3id.org/ro/crate/1.1/context",
            {"climdata": "https://climdata.readthedocs.io/", "sha256": "https://w3id.org/ro/terms/workflow-run#sha256"},
        ],
        "@graph": graph,
    }


# ---------------------------------------------------------------------------
# HTML preview
# ---------------------------------------------------------------------------
_CSS = """
:root{--fg:#1f2328;--muted:#59636e;--bg:#fff;--line:#d1d9e0;--head:#f6f8fa;--accent:#0969da;
  --ok:#1a7f37;--fail:#cf222e;--skip:#9a6700}
@media (prefers-color-scheme:dark){:root{--fg:#e6edf3;--muted:#9198a1;--bg:#0d1117;--line:#3d444d;
  --head:#151b23;--accent:#4493f8;--ok:#3fb950;--fail:#f85149;--skip:#d29922}}
body{font:15px/1.5 system-ui,-apple-system,"Segoe UI",sans-serif;color:var(--fg);background:var(--bg);
  max-width:1000px;margin:2rem auto;padding:0 16px}
h1{font-size:1.6rem;margin-bottom:.2rem}h2{font-size:1.15rem;margin-top:2rem;border-bottom:1px solid var(--line);padding-bottom:.3rem}
.sub{color:var(--muted);margin-top:0}
table{border-collapse:collapse;width:100%;font-size:.9rem;display:block;overflow-x:auto}
th,td{border:1px solid var(--line);padding:.35rem .6rem;text-align:left;vertical-align:top}
th{background:var(--head);font-weight:600}
code{font:12.5px ui-monospace,SFMono-Regular,Menlo,monospace;word-break:break-all}
.ok{color:var(--ok)}.failed{color:var(--fail)}.skipped{color:var(--skip)}
a{color:var(--accent)}ul.p{margin:0;padding-left:1rem}
details summary{cursor:pointer;color:var(--accent)}
pre{background:var(--head);border:1px solid var(--line);padding:.8rem;overflow-x:auto;font-size:12px}
"""


def _e(x) -> str:
    return html.escape("" if x is None else str(x))


def _fmt_size(n) -> str:
    n = float(n or 0)
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} GB"


def _fmt_output(out: Optional[Dict[str, Any]]) -> str:
    if not out:
        return ""
    if out.get("type") == "xarray.Dataset":
        vars_ = ", ".join(v["name"] for v in out.get("variables", []))
        dims = " × ".join(f"{k}={v}" for k, v in out.get("dims", {}).items())
        tc = out.get("temporal_coverage")
        period = f"<br>{_e(tc['start'])} → {_e(tc['end'])}" if tc else ""
        return f"<code>{_e(vars_)}</code><br>{_e(dims)}{period}"
    if out.get("type") == "pandas.DataFrame":
        return f"DataFrame, {out['rows']:,} rows"
    if out.get("type") == "file":
        if "files" in out:
            return f"<code>{_e(out['path'])}</code> ({len(out['files'])} files)"
        return f"<code>{_e(out['path'])}</code><br>{_fmt_size(out.get('size_bytes'))}"
    return _e(out.get("type"))


def render_preview(meta: Dict[str, Any]) -> str:
    """Render ``meta.json`` content as a standalone HTML summary page."""
    data, sw = meta["data"], meta["software"]
    tc, sc = data.get("temporal_coverage"), data.get("spatial_coverage")
    period = f"{tc['start']} → {tc['end']}" if tc else "—"
    extent = (f"{sc['lat_min']}–{sc['lat_max']} °N, {sc['lon_min']}–{sc['lon_max']} °E"
              if sc else "—")
    dims = ", ".join(f"{k}: {v}" for k, v in data.get("dims", {}).items())

    overview = [
        ("Source dataset", meta.get("source_dataset") or "—"),
        ("Period", period),
        ("Spatial extent", extent),
        ("Dimensions", dims or "—"),
        ("Created", meta["created"]),
        ("Updated", meta.get("updated") or meta["created"]),
        ("Creator", meta.get("creator") or "—"),
        ("License", f'<a href="{_e(meta["license"])}">{_e(meta["license"])}</a>'
                    if _is_url(meta["license"]) else _e(meta["license"])),
        ("Produced by", f"climdata {sw['packages'].get('climdata', '')}"),
        ("Crate ID", f"<code>{_e(meta['crate_id'])}</code>"),
    ]
    overview_rows = "".join(
        f"<tr><th>{_e(k)}</th><td>{v if k in ('License', 'Crate ID') else _e(v)}</td></tr>"
        for k, v in overview
    )

    var_rows = "".join(
        f"<tr><td><code>{_e(v['name'])}</code></td><td>{_e(v.get('long_name'))}</td>"
        f"<td>{_e(v.get('units'))}</td><td>{_e(', '.join(v.get('dims', [])))}</td>"
        f"<td>{_e(v.get('dtype'))}</td></tr>"
        for v in data["variables"]
    )

    step_rows = []
    for s in meta["processing"]:
        params = "".join(
            f"<li><code>{_e(k)}</code> = <code>{_e(json.dumps(v) if not isinstance(v, str) else v)}</code></li>"
            for k, v in (s.get("parameters") or {}).items()
        )
        inputs = "".join(
            f"<li>input <code>{_e(f['path'])}</code></li>" for f in (s.get("inputs") or [])
        ) + "".join(
            f"<li>wrote <a href=\"{_e(r)}\"><code>{_e(r)}</code></a></li>" for r in (s.get("results") or [])
        )
        status = s.get("status", "")
        status_html = f'<span class="{_e(status) if status != "completed" else "ok"}">{_e(status)}</span>'
        if s.get("error"):
            status_html += f"<br><code>{_e(s['error'])}</code>"
        step_rows.append(
            f"<tr><td>{s['step']}</td><td><b>{_e(s['action'])}</b><br>{status_html}</td>"
            f"<td>{_e(s.get('started'))}<br>{_e(s.get('duration_s'))} s</td>"
            f"<td><ul class='p'>{params}{inputs}</ul></td>"
            f"<td>{_fmt_output(s.get('output'))}</td></tr>"
        )

    file_rows = "".join(
        f"<tr><td><a href=\"{_e(f['path'])}\">{_e(f['path'])}</a></td><td>{_e(f.get('description'))}</td>"
        f"<td>{_fmt_size(f.get('size_bytes')) if f.get('size_bytes') else ''}</td>"
        f"<td><code>{_e(f.get('sha256', ''))}</code></td></tr>"
        for f in meta.get("files", [])
    )

    pkg_rows = "".join(
        f"<tr><td>{_e(k)}</td><td>{_e(v)}</td></tr>" for k, v in sw["packages"].items()
    )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{_e(meta['title'])}</title>
<style>{_CSS}</style>
</head>
<body>
<h1>{_e(meta['title'])}</h1>
<p class="sub">{_e(meta['description'])}</p>

<h2>Overview</h2>
<table>{overview_rows}</table>

<h2>Variables</h2>
<table><tr><th>Name</th><th>Long name</th><th>Units</th><th>Dimensions</th><th>Type</th></tr>{var_rows}</table>

<h2>Processing history</h2>
<p class="sub">Every step that produced the files in this crate, in the order it ran.</p>
<table><tr><th>#</th><th>Step</th><th>Started (UTC) / duration</th><th>Parameters</th><th>Output</th></tr>
{''.join(step_rows)}</table>

<h2>Files</h2>
<table><tr><th>File</th><th>Content</th><th>Size</th><th>SHA-256</th></tr>{file_rows}</table>

<h2>Software environment</h2>
<p class="sub">Python {_e(sw['python'])} on {_e(sw['platform'])}</p>
<table><tr><th>Package</th><th>Version</th></tr>{pkg_rows}</table>

<h2>Configuration</h2>
<details><summary>Show the configuration this run used</summary>
<pre>{_e(json.dumps(meta.get('config', {}), indent=2))}</pre></details>
</body>
</html>
"""
