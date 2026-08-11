"""Map visualisation helpers for climdata xarray datasets.

Provides :func:`plot_map`, a thin wrapper around xarray's ``.plot`` that draws
the field on Cartopy ``GeoAxes`` with coastlines and country boundaries.

Cartopy and Matplotlib are imported lazily so that ``import climdata`` does not
require them; a clear error is raised only if a plotting function is called
without them installed.
"""
from __future__ import annotations

import xarray as xr


def _require_plotting():
    """Import Matplotlib/Cartopy on demand, with a helpful error if missing."""
    try:
        import matplotlib.pyplot as plt
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
    except ImportError as e:  # pragma: no cover - depends on environment
        raise ImportError(
            "Plotting requires matplotlib and cartopy. "
            "Install with: conda install -c conda-forge cartopy matplotlib"
        ) from e
    return plt, ccrs, cfeature


def _as_dataarray(data, variable=None) -> xr.DataArray:
    """Return a DataArray from a DataArray or a (Dataset, variable) pair."""
    if isinstance(data, xr.Dataset):
        if variable is None:
            data_vars = list(data.data_vars)
            if len(data_vars) != 1:
                raise ValueError(
                    f"Dataset has {len(data_vars)} variables {data_vars}; "
                    "pass `variable=` to choose which to plot."
                )
            variable = data_vars[0]
        return data[variable]
    if isinstance(data, xr.DataArray):
        return data
    raise TypeError(f"Expected an xarray Dataset or DataArray, got {type(data)!r}.")


def _spatial_dim_names(da: xr.DataArray):
    """Find (x, y) coordinate names, preferring cf_xarray, then common aliases."""
    x = y = None
    try:
        x = da.cf["longitude"].name
        y = da.cf["latitude"].name
    except (KeyError, AttributeError):
        for cand in ("lon", "longitude", "x"):
            if cand in da.coords or cand in da.dims:
                x = cand
                break
        for cand in ("lat", "latitude", "y"):
            if cand in da.coords or cand in da.dims:
                y = cand
                break
    if x is None or y is None:
        raise ValueError(
            "Could not identify longitude/latitude coordinates on the data. "
            f"Available coords: {list(da.coords)}."
        )
    return x, y


def _reduce_to_2d(da: xr.DataArray, x, y, time, reduce):
    """Collapse non-spatial dims to obtain a single 2D (y, x) field."""
    spatial = {x, y}
    other_dims = [d for d in da.dims if d not in spatial]

    for dim in other_dims:
        if dim == "time":
            if reduce in ("mean", "sum", "min", "max", "median", "std"):
                da = getattr(da, reduce)(dim="time")
            elif isinstance(time, int):
                da = da.isel(time=time)
            else:
                da = da.sel(time=time, method="nearest")
        else:
            # Any remaining extra dim (e.g. geom_id): take the first slice.
            da = da.isel({dim: 0})
    return da


def plot_map(
    data,
    *,
    variable=None,
    time=0,
    reduce=None,
    ax=None,
    projection=None,
    figsize=(10, 5),
    coastlines=True,
    borders=True,
    gridlines=True,
    land=False,
    ocean=False,
    add_colorbar=True,
    robust=True,
    cmap=None,
    title=None,
    **plot_kwargs,
):
    """Plot a 2D geographic field with coastlines and country boundaries.

    Args:
        data: An ``xr.DataArray`` or ``xr.Dataset``. For a Dataset pass ``variable``
            (or it must contain exactly one data variable).
        variable (str, optional): Variable name when ``data`` is a Dataset.
        time (int | str | datetime, optional): Time slice to plot when a ``time``
            dimension is present. An ``int`` selects by position (``isel``);
            anything else selects by label (nearest). Ignored if ``reduce`` is set.
        reduce (str, optional): Reduce the time dimension instead of slicing it —
            one of ``"mean"``, ``"sum"``, ``"min"``, ``"max"``, ``"median"``, ``"std"``.
        ax: Existing Cartopy ``GeoAxes`` to draw on. If ``None``, one is created.
        projection: A ``cartopy.crs`` projection for a new axis (default ``PlateCarree``).
        figsize (tuple): Figure size when creating a new figure.
        coastlines, borders, gridlines, land, ocean (bool): Map features to add.
        add_colorbar (bool): Whether to draw a colorbar.
        robust (bool): Use 2nd/98th percentiles for color limits (passed to xarray).
        cmap (str, optional): Matplotlib colormap name.
        title (str, optional): Axis title. Defaults to the variable name + time.
        **plot_kwargs: Forwarded to ``xarray.DataArray.plot`` (e.g. ``vmin``, ``vmax``).

    Returns:
        The Cartopy ``GeoAxes`` the field was drawn on.
    """
    plt, ccrs, cfeature = _require_plotting()

    da = _as_dataarray(data, variable)
    x, y = _spatial_dim_names(da)
    da = _reduce_to_2d(da, x, y, time, reduce)

    if da.ndim != 2:
        raise ValueError(
            f"Expected a 2D field after reduction but got dims {da.dims}. "
            "This can happen for point/station extractions that have no spatial grid."
        )

    data_crs = ccrs.PlateCarree()
    if ax is None:
        if projection is None:
            projection = ccrs.PlateCarree()
        fig, ax = plt.subplots(figsize=figsize, subplot_kw={"projection": projection})

    da.plot(
        ax=ax,
        x=x,
        y=y,
        transform=data_crs,
        robust=robust,
        cmap=cmap,
        add_colorbar=add_colorbar,
        **plot_kwargs,
    )

    if land:
        ax.add_feature(cfeature.LAND, facecolor="none", zorder=0)
    if ocean:
        ax.add_feature(cfeature.OCEAN, zorder=0)
    if coastlines:
        ax.coastlines(resolution="50m", linewidth=0.6)
    if borders:
        ax.add_feature(cfeature.BORDERS, linewidth=0.4, edgecolor="gray")
    if gridlines:
        gl = ax.gridlines(draw_labels=True, linewidth=0.3, color="gray", alpha=0.5)
        gl.top_labels = False
        gl.right_labels = False

    if title is None:
        name = da.name or (variable if variable else "")
        title = str(name)
        if "time" in da.coords and da["time"].ndim == 0:
            try:
                import pandas as pd

                title = f"{title} @ {pd.to_datetime(da['time'].values):%Y-%m-%d}"
            except Exception:
                title = f"{title} @ {da['time'].values}"
        elif reduce:
            title = f"{title} ({reduce} over time)"
    ax.set_title(title)

    return ax
