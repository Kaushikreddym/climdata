"""Thin wrapper around ``climdata.ClimData`` for GUI consumption."""

from __future__ import annotations

from typing import Any, List, Optional, Sequence


def get_datasets() -> List[str]:
    """Return available dataset names from ClimData configuration.

    Falls back to an empty list if the climdata stack is not installed.
    """
    try:
        from climdata import ClimData
        return ClimData().get_datasets()
    except Exception:
        return []


def run_pipeline(overrides: Sequence[str], seq: Optional[Sequence[str]] = None) -> Any:
    """Run a climdata workflow with the given Hydra overrides.

    Args:
        overrides: Hydra override strings, e.g. ``["dataset=mswx", "lat=52.5"]``.
        seq: Ordered list of workflow actions. Defaults to
            ``["extract", "to_csv"]``.

    Returns:
        The ``WorkflowResult`` from ``ClimData.run_workflow``.

    Notes:
        Imports ``climdata`` lazily so the GUI module can be imported (and
        the scaffold inspected) even when the scientific stack is not yet
        installed.
    """
    from climdata import ClimData  # lazy import

    actions: List[str] = list(seq) if seq else ["extract", "to_csv"]
    extractor = ClimData(overrides=list(overrides))
    return extractor.run_workflow(actions=actions)


def _geo_bounds(da):
    """Return ``(lat_min, lat_max, lon_min, lon_max)`` in WGS84 degrees.

    Prefers 2-D ``lat``/``lon`` coordinate arrays (projected grids such as
    HYRAS/DWD that use LAEA x/y dims but attach geographic coords).
    Falls back to 1-D dimension coordinates named ``lat``/``lon``/``latitude``/
    ``longitude``.
    """
    import numpy as np

    # 2-D lat/lon coordinate variables (e.g. HYRAS where dims are y, x)
    lat_c = next((c for c in da.coords
                  if c in ("lat", "latitude") and len(da[c].dims) > 1), None)
    lon_c = next((c for c in da.coords
                  if c in ("lon", "longitude") and len(da[c].dims) > 1), None)
    if lat_c and lon_c:
        lat_arr = np.asarray(da[lat_c]).ravel()
        lon_arr = np.asarray(da[lon_c]).ravel()
        return (float(np.nanmin(lat_arr)), float(np.nanmax(lat_arr)),
                float(np.nanmin(lon_arr)), float(np.nanmax(lon_arr)))

    # 1-D dimension coordinates
    lat_d = next((d for d in da.dims if d in ("lat", "latitude")), None)
    lon_d = next((d for d in da.dims if d in ("lon", "longitude")), None)
    if lat_d and lon_d:
        lats = da[lat_d].values.astype(float)
        lons = da[lon_d].values.astype(float)
        return (float(lats.min()), float(lats.max()),
                float(lons.min()), float(lons.max()))

    return None


def get_overlay_meta(result):
    try:
        import numpy as np
        ds = result.dataset
        print(f"[DEBUG get_overlay_meta] dataset={ds}", flush=True)
        if ds is None:
            print("[DEBUG get_overlay_meta] dataset is None", flush=True)
            return None
        vars_list = [v for v in ds.data_vars if not v.startswith("quality_")]
        print(f"[DEBUG get_overlay_meta] vars={vars_list}", flush=True)
        if not vars_list:
            return None
        da0 = ds[vars_list[0]]
        print(f"[DEBUG get_overlay_meta] da0.dims={da0.dims}  da0.coords={list(da0.coords)}", flush=True)
        bounds = _geo_bounds(da0)
        print(f"[DEBUG get_overlay_meta] WGS84 bounds={bounds}", flush=True)
        if bounds is None:
            return None
        lat_min, lat_max, lon_min, lon_max = bounds
        meta = {
            "vars": vars_list,
            "lat_min": lat_min, "lat_max": lat_max,
            "lon_min": lon_min, "lon_max": lon_max,
        }
        print(f"[DEBUG get_overlay_meta] meta={meta}", flush=True)
        return meta
    except Exception as exc:
        print(f"[DEBUG get_overlay_meta] EXCEPTION: {exc}", flush=True)
        return None


def render_variable(ds, var_name: str, colormap: str = "RdYlBu_r",
                    vmin: Optional[float] = None,
                    vmax: Optional[float] = None) -> Optional[dict]:
    """Render the time-mean of *var_name* as a transparent PNG + colorbar.

    Args:
        vmin: Colour-scale minimum.  ``None`` → auto (data minimum).
        vmax: Colour-scale maximum.  ``None`` → auto (data maximum).

    Returns ``{"b64", "colorbar_b64", "vmin", "vmax", "lat_min", "lon_min",
    "lat_max", "lon_max"}`` or ``None`` on failure.
    """
    import base64, io, pathlib, tempfile

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        import numpy as np
    except ImportError:
        return None

    try:
        if var_name not in ds:
            return None
        da = ds[var_name]
        lat_dim = next((d for d in da.dims if d in ("lat", "latitude", "y")), None)
        lon_dim = next((d for d in da.dims if d in ("lon", "longitude", "x")), None)
        if lat_dim is None or lon_dim is None:
            return None

        mean = da.mean(dim="time") if "time" in da.dims else da
        # projected dim coords (used only for image orientation)
        lats_proj = mean[lat_dim].values.astype(float)
        lons_proj = mean[lon_dim].values.astype(float)

        vals = mean.values
        if mean.dims[0] != lat_dim:
            vals = vals.T
        if lats_proj[0] < lats_proj[-1]:   # y increases northward → flip
            vals = np.flipud(vals)

        auto_vmin = float(np.nanmin(vals))
        auto_vmax = float(np.nanmax(vals))
        plot_vmin = vmin if vmin is not None else auto_vmin
        plot_vmax = vmax if vmax is not None else auto_vmax
        units = da.attrs.get("units", "")

        # WGS84 bounds for Leaflet overlay
        bounds = _geo_bounds(da)
        if bounds is None:
            return None
        lat_min, lat_max, lon_min, lon_max = bounds
        print(f"[DEBUG render_variable] '{var_name}' WGS84 bounds: lat=[{lat_min},{lat_max}] lon=[{lon_min},{lon_max}]", flush=True)

        # ---- raster image (axes-less, transparent background) ----
        fig, ax = plt.subplots(figsize=(6, 5), dpi=100)
        ax.axis("off")
        ax.imshow(vals, cmap=colormap, aspect="auto",
                  extent=[lon_min, lon_max, lat_min, lat_max],
                  vmin=plot_vmin, vmax=plot_vmax)
        buf = io.BytesIO()
        fig.savefig(buf, format="png", transparent=True, bbox_inches="tight", dpi=100)
        plt.close(fig)
        buf.seek(0)
        raw = buf.read()

        # DEBUG: save to /tmp for inspection
        dbg = pathlib.Path(tempfile.gettempdir()) / f"climdata_overlay_{var_name}_debug.png"
        dbg.write_bytes(raw)
        print(f"[DEBUG] Overlay PNG for '{var_name}' ({colormap}) → {dbg}")

        b64 = base64.b64encode(raw).decode()

        # ---- colorbar image ----
        cbar_b64 = _render_colorbar(colormap, plot_vmin, plot_vmax, var_name, units)

        return {
            "b64": b64,
            "colorbar_b64": cbar_b64,
            "vmin": plot_vmin, "vmax": plot_vmax,
            "lat_min": lat_min, "lon_min": lon_min,
            "lat_max": lat_max, "lon_max": lon_max,
        }
    except Exception:
        return None


def _render_colorbar(cmap_name: str, vmin: float, vmax: float,
                     label: str, units: str) -> str:
    """Render a horizontal colorbar strip as a base64 PNG."""
    import base64, io
    try:
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(3.6, 0.38), dpi=90)
        fig.patch.set_alpha(0)
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap_name, norm=norm,
                                       orientation="horizontal")
        lbl = f"{label} [{units}]" if units else label
        cb.set_label(lbl, fontsize=7)
        cb.ax.tick_params(labelsize=6)
        fig.tight_layout(pad=0.2)
        buf = io.BytesIO()
        fig.savefig(buf, format="png", transparent=True, bbox_inches="tight", dpi=90)
        plt.close(fig)
        buf.seek(0)
        return base64.b64encode(buf.read()).decode()
    except Exception:
        return ""
