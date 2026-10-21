"""Bias correction and statistical downscaling, wrapping ISIMIP3BASD for xarray.

Two stages, usable separately or chained by :class:`BCSD`:

* :class:`BiasCorrection` removes systematic bias from a model simulation on the
  *coarse* grid, by trend-preserving quantile mapping against observations
  regridded to that grid.
* :class:`StatisticalDownscaling` maps the corrected coarse field onto the fine
  observational grid with a modified MBCn algorithm.

The order is deliberate: correcting before downscaling keeps the expensive
quantile mapping on the smaller grid, and gives the downscaling step a field
whose distribution already matches the observations.

The vendored ISIMIP3BASD code works in ``iris`` cubes on a proleptic Gregorian
calendar, so this module also carries the conversions that get real GCM output
into that shape — notably :func:`_to_proleptic_gregorian`, which inserts a 29
February into ``365_day`` model calendars rather than relabelling and silently
shifting every date.

Requires ``iris``; ``xesmf`` or a CDO installation is needed for the regridding
step.

Reference:
    Lange, S. (2019). Trend-preserving bias adjustment and statistical
    downscaling with ISIMIP3BASD (v1.0). Geoscientific Model Development,
    12(7), 3055-3070. https://doi.org/10.5194/gmd-12-3055-2019

Example:
    >>> from climdata.sdba import BCSD                             # doctest: +SKIP
    >>> bcsd = BCSD(variable="tas")                                # doctest: +SKIP
    >>> paths = bcsd.run(obs_fine, sim_hist_coarse, sim_fut_coarse)  # doctest: +SKIP
"""

import xarray as xr
import numpy as np
from pathlib import Path
from typing import Optional, Union, Dict, List
import shutil
import warnings
import subprocess
import os

try:
    import xesmf as xe
    XESMF_AVAILABLE = True
except ImportError:
    XESMF_AVAILABLE = False
    xe = None

try:
    import iris
    IRIS_AVAILABLE = True
except ImportError:
    IRIS_AVAILABLE = False
    iris = None

warnings.filterwarnings("ignore", category=Warning)


# ---------------------------------------------------------------------------
# Iris / calendar helpers
# ---------------------------------------------------------------------------

def _to_proleptic_gregorian(cube):
    """Convert an iris cube time coord to proleptic_gregorian (ISIMIP3BASD requirement).

    The point values count time steps *in the source calendar*, so a 365_day
    (noleap) coord cannot simply be relabelled: every date would slip backwards by
    the number of leap days between the unit origin and the data (40 days for
    GFDL-ESM4 / CanESM5 on 'days since 1850-01-01').  Convert points → dates →
    points instead, so each step keeps its real calendar date.  For
    standard/gregorian/julian a relabel is exact over the CF data range.

    365_day (noleap) GCMs never have a 29 Feb, which would otherwise leave a
    day-of-year gap in every real-calendar leap year and misalign day-of-year
    based steps (e.g. halfwin_upper_bound_climatology) against obs data.  See
    `_insert_feb29`, which duplicates 28 Feb as 29 Feb for those years.
    """
    import cf_units
    from datetime import datetime as _dt

    t = cube.coord('time')
    src = t.units
    if src.calendar == 'proleptic_gregorian':
        return cube

    tgt = cf_units.Unit(src.origin, calendar='proleptic_gregorian')

    if src.calendar in ('standard', 'gregorian', 'julian'):
        t.units = tgt
        return cube

    if src.calendar in ('360_day', '360'):
        raise NotImplementedError(
            f"calendar '{src.calendar}' has dates (e.g. Feb 30) with no "
            "proleptic_gregorian equivalent; resample the input to a real "
            "calendar before bias correction"
        )

    if src.calendar in ('365_day', 'noleap'):
        return _insert_feb29(cube, src, tgt)

    def _renumber(points):
        dates = np.atleast_1d(src.num2date(points))
        return tgt.date2num(
            [_dt(d.year, d.month, d.day, d.hour, d.minute, d.second) for d in dates]
        )

    new_points = _renumber(t.points)
    new_bounds = (
        _renumber(t.bounds.ravel()).reshape(t.bounds.shape)
        if t.bounds is not None else None
    )
    t.units = tgt
    t.points = new_points
    t.bounds = new_bounds
    return cube


def _insert_feb29(cube, src, tgt):
    """Duplicate 28 Feb as 29 Feb in every Gregorian leap year, for a 365_day source.

    The source calendar has no 29 Feb, so a new time step has to be inserted
    (not just relabelled) to get a gapless proleptic_gregorian axis. We reuse
    28 Feb's data for the inserted 29 Feb step, then rebuild the cube since its
    time dimension length changes. Time bounds are dropped rather than carried
    over incorrectly, since downstream BASD code only reads time points.
    """
    import calendar as _calendar
    from iris.coords import DimCoord
    from iris.cube import Cube
    from datetime import datetime as _dt

    t = cube.coord('time')
    dates = np.atleast_1d(src.num2date(t.points))

    gather_idx = []
    new_dates = []
    for i, d in enumerate(dates):
        gather_idx.append(i)
        new_dates.append(_dt(d.year, d.month, d.day, d.hour, d.minute, d.second))
        if d.month == 2 and d.day == 28 and _calendar.isleap(d.year):
            gather_idx.append(i)
            new_dates.append(_dt(d.year, 2, 29, d.hour, d.minute, d.second))

    new_time_coord = DimCoord(
        tgt.date2num(new_dates), standard_name='time', units=tgt
    )

    if not any(g != i for i, g in enumerate(gather_idx)):
        # no leap years in range: nothing to insert, just relabel
        t.units = tgt
        t.points = new_time_coord.points
        return cube

    time_dim = cube.coord_dims(t)[0]
    new_data = np.take(cube.core_data(), gather_idx, axis=time_dim)

    dim_coords_and_dims = [
        (new_time_coord if c is t else c.copy(), cube.coord_dims(c)[0])
        for c in cube.coords(dim_coords=True)
    ]
    aux_coords_and_dims = [
        (c.copy(), cube.coord_dims(c)) for c in cube.coords(dim_coords=False) if c is not t
    ]

    return Cube(
        new_data,
        standard_name=cube.standard_name,
        long_name=cube.long_name,
        var_name=cube.var_name,
        units=cube.units,
        attributes=cube.attributes,
        cell_methods=cube.cell_methods,
        dim_coords_and_dims=dim_coords_and_dims,
        aux_coords_and_dims=aux_coords_and_dims,
    )


def _stamp_geogcs(cube):
    """
    Attach a GeogCS and set units='degrees' on horizontal DimCoords.

    iris.analysis.Linear requires both a CRS *and* units == 'degrees' (not
    'degrees_north' / 'degrees_east') on every lat/lon DimCoord.
    """
    import iris.coord_systems as ics
    import cf_units
    geog_cs = ics.GeogCS(6371229.0)
    deg = cf_units.Unit('degrees')
    for coord in cube.coords(dim_coords=True):
        if iris.util.guess_coord_axis(coord) in ('X', 'Y'):
            coord.coord_system = geog_cs
            coord.units = deg
    return cube


def _ensure_rectilinear_iris_cube(cube):
    """
    Guarantee that the cube has 1-D latitude/longitude DimCoords.

    ISIMIP3BASD vendor code requires:
    - 1-D DimCoords (not 2-D AuxCoords)
    - standard_name='latitude' / 'longitude'
    - units='degrees'
    Idempotent — safe to call multiple times.
    """
    import cf_units
    from iris.coords import DimCoord

    x_dim = cube.coords(axis="x", dim_coords=True)
    y_dim = cube.coords(axis="y", dim_coords=True)
    if (len(x_dim) == 1 and len(y_dim) == 1 and
            x_dim[0].ndim == 1 and y_dim[0].ndim == 1 and
            "latitude"  in str(y_dim[0].standard_name or "").lower() and
            "longitude" in str(x_dim[0].standard_name or "").lower()):
        for c in (x_dim[0], y_dim[0]):
            c.units = cf_units.Unit('degrees')
        return cube

    lat_c = lon_c = None
    for c in cube.coords():
        name = (str(c.standard_name or "") + str(c.name() or "")).lower()
        if "lat" in name and lat_c is None:
            lat_c = c
        if "lon" in name and lon_c is None:
            lon_c = c
    if lat_c is None or lon_c is None:
        return cube

    lat1d = lat_c.points if lat_c.ndim == 1 else np.mean(lat_c.points, axis=1)
    lon1d = lon_c.points if lon_c.ndim == 1 else np.mean(lon_c.points, axis=0)

    nd = cube.ndim
    for c in (lat_c, lon_c):
        try:
            cube.remove_coord(c)
        except Exception:
            pass
    cube.add_dim_coord(DimCoord(lat1d, standard_name="latitude",  units='degrees'), nd - 2)
    cube.add_dim_coord(DimCoord(lon1d, standard_name="longitude", units='degrees'), nd - 1)
    return cube


def _xda_to_cube(xda):
    """xarray DataArray → dask-backed proleptic_gregorian iris cube."""
    import dask.array as da
    cube = _to_proleptic_gregorian(xda.to_iris())
    cube.data = da.ma.masked_invalid(cube.core_data())
    return cube


# ---------------------------------------------------------------------------
# Regrid fine → coarse  (used to build obs_hist at GCM resolution)
# ---------------------------------------------------------------------------

def regrid_to_coarse(
    fine_data: xr.Dataset,
    coarse_template: xr.Dataset,
    method: str = 'conservative',
    regridding_tool: str = 'xesmf',
    cdo_method: str = 'remapcon',
    cdo_env: str = 'cdo_stable',
    weights_dir: Optional[str] = None
) -> xr.Dataset:
    """Regrid fine-resolution data to coarse resolution.

    Args:
        fine_data (xr.Dataset): Fine-resolution dataset to regrid.
        coarse_template (xr.Dataset): Coarse-resolution dataset used as the target grid.
        method (str): xESMF method: 'conservative' (default), 'bilinear', or 'nearest'.
        regridding_tool (str): 'xesmf' (default) or 'cdo'.
        cdo_method (str): CDO operator string when regridding_tool='cdo'.
        cdo_env (str): Conda environment that has CDO installed.
        weights_dir (str, optional): Directory to cache regridding weight files.

    Returns:
        xr.Dataset: Regridded dataset on the coarse grid.
    """
    print(f"Regridding fine→coarse using {regridding_tool}...")

    if regridding_tool == 'xesmf':
        if not XESMF_AVAILABLE:
            raise ImportError("xESMF is not installed. Use: pip install xesmf")

        weight_file = None
        if weights_dir:
            os.makedirs(weights_dir, exist_ok=True)
            weight_file = os.path.join(weights_dir, f'weights_fine_to_coarse_{method}.nc')

        regridder = xe.Regridder(
            fine_data, coarse_template, method=method,
            periodic=False, ignore_degenerate=True,
            filename=weight_file,
            reuse_weights=bool(weight_file and os.path.exists(weight_file))
        )

        # Build a valid-source mask so border cells with no valid source → NaN
        first_var = next(iter(fine_data.data_vars))
        fine_valid      = xr.where(fine_data[first_var].isel(time=0).notnull(), 1.0, 0.0)
        coarse_valid_frac = regridder(fine_valid)
        valid_mask      = coarse_valid_frac > 0

        result_vars = {}
        for var in fine_data.data_vars:
            regridded = regridder(fine_data[var].where(fine_data[var].notnull(), other=0.0))
            regridded = (regridded / coarse_valid_frac).where(valid_mask)
            regridded.attrs = fine_data[var].attrs
            result_vars[var] = regridded

        result = xr.Dataset(result_vars, attrs=fine_data.attrs)

    elif regridding_tool == 'cdo':
        import tempfile
        with tempfile.TemporaryDirectory() as tmpdir:
            fine_file     = os.path.join(tmpdir, 'fine_input.nc')
            template_file = os.path.join(tmpdir, 'coarse_template.nc')
            output_file   = os.path.join(tmpdir, 'regridded_output.nc')
            fine_data.to_netcdf(fine_file)
            coarse_template.to_netcdf(template_file)
            cmd = (f"conda run -n {cdo_env} "
                   f"cdo {cdo_method},{template_file} {fine_file} {output_file}")
            print(f"   Running: {cmd}")
            try:
                subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
                result = xr.open_dataset(output_file).load()
            except subprocess.CalledProcessError as e:
                raise RuntimeError(f"CDO regridding failed: {e.stderr}") from e
    else:
        raise ValueError(f"Unknown regridding_tool: {regridding_tool!r}. Use 'xesmf' or 'cdo'.")

    print("   Regridding complete.")
    return result


# ---------------------------------------------------------------------------
# Bias Correction
# ---------------------------------------------------------------------------

class BiasCorrection:
    """Bias correction using ISIMIP3BASD trend-preserving quantile mapping.

    Applies bias correction at the coarse GCM resolution before any spatial
    downscaling.  Supports single-variable and multivariate (MBCn) correction.

    Args:
        variable (str or list of str): Climate variable name(s). When a list is given
            MBCn is used.
        method (str): 'parametric' (default).
        distribution (str, optional): Parametric distribution override. Defaults come
            from DEFAULT_CONFIGS.
        trend_preservation (str, optional): Override for trend preservation ('additive',
            'multiplicative', 'mixed', 'bounded').
        detrend (bool): Detrend before correction.
        adjust_p_values (bool): Adjust p-values for perfect match in the reference
            period.
        lower_bound, lower_threshold, upper_bound, upper_threshold (float, optional):
            Physical bounds / thresholds.
        halfwin_upper_bound_climatology (int): Half-window for upper bound climatology
            (0 = disabled).
        n_processes (int, optional): Parallel worker count. Defaults to all CPUs.
        n_quantiles (int): Quantiles for non-parametric QM (default 50).
        n_iterations (int): MBCn iterations (default 20; 0 = univariate, no MBCn).
        randomization_seed (int, optional): RNG seed for reproducibility.

    Example:
        >>> bc = BiasCorrection(variable='tas', distribution='normal',
        ...                     trend_preservation='additive', detrend=True)
        >>> tas_corrected = bc.correct(obs_hist=obs_ds, sim_hist=gcm_hist_ds,
        ...                            sim_fut=gcm_fut_ds)
    """

    DEFAULT_CONFIGS = {
        'tas': {
            'distribution': 'normal', 'trend_preservation': 'additive',
            'detrend': True, 'adjust_p_values': False,
        },
        'tasmax': {
            'distribution': 'normal', 'trend_preservation': 'additive',
            'detrend': True, 'adjust_p_values': False,
        },
        'tasmin': {
            'distribution': 'normal', 'trend_preservation': 'additive',
            'detrend': True, 'adjust_p_values': False,
        },
        'pr': {
            'distribution': 'gamma', 'trend_preservation': 'mixed',
            'lower_bound': 0, 'lower_threshold': 0.1, 'adjust_p_values': True,
        },
        'rsds': {
            'distribution': 'beta', 'trend_preservation': 'bounded',
            'lower_bound': 0, 'lower_threshold': 0.01,
            'upper_bound': 1, 'upper_threshold': 0.9999,
            'adjust_p_values': True, 'halfwin_upper_bound_climatology': 15,
        },
        'hurs': {
            'distribution': 'beta', 'trend_preservation': 'bounded',
            'lower_bound': 0, 'lower_threshold': 0.01,
            'upper_bound': 100, 'upper_threshold': 99.99,
            'adjust_p_values': True,
        },
        'sfcWind': {
            'distribution': 'weibull', 'trend_preservation': 'mixed',
            'lower_bound': 0, 'lower_threshold': 0.01, 'adjust_p_values': True,
        },
        'psl': {
            'distribution': 'normal', 'trend_preservation': 'additive',
            'adjust_p_values': True, 'detrend': True,
        },
        'rlds': {
            'distribution': 'normal', 'trend_preservation': 'additive',
            'adjust_p_values': True, 'detrend': True,
        },
    }

    def __init__(
        self,
        variable: Union[str, List[str]],
        method: str = 'parametric',
        distribution: Optional[str] = None,
        trend_preservation: Optional[str] = None,
        detrend: bool = False,
        adjust_p_values: bool = False,
        lower_bound: Optional[float] = None,
        lower_threshold: Optional[float] = None,
        upper_bound: Optional[float] = None,
        upper_threshold: Optional[float] = None,
        halfwin_upper_bound_climatology: int = 0,
        n_processes: Optional[int] = None,
        n_quantiles: int = 50,
        n_iterations: int = 20,
        randomization_seed: Optional[int] = None,
        **kwargs
    ):
        """Configure a bias-correction run.

        Every argument is documented on the class. Where an argument is left as
        ``None``, the per-variable entry in :attr:`DEFAULT_CONFIGS` supplies the
        value, so a bare ``BiasCorrection("pr")`` already uses the distribution,
        bounds and trend-preservation mode appropriate to precipitation.

        Raises:
            ImportError: If ``iris`` is not installed.
        """
        self.variables = [variable] if isinstance(variable, str) else list(variable)
        self.variable  = self.variables[0] if len(self.variables) == 1 else None
        self.method    = method
        self.n_processes       = n_processes if n_processes is not None else os.cpu_count()
        self.n_quantiles       = n_quantiles
        self.n_iterations      = n_iterations
        self.randomization_seed = randomization_seed
        self.kwargs    = kwargs

        self._var_configs = []
        for v in self.variables:
            cfg = self.DEFAULT_CONFIGS.get(v, {})
            self._var_configs.append({
                'distribution':                  distribution       if distribution       is not None  else cfg.get('distribution'),
                'trend_preservation':            trend_preservation if trend_preservation is not None  else cfg.get('trend_preservation', 'additive'),
                'detrend':                       detrend            if detrend            is not False else cfg.get('detrend', False),
                'adjust_p_values':               adjust_p_values    if adjust_p_values    is not False else cfg.get('adjust_p_values', False),
                'lower_bound':                   lower_bound        if lower_bound        is not None  else cfg.get('lower_bound'),
                'lower_threshold':               lower_threshold    if lower_threshold    is not None  else cfg.get('lower_threshold'),
                'upper_bound':                   upper_bound        if upper_bound        is not None  else cfg.get('upper_bound'),
                'upper_threshold':               upper_threshold    if upper_threshold    is not None  else cfg.get('upper_threshold'),
                'halfwin_upper_bound_climatology': halfwin_upper_bound_climatology if halfwin_upper_bound_climatology != 0
                                                   else cfg.get('halfwin_upper_bound_climatology', 0),
            })

        # Scalar attrs for single-variable backward compatibility
        if len(self.variables) == 1:
            c = self._var_configs[0]
            self.distribution                    = c['distribution']
            self.trend_preservation              = c['trend_preservation']
            self.detrend                         = c['detrend']
            self.adjust_p_values                 = c['adjust_p_values']
            self.lower_bound                     = c['lower_bound']
            self.lower_threshold                 = c['lower_threshold']
            self.upper_bound                     = c['upper_bound']
            self.upper_threshold                 = c['upper_threshold']
            self.halfwin_upper_bound_climatology = c['halfwin_upper_bound_climatology']

        mv_label = ('(MBCn multivariate)' if self.n_iterations > 0 and len(self.variables) > 1
                    else '(MBCn univariate)' if self.n_iterations > 0
                    else '(univariate, no MBCn)')
        print(f"BiasCorrection initialized for {self.variables}")
        print(f"  n_iterations={self.n_iterations} {mv_label}  n_processes={self.n_processes}  seed={self.randomization_seed}")
        for v, c in zip(self.variables, self._var_configs):
            print(f"  [{v}] dist={c['distribution']}  trend={c['trend_preservation']}  "
                  f"detrend={c['detrend']}  adjust_p={c['adjust_p_values']}")

    # ------------------------------------------------------------------ helpers

    def _to_iris_cubes(self, obs_hist, sim_hist, sim_fut):
        """Convert three xarray Datasets to per-variable iris cube lists."""
        obs_cubes, sh_cubes, sf_cubes = [], [], []
        for var in self.variables:
            for ds, out in [(obs_hist, obs_cubes), (sim_hist, sh_cubes), (sim_fut, sf_cubes)]:
                out.append(_xda_to_cube(ds[var]))
        return obs_cubes, sh_cubes, sf_cubes

    def _ba_vendor_kwargs(self) -> dict:
        """Build per-variable kwarg lists consumed by ba.adjust_bias."""
        return {
            'distribution':                    [c['distribution']                   for c in self._var_configs],
            'trend_preservation':              [c['trend_preservation']              for c in self._var_configs],
            'detrend':                         [c['detrend']                         for c in self._var_configs],
            'adjust_p_values':                 [c['adjust_p_values']                 for c in self._var_configs],
            'halfwin_upper_bound_climatology': [c['halfwin_upper_bound_climatology'] for c in self._var_configs],
            'lower_bound':                     [c['lower_bound']                     for c in self._var_configs],
            'lower_threshold':                 [c['lower_threshold']                 for c in self._var_configs],
            'upper_bound':                     [c['upper_bound']                     for c in self._var_configs],
            'upper_threshold':                 [c['upper_threshold']                 for c in self._var_configs],
            'if_all_invalid_use':              [np.nan  for _ in self._var_configs],
            'unconditional_ccs_transfer':      [False   for _ in self._var_configs],
            'trendless_bound_frequency':       [False   for _ in self._var_configs],
        }

    # --------------------------------------------------- low-level hooks for
    # custom per-location loops (e.g. nexgddp_hyras_bc_zarr.py)

    def _prepare_ba_state(self, obs_cubes, sh_cubes, sf_cubes, step_size=0):
        """
        Set up global vendor BA state needed by adjust_bias_one_location.

        Returns (lazy_data, month_numbers, years, doys, space_shape, rotation_matrices).
        Call ba.initializer internally, then iterate locations with correct_one_location.
        """
        from climdata._vendor.isimip3basd import bias_adjustment as ba
        from climdata._vendor.isimip3basd import utility_functions as uf

        detrend = [c['detrend']                        for c in self._var_configs]
        halfwin = [c['halfwin_upper_bound_climatology'] for c in self._var_configs]

        cubes_map   = {'obs_hist': obs_cubes, 'sim_hist': sh_cubes, 'sim_fut': sf_cubes}
        lazy_data   = {}
        space_shape = None
        month_numbers = {k: None for k in cubes_map}
        years         = {k: None for k in cubes_map}
        doys          = {k: None for k in cubes_map}

        for key, cube_list in cubes_map.items():
            lazy_data[key] = [c.core_data() for c in cube_list]
            for i, cube in enumerate(cube_list):
                if space_shape is None:
                    space_shape = cube.shape[1:]
                uf.assert_calendar(cube, 'proleptic_gregorian')
                uf.assert_coord_axis(cube, 'time', 0)
                datetimes = cube.coord('time').units.num2date(cube.coord('time').points)
                if not step_size:
                    j = uf.convert_datetimes(datetimes, 'month_number')
                    if month_numbers[key] is None:
                        month_numbers[key] = j
                if step_size or detrend[i]:
                    j = uf.convert_datetimes(datetimes, 'year')
                    if years[key] is None:
                        years[key] = j
                if step_size or halfwin[i]:
                    j = uf.convert_datetimes(datetimes, 'day_of_year')
                    if doys[key] is None:
                        doys[key] = j
            if step_size:
                uf.assert_full_period_coverage(years[key], doys[key], key)

        if self.randomization_seed is not None:
            np.random.seed(self.randomization_seed)
        rotation_matrices = [
            uf.generateCREmatrix(len(self.variables))
            for _ in range(self.n_iterations)
        ]

        ba.initializer(lazy_data, month_numbers, years, doys)
        return lazy_data, month_numbers, years, doys, space_shape, rotation_matrices

    def correct_one_location(self, i_loc, space_shape, npy_paths, rotation_matrices=None, step_size=0):
        """
        Bias-correct one grid cell via ba.adjust_bias_one_location.

        Requires _prepare_ba_state to have been called first (sets vendor globals).
        """
        from climdata._vendor.isimip3basd import bias_adjustment as ba

        vkw = self._ba_vendor_kwargs()
        ba.adjust_bias_one_location(
            i_loc=i_loc,
            sim_fut_ba_path=[str(p) for p in npy_paths],
            space_shape=space_shape,
            halfwin_upper_bound_climatology=vkw.pop('halfwin_upper_bound_climatology'),
            lower_bound=vkw.pop('lower_bound'),
            lower_threshold=vkw.pop('lower_threshold'),
            upper_bound=vkw.pop('upper_bound'),
            upper_threshold=vkw.pop('upper_threshold'),
            if_all_invalid_use=vkw.pop('if_all_invalid_use'),
            step_size=step_size,
            rotation_matrices=rotation_matrices or [],
            randomization_seed=self.randomization_seed,
            n_quantiles=self.n_quantiles,
            **vkw,
        )

    # ------------------------------------------------------------------ public

    def correct(
        self,
        obs_hist: xr.Dataset,
        sim_hist: xr.Dataset,
        sim_fut: xr.Dataset,
        output_path: Optional[Union[str, Path]] = None,
        **kwargs,
    ) -> xr.Dataset:
        """Apply bias correction to climate model data.

        Args:
            obs_hist (xr.Dataset): Historical observations (e.g. 1980–2014).
            sim_hist (xr.Dataset): Historical model simulations matching the obs_hist
                period.
            sim_fut (xr.Dataset): Future model simulations to be corrected.
            output_path (str or Path, optional): Save the merged result to this NetCDF
                path.

        Returns:
            xr.Dataset: Bias-corrected future simulations.
        """
        if not IRIS_AVAILABLE:
            raise ImportError("iris is required: pip install scitools-iris")
        try:
            from climdata._vendor.isimip3basd import bias_adjustment as ba
            from climdata._vendor.isimip3basd import utility_functions as uf
        except ImportError:
            raise ImportError(
                "ISIMIP3BASD vendor code not found in climdata/_vendor/isimip3basd/."
            )

        mv_tag = (f"{len(self.variables)} variables (MBCn multivariate)"
                  if len(self.variables) > 1 else self.variables[0])
        print(f"\nBias correction: {mv_tag}")
        print(f"  obs_hist : {obs_hist.time.values[0]} → {obs_hist.time.values[-1]}")
        print(f"  sim_hist : {sim_hist.time.values[0]} → {sim_hist.time.values[-1]}")
        print(f"  sim_fut  : {sim_fut.time.values[0]} → {sim_fut.time.values[-1]}")

        tmpdir = Path.cwd() / "tmp_bcsd"
        tmpdir.mkdir(parents=True, exist_ok=True)
        try:
            print("  Converting to iris cubes...")
            obs_cubes, sh_cubes, sf_cubes = self._to_iris_cubes(obs_hist, sim_hist, sim_fut)

            # Resolve one output path per variable
            if output_path is not None and len(self.variables) == 1:
                ba_paths = [Path(output_path).resolve()]
            else:
                ba_paths = [
                    (Path(output_path).parent / f"{Path(output_path).stem}_{v}{Path(output_path).suffix}").resolve()
                    if output_path else (tmpdir / f"sim_fut_ba_{v}.nc").resolve()
                    for v in self.variables
                ]

            vkw = self._ba_vendor_kwargs()
            print(f"  Running ba.adjust_bias (n_iterations={self.n_iterations}) ...")
            ba.adjust_bias(
                obs_hist=obs_cubes,
                sim_hist=sh_cubes,
                sim_fut=sf_cubes,
                sim_fut_ba_path=[str(p) for p in ba_paths],
                n_processes=self.n_processes,
                n_quantiles=self.n_quantiles,
                n_iterations=self.n_iterations,
                randomization_seed=self.randomization_seed,
                detrend=vkw.pop('detrend'),
                halfwin_upper_bound_climatology=vkw.pop('halfwin_upper_bound_climatology'),
                **vkw, **kwargs, **self.kwargs,
            )

            print("  Collecting results from npy_stack...")
            result_datasets = []
            for var, ba_path, sf_cube in zip(self.variables, ba_paths, sf_cubes):
                npy_dir = uf.npy_stack_dir(str(ba_path))
                # Assemble a contiguous array from the per-location npy stack. Each
                # {flat}.npy holds one location's (n_time, 1) series; flat is the
                # row-major index over the spatial dims. This replaces
                # da.from_npy_stack(...).reshape(...), whose lazy graph yielded
                # all-NaN with this dask/iris stack. Building the xarray EAGERLY
                # from the numpy array (no iris.save -> lazy load_cube roundtrip)
                # also avoids returning data bound to the temp file that the
                # `finally` block below deletes.
                arr = np.full(sf_cube.shape, np.nan, dtype=np.float32)
                flat_view = arr.reshape(sf_cube.shape[0], -1)
                for flat in range(flat_view.shape[1]):
                    fp = f"{npy_dir}{flat}.npy"
                    if os.path.exists(fp):
                        flat_view[:, flat] = np.load(fp).squeeze()
                out_cube = sf_cube.copy()
                out_cube.data = np.ma.masked_invalid(arr)
                result_datasets.append(
                    xr.DataArray.from_iris(out_cube).to_dataset(name=var))

            result = xr.merge(result_datasets)
            result.attrs['bias_correction_method'] = (
                'ISIMIP3BASD-MBCn' if self.n_iterations > 0 else 'ISIMIP3BASD')
            result.attrs['n_iterations'] = self.n_iterations
            result.attrs['variables']    = str(self.variables)

            if output_path is not None:
                result.to_netcdf(output_path)
                print(f"  Saved to: {output_path}")

        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

        return result


# ---------------------------------------------------------------------------
# Statistical Downscaling
# ---------------------------------------------------------------------------

class StatisticalDownscaling:
    """Statistical downscaling using the ISIMIP3BASD modified MBCn algorithm.

    Downscales bias-corrected coarse-resolution data to the fine resolution
    of the reference observations, fitting multi-decade windows (``window_years``)
    and writing one file per year.  Already-completed years are skipped.

    Args:
        variable (str): Climate variable name.
        n_processes (int, optional): Parallel worker count.
        n_iterations (int): MBCn iterations (default 20).
        randomization_seed (int, optional): RNG seed.
        lower_bound, lower_threshold, upper_bound, upper_threshold (float, optional):
            Physical bounds / thresholds (defaults from BiasCorrection.DEFAULT_CONFIGS).
        load_obs_in_memory (bool): Read ``obs_fine`` into memory once before
            the year loop instead of re-reading it lazily for every coarse
            cell of every year (default True). Needs RAM for one variable's
            full obs_fine array; set False for domains too large for that.
        window_years (int): Length of the fitting window (default 30). The
            years are split into equal chunks of about this length and each
            chunk is downscaled in one call; output is still one file per year.
            Short windows (e.g. 1 year) bias pr high and inflate variance.

    Example:
        >>> sd = StatisticalDownscaling(variable='tas', n_iterations=20)
        >>> paths = sd.downscale(obs_fine=obs_fine_ds, sim_coarse=gcm_corrected_ds,
        ...                      output_path='./downscaled/', model='ACCESS', scenario='ssp585')
    """

    DEFAULT_CONFIGS = BiasCorrection.DEFAULT_CONFIGS

    def __init__(
        self,
        variable: str,
        n_processes: Optional[int] = None,
        n_iterations: int = 20,
        randomization_seed: Optional[int] = None,
        lower_bound: Optional[float] = None,
        lower_threshold: Optional[float] = None,
        upper_bound: Optional[float] = None,
        upper_threshold: Optional[float] = None,
        load_obs_in_memory: bool = True,
        window_years: int = 30,
        **kwargs
    ):
        """Configure a downscaling run.

        Every argument is documented on the class. Arguments left as ``None``
        fall back to the variable's entry in :attr:`DEFAULT_CONFIGS`.

        Raises:
            ImportError: If ``iris`` is not installed.
        """
        self.variable        = variable
        self.n_processes     = n_processes if n_processes is not None else os.cpu_count()
        self.n_iterations    = n_iterations
        self.randomization_seed = randomization_seed
        self.load_obs_in_memory = load_obs_in_memory
        self.window_years    = window_years
        self.kwargs          = kwargs

        cfg = self.DEFAULT_CONFIGS.get(variable, {})
        self.lower_bound     = lower_bound     if lower_bound     is not None else cfg.get('lower_bound')
        self.lower_threshold = lower_threshold if lower_threshold is not None else cfg.get('lower_threshold')
        self.upper_bound     = upper_bound     if upper_bound     is not None else cfg.get('upper_bound')
        self.upper_threshold = upper_threshold if upper_threshold is not None else cfg.get('upper_threshold')
        # With halfwin_upper_bound_climatology (rsds) the bias-adjustment upper
        # bound is a factor on a running climatology of daily maxima, not a value.
        # Downscaling has no such climatology and would clip at the factor itself
        # (rsds at 1 W m-2), so the default upper bound is dropped here.
        if cfg.get('halfwin_upper_bound_climatology') and upper_bound is None:
            self.upper_bound = self.upper_threshold = None

        print(f"StatisticalDownscaling initialized for {variable}")
        print(f"  n_iterations={n_iterations}  n_processes={self.n_processes}  seed={self.randomization_seed}")

    def downscale(
        self,
        obs_fine: xr.Dataset,
        sim_coarse: xr.Dataset,
        output_path: Optional[Union[str, Path]] = None,
        model: Optional[str] = None,
        scenario: Optional[str] = None,
        **kwargs
    ) -> list:
        """Apply statistical downscaling to coarse-resolution data.

        Args:
            obs_fine (xr.Dataset or DataArray): Fine-resolution observations.
            sim_coarse (xr.Dataset or DataArray): Coarse-resolution bias-corrected
                simulations.
            output_path (str or Path, optional): Directory for per-year output NetCDF
                files.
            model (str, optional): Model tag used in output filenames.
            scenario (str, optional): Scenario tag used in output filenames.

        Returns:
            list of Path: Per-year NetCDF paths
                (BCSD_{var}_{model}_{scenario}_{year}.nc).
        """
        if not IRIS_AVAILABLE:
            raise ImportError("iris is required: pip install scitools-iris")
        try:
            from climdata._vendor.isimip3basd import statistical_downscaling as sd
            from climdata._vendor.isimip3basd import utility_functions as uf
        except ImportError:
            raise ImportError("ISIMIP3BASD vendor code not found in climdata/_vendor/isimip3basd/.")
        import dask.array as da

        if isinstance(obs_fine,   xr.DataArray):
            obs_fine   = obs_fine.to_dataset(name=self.variable)
        if isinstance(sim_coarse, xr.DataArray):
            sim_coarse = sim_coarse.to_dataset(name=self.variable)

        print(f"\nStatistical downscaling: {self.variable}")

        tmpdir = Path.cwd() / "tmp_bcsd_downscale"
        tmpdir.mkdir(parents=True, exist_ok=True)
        try:
            # Convert to iris cubes
            obs_fine_cube   = _xda_to_cube(obs_fine[self.variable])
            sim_coarse_cube = _xda_to_cube(sim_coarse[self.variable])
            _stamp_geogcs(obs_fine_cube)
            _stamp_geogcs(sim_coarse_cube)
            obs_fine_cube = _ensure_rectilinear_iris_cube(obs_fine_cube)
            print(f"  obs_fine {obs_fine_cube.shape}  sim_coarse {sim_coarse_cube.shape}")

            # Crop obs_fine to integer multiple of the coarse grid
            lat_name = 'lat' if 'lat' in sim_coarse.coords else 'latitude'
            lon_name = 'lon' if 'lon' in sim_coarse.coords else 'longitude'
            nlat_c, nlon_c = sim_coarse.sizes[lat_name], sim_coarse.sizes[lon_name]
            nlat_f, nlon_f = obs_fine_cube.shape[1], obs_fine_cube.shape[2]
            nlat_crop = (nlat_f // nlat_c) * nlat_c
            nlon_crop = (nlon_f // nlon_c) * nlon_c
            if nlat_crop != nlat_f or nlon_crop != nlon_f:
                r = (nlat_f - nlat_crop) // 2
                c = (nlon_f - nlon_crop) // 2
                obs_fine_cube = obs_fine_cube[:, r:r + nlat_crop, c:c + nlon_crop]
                print(f"  obs_fine cropped to ({nlat_crop}×{nlon_crop}) "
                      f"[was {nlat_f}×{nlon_f}, coarse {nlat_c}×{nlon_c}]")

            # Realize obs_fine once, up front. The vendor reads every coarse
            # cell's obs block with .compute() — once per cell per year — and
            # a lazy, store-backed obs_fine makes each of those reads decompress
            # (and masked_invalid) whole Zarr chunks, e.g. (13148, 100, 100) =
            # 526 MB for a 10×12 block: ~60% of the per-cell time, repeated for
            # every year. Held in memory, a read is just a slice. Still wrapped
            # in dask because the vendor calls .compute(); worker processes are
            # forked, so they share this array rather than copying it.
            if self.load_obs_in_memory:
                _t0 = __import__("time").time()
                obs_arr = np.ma.masked_invalid(obs_fine_cube.core_data().compute())
                obs_fine_cube.data = da.from_array(obs_arr, chunks=-1, asarray=False)
                print(f"  obs_fine loaded into memory: {obs_arr.nbytes / 1e9:.1f} GB "
                      f"({__import__('time').time() - _t0:.0f}s)")

            # Output directory and filename template
            out_dir = Path(output_path).resolve() if output_path else tmpdir
            out_dir.mkdir(parents=True, exist_ok=True)
            _model_tag    = model    or "model"
            _scenario_tag = scenario or "scenario"

            # Fit window. The vendor algorithm quantile-maps the simulated
            # fine-cell series against ~36 years of obs per calendar month and
            # clips values at the bounds; with only one year (~30 days per
            # month) per call it badly inflates pr (+13…+118 % in a
            # perfect-model test) and temperature variance. Each call therefore
            # sees a multi-decade window (years split into equal chunks of
            # about window_years), and results are still written one file per
            # year so the on-disk layout and resume logic are unchanged.
            tc    = sim_coarse_cube.coord("time")
            years = np.array([d.year for d in tc.units.num2date(tc.points)])
            unique_years = np.unique(years)
            n_win   = max(1, int(np.ceil(len(unique_years) / self.window_years)))
            windows = np.array_split(unique_years, n_win)
            print(f"  Processing {len(unique_years)} years: {unique_years[0]}–{unique_years[-1]} "
                  f"in {n_win} windows of ~{len(windows[0])} years")

            def _yr_path(yr):
                return (out_dir / f"BCSD_{self.variable}_{_model_tag}_{_scenario_tag}_{yr}.nc").resolve()

            try:
                from tqdm import tqdm as _tqdm
                win_iter = _tqdm(windows, desc="  Downscaling", unit="win", dynamic_ncols=True)
            except Exception:
                win_iter = windows

            yearly_result_paths = []
            for win in win_iter:
                todo = [yr for yr in win if not _yr_path(yr).exists()]
                if not todo:
                    yearly_result_paths += [_yr_path(yr) for yr in win]
                    continue

                idx = np.where(np.isin(years, win))[0]
                win_npy_root     = tmpdir / f"npy_{win[0]}_{win[-1]}"
                win_npy_root.mkdir(parents=True, exist_ok=True)
                win_npy_sentinel = str(win_npy_root / "sim_fine.nc")
                uf.setup_npy_stack(win_npy_sentinel, (len(idx),) + obs_fine_cube.shape[1:])

                sc_win   = _ensure_rectilinear_iris_cube(sim_coarse_cube[idx])
                # Vendor uf.remapbil bilinearly interpolates coarse → fine grid.
                # iris Linear regrid realizes the data to a numpy MaskedArray, but
                # the vendor sd.downscale calls .compute() on every cube's
                # core_data(), so re-wrap as a dask masked array (as _xda_to_cube does).
                #
                # Both arrays are realized and masked once, here, then wrapped
                # with da.from_array: a worker's .compute() of its block is then
                # just a slice. Wrapping with da.ma.masked_invalid instead leaves
                # one whole-window task in the graph, so every forked worker
                # re-masks the full window (~10 GB for 28 years at 0.1°) to read
                # its block — 40 workers OOM-killed and the pool hung.
                sc_win.data = da.from_array(
                    np.ma.masked_invalid(np.ma.asarray(sc_win.core_data().compute()
                                                       if hasattr(sc_win.core_data(), "compute")
                                                       else sc_win.core_data())),
                    chunks=-1, asarray=False)
                remap_win = uf.remapbil(sc_win, obs_fine_cube)
                _r = remap_win.core_data()
                _r = _r.compute() if hasattr(_r, "compute") else _r
                remap_win.data = da.from_array(np.ma.masked_invalid(_r), chunks=-1, asarray=False)
                del _r
                _stamp_geogcs(remap_win)

                t0 = __import__("time").time()
                sd.downscale(
                    obs_fine=obs_fine_cube,
                    sim_coarse=sc_win,
                    sim_coarse_remapbil=remap_win,
                    sim_fine_path=win_npy_sentinel,
                    n_processes=self.n_processes,
                    n_iterations=self.n_iterations,
                    randomization_seed=self.randomization_seed,
                    lower_bound=self.lower_bound,
                    lower_threshold=self.lower_threshold,
                    upper_bound=self.upper_bound,
                    upper_threshold=self.upper_threshold,
                    **kwargs, **self.kwargs,
                )
                elapsed = __import__("time").time() - t0
                del remap_win

                # npy_stack → contiguous array. Each {flat}.npy holds one fine
                # location's (n_time, 1) series; flat is the row-major spatial
                # index. A contiguous array writes compact files quickly, unlike
                # da.from_npy_stack fed straight to iris.save.
                npy_dir = uf.npy_stack_dir(win_npy_sentinel)
                arr = np.full((len(idx),) + obs_fine_cube.shape[1:], np.nan, dtype=np.float32)
                flat_view = arr.reshape(len(idx), -1)
                for flat in range(flat_view.shape[1]):
                    fp = f"{npy_dir}{flat}.npy"
                    if os.path.exists(fp):
                        flat_view[:, flat] = np.load(fp).squeeze()
                shutil.rmtree(str(win_npy_root), ignore_errors=True)

                # Split the window into yearly NetCDFs (only years still missing).
                win_years = years[idx]
                for yr in win:
                    yr_path = _yr_path(yr)
                    yearly_result_paths.append(yr_path)
                    if yr_path.exists():
                        continue
                    k = np.where(win_years == yr)[0]
                    sim_fine_yr = obs_fine_cube[:len(k)].copy()
                    sim_fine_yr.data = np.ma.masked_invalid(arr[k])
                    # Output carries the (application) sim_coarse times, not obs_fine's.
                    _t, _src = sim_fine_yr.coord("time"), sc_win.coord("time")
                    _t.points, _t.units = _src.points[k], _src.units
                    iris.save(sim_fine_yr, str(yr_path),
                              saver=iris.fileformats.netcdf.save,
                              unlimited_dimensions=["time"], zlib=True, complevel=1)
                del arr
                print(f"  {win[0]}–{win[-1]} → {len(todo)} yearly files  ({elapsed:.0f}s)", flush=True)

            print("  Statistical downscaling complete.")

        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

        return yearly_result_paths


# ---------------------------------------------------------------------------
# Combined BCSD workflow
# ---------------------------------------------------------------------------

class BCSD:
    """Complete Bias Correction and Statistical Downscaling (BCSD) workflow.

    Combines BiasCorrection and StatisticalDownscaling in a single pipeline
    following the ISIMIP3BASD methodology.

    Passing a list of variables runs the bias-correction stage **multivariately**
    (MBCn over the joint distribution, so inter-variable dependence — e.g. the
    temperature/precipitation correlation, or tasmax ≥ tasmin — is carried into
    the corrected fields), then downscales each corrected variable
    **univariately**, one at a time.  This asymmetry is deliberate: the vendored
    ISIMIP3BASD downscaling routine is single-variable by construction.

    Args:
        variable (str or list of str): Climate variable name(s). A list enables
            multivariate (MBCn) bias correction across all of them at once.
        bias_correction_kwargs (dict, optional): Forwarded to BiasCorrection.
        downscaling_kwargs (dict, optional): Forwarded to every per-variable
            StatisticalDownscaling. Variable-specific bounds and thresholds are
            not needed here — they come from
            :attr:`BiasCorrection.DEFAULT_CONFIGS` per variable.
        regridding_method (str): xESMF method for obs_fine→coarse regridding (default
            'conservative').
        regridding_tool (str): 'xesmf' (default) or 'cdo'.
        cdo_method (str): CDO operator when regridding_tool='cdo'.
        cdo_env (str): Conda environment with CDO.
        weights_dir (str, optional): Directory to cache regridding weights.

    Attributes:
        variables (list of str): Every variable this pipeline handles.
        variable (str or None): The single variable, or None for a multivariate
            pipeline.
        downscalers (dict): Variable name → its StatisticalDownscaling instance.
        downscaling (StatisticalDownscaling or None): The single downscaler, or
            None for a multivariate pipeline (use :attr:`downscalers`).

    Example:
        >>> bcsd = BCSD(variable='tas',
        ...             bias_correction_kwargs={'detrend': True},
        ...             downscaling_kwargs={'n_iterations': 20})
        >>> paths = bcsd.run(obs_fine=obs_025deg,
        ...                  sim_hist_coarse=gcm_hist_1deg,
        ...                  sim_fut_coarse=gcm_fut_1deg)

        Multivariate bias correction, univariate downscaling:

        >>> bcsd = BCSD(variable=['pr', 'tas', 'tasmax', 'tasmin'],
        ...             bias_correction_kwargs={'n_iterations': 20})
        >>> paths = bcsd.run(obs_fine=obs_025deg,
        ...                  sim_hist_coarse=gcm_hist_1deg,
        ...                  sim_fut_coarse=gcm_fut_1deg)     # doctest: +SKIP
        >>> sorted(paths)                                      # doctest: +SKIP
        ['pr', 'tas', 'tasmax', 'tasmin']
    """

    def __init__(
        self,
        variable: Union[str, List[str]],
        bias_correction_kwargs: Optional[Dict] = None,
        downscaling_kwargs: Optional[Dict] = None,
        regridding_method: str = 'conservative',
        regridding_tool: str = 'xesmf',
        cdo_method: str = 'remapcon',
        cdo_env: str = 'cdo_stable',
        weights_dir: Optional[str] = None
    ):
        """Configure the two stages of a BCSD run.

        Every argument is documented on the class. The stage objects are built
        here, so an invalid ``bias_correction_kwargs`` or ``downscaling_kwargs``
        is rejected now rather than partway through :meth:`run`.

        Raises:
            ValueError: If no variable, or a duplicated variable, is given.
            ImportError: If ``iris`` is not installed.
        """
        self.variables = [variable] if isinstance(variable, str) else list(variable)
        if not self.variables:
            raise ValueError("BCSD needs at least one variable.")
        if len(set(self.variables)) != len(self.variables):
            raise ValueError(f"Duplicate variables in {self.variables}.")
        # None for a multivariate pipeline, mirroring BiasCorrection.variable
        self.variable          = self.variables[0] if len(self.variables) == 1 else None
        self.regridding_method = regridding_method
        self.regridding_tool   = regridding_tool
        self.cdo_method        = cdo_method
        self.cdo_env           = cdo_env
        self.weights_dir       = weights_dir

        # One bias corrector over all variables jointly (MBCn when >1) ...
        self.bias_correction = BiasCorrection(variable=self.variables,
                                              **(bias_correction_kwargs or {}))
        # ... but one downscaler per variable: ISIMIP3BASD downscaling is
        # single-variable, and each needs its own bounds from DEFAULT_CONFIGS.
        self.downscalers = {
            v: StatisticalDownscaling(variable=v, **(downscaling_kwargs or {}))
            for v in self.variables
        }
        self.downscaling = self.downscalers[self.variable] if self.variable else None

        mode = ('multivariate (MBCn) bias correction + univariate downscaling'
                if len(self.variables) > 1 else 'univariate')
        print(f"\n{'='*60}")
        print(f"BCSD Pipeline initialized for {self.variables}")
        print(f"Mode: {mode}")
        print(f"Regridding: {regridding_tool} ({regridding_method})")
        print(f"{'='*60}")

    # ------------------------------------------------------------------ helpers

    def _as_dataset(self, data, name: str) -> xr.Dataset:
        """Accept a DataArray for a single-variable pipeline; require a Dataset
        otherwise, since a bare array cannot name several variables."""
        if isinstance(data, xr.DataArray):
            if self.variable is None:
                raise TypeError(
                    f"{name} must be an xr.Dataset holding {self.variables} — a "
                    "DataArray only carries one variable."
                )
            return data.to_dataset(name=self.variable)
        return data

    def _select(self, ds: xr.Dataset, name: str) -> xr.Dataset:
        """Subset to this pipeline's variables, failing loudly on a missing one."""
        missing = [v for v in self.variables if v not in ds.data_vars]
        if missing:
            raise KeyError(
                f"{name} is missing variable(s) {missing} — BCSD is configured for "
                f"{self.variables}, and multivariate bias correction needs every one "
                f"of them on the same time axis. Present: {sorted(ds.data_vars)}."
            )
        return ds[self.variables]

    @property
    def _tag(self) -> str:
        """Filename stem for intermediate files: the variable, or all of them."""
        return "_".join(self.variables)

    def run(
        self,
        obs_fine: xr.Dataset,
        sim_hist_coarse: xr.Dataset,
        sim_fut_coarse: xr.Dataset,
        obs_hist_coarse: Optional[xr.Dataset] = None,
        output_path: Optional[Union[str, Path]] = None,
        save_intermediate: bool = False,
        model: Optional[str] = None,
        scenario: Optional[str] = None,
    ) -> Union[list, Dict[str, list]]:
        """Run the complete BCSD workflow.

        Three steps, of which the first is conditional:

        0. If ``obs_hist_coarse`` is not supplied, regrid ``obs_fine`` onto the
           simulation's coarse grid by conservative remapping.
        1. Bias-correct the coarse simulation against the coarse observations
           using ISIMIP3BASD — jointly over every variable (MBCn) when the
           pipeline was built with a list.
        2. Downscale each corrected coarse field onto the fine grid, again with
           ISIMIP3BASD, writing one NetCDF per variable per year.

        Args:
            obs_fine (xr.Dataset): Historical observations at fine resolution.
            sim_hist_coarse (xr.Dataset): Historical GCM simulations at coarse
                resolution.
            sim_fut_coarse (xr.Dataset): Future GCM simulations at coarse resolution.
            obs_hist_coarse (xr.Dataset, optional): Historical observations already at
                coarse resolution. When None, derived automatically via conservative
                regridding.
            output_path (str or Path, optional): Directory for per-year downscaled
                NetCDF files. Per-year filenames carry the variable name, so all
                variables of a multivariate run can share one directory.
            save_intermediate (bool): Save intermediate bias-corrected output alongside
                final results.
            model (str, optional): Model tag used in output filenames.
            scenario (str, optional): Scenario tag used in output filenames.

        Returns:
            list[Path] for a single-variable pipeline: the downscaled NetCDF
            written for each year. dict[str, list[Path]] for a multivariate one:
            that same per-year list, keyed by variable.

        Raises:
            KeyError: If any input is missing one of the pipeline's variables.
        """
        print(f"\n{'='*60}\nStarting BCSD Workflow\n{'='*60}\n")

        obs_fine        = self._select(self._as_dataset(obs_fine, "obs_fine"), "obs_fine")
        sim_hist_coarse = self._select(self._as_dataset(sim_hist_coarse, "sim_hist_coarse"),
                                       "sim_hist_coarse")
        sim_fut_coarse  = self._select(self._as_dataset(sim_fut_coarse, "sim_fut_coarse"),
                                       "sim_fut_coarse")

        # Step 0: Derive coarse observations if not supplied
        if obs_hist_coarse is None:
            print("STEP 0: Regridding obs_fine → coarse resolution")
            print("-" * 60)
            obs_hist_coarse = regrid_to_coarse(
                obs_fine, sim_hist_coarse,
                method=self.regridding_method,
                regridding_tool=self.regridding_tool,
                cdo_method=self.cdo_method,
                cdo_env=self.cdo_env,
                weights_dir=self.weights_dir,
            )
            if save_intermediate and output_path:
                p = Path(output_path) / f"{self._tag}_obs_coarse.nc"
                obs_hist_coarse.to_netcdf(p)
                print(f"  Saved coarse obs to: {p}")
        else:
            print("Using provided coarse observations.\n")
        obs_hist_coarse = self._select(self._as_dataset(obs_hist_coarse, "obs_hist_coarse"),
                                       "obs_hist_coarse")

        # Step 1: Bias correction — one call over all variables at once, so MBCn
        # sees the joint distribution rather than each variable in isolation.
        print("\nSTEP 1: BIAS CORRECTION\n" + "-" * 60)
        bc_output = None
        if save_intermediate and output_path:
            bc_output = Path(output_path) / f"{self._tag}_bias_corrected.nc"

        sim_fut_corrected = self.bias_correction.correct(
            obs_hist=obs_hist_coarse,
            sim_hist=sim_hist_coarse,
            sim_fut=sim_fut_coarse,
            output_path=bc_output,
        )

        # Step 2: Statistical downscaling — univariate, one variable at a time.
        print("\nSTEP 2: STATISTICAL DOWNSCALING\n" + "-" * 60)
        results: Dict[str, list] = {}
        for i, var in enumerate(self.variables, start=1):
            if len(self.variables) > 1:
                print(f"\n  [{i}/{len(self.variables)}] {var}")
            results[var] = self.downscalers[var].downscale(
                obs_fine=obs_fine[[var]],
                sim_coarse=sim_fut_corrected[[var]],
                output_path=output_path,
                model=model,
                scenario=scenario,
            )

        print(f"\n{'='*60}\nBCSD Workflow Complete!\n{'='*60}\n")
        return results[self.variable] if self.variable else results
