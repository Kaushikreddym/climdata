"""
Bias Correction and Statistical Downscaling (BCSD)

This module provides a high-level interface to the ISIMIP3BASD methods
for climate data downloaded via climdata. It wraps the ISIMIP3BASD
functionality to work seamlessly with xarray datasets.

The BCSD method performs:
1. Bias Correction: Adjusts systematic biases in climate model outputs
2. Statistical Downscaling: Increases spatial resolution using fine-scale observations

Reference:
    Lange, S. (2019). Trend-preserving bias adjustment and statistical downscaling
    with ISIMIP3BASD (v1.0). Geoscientific Model Development, 12(7), 3055-3070.
    https://doi.org/10.5194/gmd-12-3055-2019
"""

import xarray as xr
import numpy as np
from pathlib import Path
from typing import Optional, Union, Dict, List, Literal
import warnings
import tempfile
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


def _to_proleptic_gregorian(cube):
    """
    Reinterpret an iris cube's time coordinate as proleptic_gregorian calendar.

    ISIMIP3BASD asserts this calendar on every input cube. Most climate datasets
    use 'gregorian' or 'standard', which are functionally equivalent for dates
    after 1582, so a simple unit-string replacement is safe.
    """
    import cf_units
    t = cube.coord('time')
    if t.units.calendar != 'proleptic_gregorian':
        new_units = cf_units.Unit(
            t.units.origin,          # keep the epoch string (e.g. 'days since ...')
            calendar='proleptic_gregorian'
        )
        t.units = new_units
    return cube


def _is_curvilinear(cube):
    """Return True if cube's horizontal coords are 2D AuxCoords (curvilinear/rotated-pole)."""
    try:
        x_dim_coords = cube.coords(axis='x', dim_coords=True)
        y_dim_coords = cube.coords(axis='y', dim_coords=True)
        return len(x_dim_coords) != 1 or len(y_dim_coords) != 1
    except Exception:
        return True


def _stamp_geogcs(cube):
    """
    Attach a GeogCS to every 1-D latitude/longitude DimCoord and normalise
    their units to 'degrees'.

    iris.analysis.Linear (RectilinearRegridder) has two distinct checks:

    1. When *neither* cube has a CRS it compares metadata strings
       (standard_name, units, var_name) exactly.  xarray→iris produces
       'degrees_north'/'degrees_east', while _regrid_curvilinear_to_rectilinear
       produces 'degrees' — the mismatch raises a ValueError.

    2. Once a coord *has* a GeogCS, RectilinearRegridder._check_units strictly
       requires units == 'degrees' (not 'degrees_north', not 'degrees_east').
       Leaving the original units in place therefore causes a second ValueError.

    This function fixes both by stamping the same GeogCS *and* setting
    units='degrees' on every horizontal DimCoord.
    """
    import iris.coord_systems as ics
    import cf_units
    geog_cs = ics.GeogCS(6371229.0)   # standard sphere used by most climate models
    deg = cf_units.Unit('degrees')
    for coord in cube.coords(dim_coords=True):
        axis = iris.util.guess_coord_axis(coord)
        if axis in ('X', 'Y'):
            coord.coord_system = geog_cs
            coord.units = deg
    return cube


def _regrid_curvilinear_to_rectilinear(cube):
    """
    Reproject a curvilinear (e.g. rotated-pole) iris cube to a regular
    geographic lat/lon grid using scipy griddata.

    The output resolution is chosen to approximately preserve the number of
    grid cells in each direction.  A 1D latitude DimCoord and 1D longitude
    DimCoord are created so that iris.analysis.Linear can subsequently use
    the resulting cube as a regrid target.

    Parameters
    ----------
    cube : iris.cube.Cube
        Input cube with 2D AuxCoords for latitude and longitude.

    Returns
    -------
    iris.cube.Cube
        Cube on a regular lat/lon grid with 1D DimCoords.
    """
    from scipy.interpolate import griddata
    import iris.coords
    import iris.cube

    # ── Extract 2-D lat/lon values ───────────────────────────────────────────
    try:
        lat2d = cube.coord('latitude').points
        lon2d = cube.coord('longitude').points
    except iris.exceptions.CoordinateNotFoundError:
        # fall back to grid_latitude / grid_longitude (rotated pole)
        lat2d = cube.coord('grid_latitude').points
        lon2d = cube.coord('grid_longitude').points

    ny, nx = lat2d.shape

    # ── Build regular target grid ────────────────────────────────────────────
    lat_min, lat_max = lat2d.min(), lat2d.max()
    lon_min, lon_max = lon2d.min(), lon2d.max()
    reg_lat = np.linspace(lat_min, lat_max, ny)
    reg_lon = np.linspace(lon_min, lon_max, nx)
    grid_lon, grid_lat = np.meshgrid(reg_lon, reg_lat)

    # ── Interpolate each time step ───────────────────────────────────────────
    src_points = np.column_stack([lat2d.ravel(), lon2d.ravel()])
    nt = cube.shape[0]
    out = np.full((nt, ny, nx), np.nan, dtype=np.float32)
    data = cube.data if not hasattr(cube.data, 'compute') else cube.data.compute()
    data = np.ma.filled(np.ma.asarray(data), np.nan)
    for t in range(nt):
        out[t] = griddata(src_points, data[t].ravel(),
                          (grid_lat, grid_lon), method='linear')

    # ── Build output iris cube with 1D DimCoords ─────────────────────────────
    lat_coord = iris.coords.DimCoord(
        reg_lat, standard_name='latitude', units='degrees')
    lon_coord = iris.coords.DimCoord(
        reg_lon, standard_name='longitude', units='degrees')
    time_coord = cube.coord('time')

    out_masked = np.ma.masked_invalid(out)
    new_cube = iris.cube.Cube(
        out_masked,
        standard_name=cube.standard_name,
        long_name=cube.long_name,
        var_name=cube.var_name,
        units=cube.units,
        dim_coords_and_dims=[
            (time_coord.copy(), 0),
            (lat_coord, 1),
            (lon_coord, 2),
        ],
    )
    return new_cube


def regrid_to_coarse(
    fine_data: xr.Dataset,
    coarse_template: xr.Dataset,
    method: str = 'conservative',
    regridding_tool: str = 'xesmf',
    cdo_method: str = 'remapcon',
    cdo_env: str = 'cdo_stable',
    weights_dir: Optional[str] = None
) -> xr.Dataset:
    """
    Regrid fine-resolution data to coarse resolution.
    
    This function is used to create coarse-resolution observations from
    fine-resolution observations to match the GCM grid.
    
    Parameters
    ----------
    fine_data : xr.Dataset
        Fine-resolution dataset to regrid
    coarse_template : xr.Dataset
        Coarse-resolution dataset to use as template
    method : str, optional
        Regridding method:
        - 'conservative': Area-weighted conservative (recommended, default)
        - 'bilinear': Bilinear interpolation
        - 'nearest': Nearest neighbor
        For CDO: 'remapcon', 'remapbil', 'remapdis', 'remapnn'
    regridding_tool : str, optional
        Tool to use: 'xesmf' (Python, default) or 'cdo' (CDO command-line)
    cdo_method : str, optional
        CDO-specific method if regridding_tool='cdo'
    cdo_env : str, optional
        Conda environment with CDO installed
    weights_dir : str, optional
        Directory to save/load regridding weights for reuse
    
    Returns
    -------
    xr.Dataset
        Regridded coarse-resolution dataset
    """
    print(f"🔄 Regridding from fine to coarse resolution using {regridding_tool}...")
    print(f"   Fine grid: {fine_data.dims}")
    print(f"   Target coarse grid: {coarse_template.dims}")
    
    if regridding_tool == 'xesmf':
        if not XESMF_AVAILABLE:
            raise ImportError(
                "xESMF is not installed. Install it with: pip install xesmf\n"
                "Or use regridding_tool='cdo' if CDO is available."
            )
        
        # Create weights file path
        if weights_dir:
            os.makedirs(weights_dir, exist_ok=True)
            weight_file = os.path.join(weights_dir, f'weights_fine_to_coarse_{method}.nc')
        else:
            weight_file = None
        
        # Create regridder
        print(f"   Creating {method} regridder...")
        regridder = xe.Regridder(
            fine_data,
            coarse_template,
            method=method,
            periodic=False,
            ignore_degenerate=True,
            filename=weight_file,
            reuse_weights=weight_file and os.path.exists(weight_file)
        )

        # Build a time-independent valid-source mask:
        # For each coarse cell, record the fraction of source weight that came
        # from actually valid (non-NaN) fine-grid cells.  Any coarse cell that
        # had *no* valid source (fraction == 0) is masked to NaN in the output,
        # preventing xESMF's default behaviour of filling those cells with zero.
        first_var = next(iter(fine_data.data_vars))
        fine_valid = xr.where(fine_data[first_var].isel(time=0).notnull(), 1.0, 0.0)
        coarse_valid_frac = regridder(fine_valid)   # fraction of weight from valid cells
        # cells where coarse_valid_frac == 0 had only NaN sources → mask them
        valid_mask = coarse_valid_frac > 0

        # Regrid each variable and apply the mask
        result_vars = {}
        for var in fine_data.data_vars:
            print(f"   Regridding variable: {var}")
            regridded = regridder(fine_data[var].where(fine_data[var].notnull(), other=0.0))
            # Re-normalise by valid fraction to undo the zero-fill effect on
            # partially-covered border cells, then mask fully-uncovered cells
            regridded = (regridded / coarse_valid_frac).where(valid_mask)
            regridded.attrs = fine_data[var].attrs
            result_vars[var] = regridded

        result = xr.Dataset(result_vars)
        result.attrs = fine_data.attrs
        
    elif regridding_tool == 'cdo':
        # Use CDO for regridding
        with tempfile.TemporaryDirectory() as tmpdir:
            # Save inputs
            fine_file = os.path.join(tmpdir, 'fine_input.nc')
            coarse_template_file = os.path.join(tmpdir, 'coarse_template.nc')
            output_file = os.path.join(tmpdir, 'regridded_output.nc')
            
            fine_data.to_netcdf(fine_file)
            coarse_template.to_netcdf(coarse_template_file)
            
            # Run CDO regridding
            if cdo_method == 'remapcon':
                # Conservative remapping
                method_str = 'remapcon'
            elif cdo_method == 'remapbil':
                method_str = 'remapbil'
            elif cdo_method == 'remapdis':
                method_str = 'remapdis'
            elif cdo_method == 'remapnn':
                method_str = 'remapnn'
            else:
                method_str = cdo_method
            
            cdo_cmd = f"cdo {method_str},{coarse_template_file} {fine_file} {output_file}"
            cmd = f"conda run -n {cdo_env} {cdo_cmd}"
            
            print(f"   Running CDO: {cdo_cmd}")
            try:
                subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
                result = xr.open_dataset(output_file)
            except subprocess.CalledProcessError as e:
                print(f"CDO Error: {e.stderr}")
                raise RuntimeError(f"CDO regridding failed: {e.stderr}")
    
    else:
        raise ValueError(f"Unknown regridding_tool: {regridding_tool}. Use 'xesmf' or 'cdo'.")
    
    print(f"   ✅ Regridding complete!")
    return result


class BiasCorrection:
    """
    Bias correction using ISIMIP3BASD trend-preserving quantile mapping.
    
    This class applies bias correction at the coarse GCM resolution before
    any spatial downscaling.
    
    Parameters
    ----------
    variable : str
        Climate variable name (tas, pr, rsds, hurs, etc.)
    method : str, optional
        Bias correction method. Options:
        - 'parametric': Uses parametric quantile mapping (default)
        - 'non-parametric': Uses empirical quantile mapping
    distribution : str, optional
        Distribution for parametric QM. Options:
        - None: Non-parametric QM
        - 'normal': For temperature variables
        - 'gamma': For precipitation
        - 'beta': For bounded variables (humidity, radiation ratios)
        - 'weibull': For wind speed
        - 'rice': For temperature ranges
    trend_preservation : str, optional
        How to preserve climate change trends:
        - 'additive': For temperature (default for tas)
        - 'multiplicative': For ratios
        - 'mixed': For precipitation and wind
        - 'bounded': For variables with physical bounds
    detrend : bool, optional
        Whether to remove trends before correction (recommended for temperature)
    adjust_p_values : bool, optional
        Adjust p-values for perfect match in reference period
    lower_bound : float, optional
        Physical lower bound for the variable
    lower_threshold : float, optional
        Values below this are set to lower_bound
    upper_bound : float, optional
        Physical upper bound for the variable
    upper_threshold : float, optional
        Values above this are set to upper_bound
    n_processes : int, optional
        Number of parallel processes (default: 1)
    n_quantiles : int, optional
        Number of quantiles for non-parametric QM (default: 50)
    
    Examples
    --------
    >>> # Temperature bias correction
    >>> bc = BiasCorrection(
    ...     variable='tas',
    ...     distribution='normal',
    ...     trend_preservation='additive',
    ...     detrend=True
    ... )
    >>> tas_corrected = bc.correct(
    ...     obs_hist=obs_ds,
    ...     sim_hist=gcm_hist_ds,
    ...     sim_fut=gcm_fut_ds
    ... )
    
    >>> # Precipitation bias correction
    >>> bc = BiasCorrection(
    ...     variable='pr',
    ...     distribution='gamma',
    ...     trend_preservation='mixed',
    ...     lower_bound=0,
    ...     lower_threshold=0.1
    ... )
    >>> pr_corrected = bc.correct(obs_hist, sim_hist, sim_fut)
    """
    
    # Default configurations for common variables
    DEFAULT_CONFIGS = {
        'tas': {
            'distribution': 'normal',
            'trend_preservation': 'additive',
            'detrend': True,
            'adjust_p_values': False
        },
        'tasmax': {
            'distribution': 'normal',
            'trend_preservation': 'additive',
            'detrend': True,
            'adjust_p_values': False
        },
        'tasmin': {
            'distribution': 'normal',
            'trend_preservation': 'additive',
            'detrend': True,
            'adjust_p_values': False
        },
        'pr': {
            'distribution': 'gamma',
            'trend_preservation': 'mixed',
            'lower_bound': 0,
            'lower_threshold': 0.1,
            'adjust_p_values': True
        },
        'rsds': {
            'distribution': 'beta',
            'trend_preservation': 'bounded',
            'lower_bound': 0,
            'lower_threshold': 0.01,
            'upper_bound': 1,
            'upper_threshold': 0.9999,
            'adjust_p_values': True,
            'halfwin_upper_bound_climatology': 15
        },
        'hurs': {
            'distribution': 'beta',
            'trend_preservation': 'bounded',
            'lower_bound': 0,
            'lower_threshold': 0.01,
            'upper_bound': 100,
            'upper_threshold': 99.99,
            'adjust_p_values': True
        },
        'sfcWind': {
            'distribution': 'weibull',
            'trend_preservation': 'mixed',
            'lower_bound': 0,
            'lower_threshold': 0.01,
            'adjust_p_values': True
        },
        'psl': {
            'distribution': 'normal',
            'trend_preservation': 'additive',
            'adjust_p_values': True,
            'detrend': True
        },
        'rlds': {
            'distribution': 'normal',
            'trend_preservation': 'additive',
            'adjust_p_values': True,
            'detrend': True
        }
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
        n_processes: int = 1,
        n_quantiles: int = 50,
        n_iterations: int = 20,
        randomization_seed: Optional[int] = None,
        **kwargs
    ):
        # Accept a single variable name or a list for multivariate MBCn
        self.variables = [variable] if isinstance(variable, str) else list(variable)
        # Keep .variable for backward-compat (single-var case)
        self.variable = self.variables[0] if len(self.variables) == 1 else None
        self.method = method
        self.n_processes = n_processes
        self.n_quantiles = n_quantiles
        self.n_iterations = n_iterations
        self.randomization_seed = randomization_seed

        # Build per-variable config lists from DEFAULT_CONFIGS
        # (override only if an explicit scalar value was supplied)
        def _pick(key, explicit_val, default_val=None):
            """Return explicit_val when it differs from the sentinel, else default_val."""
            return explicit_val if explicit_val is not None else default_val

        self._var_configs = []
        for v in self.variables:
            cfg = self.DEFAULT_CONFIGS.get(v, {})
            self._var_configs.append({
                'distribution':                  distribution         if distribution         is not None  else cfg.get('distribution'),
                'trend_preservation':            trend_preservation   if trend_preservation   is not None  else cfg.get('trend_preservation', 'additive'),
                'detrend':                       detrend              if detrend              is not False else cfg.get('detrend', False),
                'adjust_p_values':               adjust_p_values      if adjust_p_values      is not False else cfg.get('adjust_p_values', False),
                'lower_bound':                   lower_bound          if lower_bound          is not None  else cfg.get('lower_bound'),
                'lower_threshold':               lower_threshold      if lower_threshold      is not None  else cfg.get('lower_threshold'),
                'upper_bound':                   upper_bound          if upper_bound          is not None  else cfg.get('upper_bound'),
                'upper_threshold':               upper_threshold      if upper_threshold      is not None  else cfg.get('upper_threshold'),
                'halfwin_upper_bound_climatology': halfwin_upper_bound_climatology if halfwin_upper_bound_climatology != 0 else cfg.get('halfwin_upper_bound_climatology', 0),
            })

        # Expose scalar attrs for single-variable backward compatibility
        if len(self.variables) == 1:
            c = self._var_configs[0]
            self.distribution                   = c['distribution']
            self.trend_preservation             = c['trend_preservation']
            self.detrend                        = c['detrend']
            self.adjust_p_values                = c['adjust_p_values']
            self.lower_bound                    = c['lower_bound']
            self.lower_threshold                = c['lower_threshold']
            self.upper_bound                    = c['upper_bound']
            self.upper_threshold                = c['upper_threshold']
            self.halfwin_upper_bound_climatology = c['halfwin_upper_bound_climatology']

        # Store additional kwargs
        self.kwargs = kwargs

        mv_label = '(MBCn multivariate)' if (self.n_iterations > 0 and len(self.variables) > 1) else \
                   '(MBCn univariate)'   if self.n_iterations > 0 else '(univariate, no MBCn)'
        print(f"🔧 BiasCorrection initialized for {self.variables}")
        print(f"   n_iterations : {self.n_iterations}  {mv_label}")
        print(f"   n_processes  : {self.n_processes}")
        print(f"   randomization_seed: {self.randomization_seed}")
        for v, c in zip(self.variables, self._var_configs):
            print(f"   [{v}] dist={c['distribution']}  trend={c['trend_preservation']}  "
                  f"detrend={c['detrend']}  adjust_p={c['adjust_p_values']}")
    
    def correct(
        self,
        obs_hist: xr.Dataset,
        sim_hist: xr.Dataset,
        sim_fut: xr.Dataset,
        output_path: Optional[Union[str, Path]] = None,
        **kwargs
    ) -> xr.Dataset:
        """
        Apply bias correction to climate model data.
        
        Parameters
        ----------
        obs_hist : xr.Dataset
            Historical observations (e.g., 1980-2014)
        sim_hist : xr.Dataset
            Historical model simulations matching obs_hist period
        sim_fut : xr.Dataset
            Future model simulations to be corrected
        output_path : str or Path, optional
            Path to save corrected output. If None, returns xarray Dataset
        **kwargs
            Additional arguments passed to ISIMIP3BASD
        
        Returns
        -------
        xr.Dataset
            Bias-corrected future simulations
        """
        try:
            from climdata._vendor.isimip3basd import bias_adjustment as ba
        except ImportError:
            raise ImportError(
                "ISIMIP3BASD code not found in _vendor directory. "
                "Please download bias_adjustment.py, statistical_downscaling.py, "
                "and utility_functions.py from https://github.com/ISI-MIP/isimip3basd "
                "and place them in climdata/_vendor/isimip3basd/"
            )
        
        if not IRIS_AVAILABLE:
            raise ImportError(
                "iris is required for bias correction. Install it with: pip install scitools-iris"
            )
        
        mv_tag = f"{len(self.variables)} variables (MBCn multivariate)" if len(self.variables) > 1 else self.variables[0]
        print(f"\n🔄 Starting bias correction for {mv_tag} ...")
        print(f"   Obs hist period: {obs_hist.time.values[0]} → {obs_hist.time.values[-1]}")
        print(f"   Sim hist period: {sim_hist.time.values[0]} → {sim_hist.time.values[-1]}")
        print(f"   Sim fut  period: {sim_fut.time.values[0]} → {sim_fut.time.values[-1]}")

        # Create temporary directory in current working directory
        tmpdir = Path.cwd() / "tmp_bcsd"
        tmpdir.mkdir(parents=True, exist_ok=True)

        try:
            import dask.array as da
            from climdata._vendor.isimip3basd import utility_functions as uf

            # ── Build per-variable cube lists ────────────────────────────────
            # Convert xarray → iris DIRECTLY in memory (no disk round-trip).
            # .to_iris() avoids writing/reading intermediate NetCDF files which
            # was the cause of the "stuck at converting" bottleneck for large
            # datasets like sim_fut (2015–2100, ~31k daily time steps).
            print("   Converting xarray datasets to iris cubes (in-memory)...")
            obs_hist_cubes, sim_hist_cubes, sim_fut_cubes = [], [], []
            sim_fut_ba_paths = []

            for i, (var, cfg) in enumerate(zip(self.variables, self._var_configs)):
                # Convert to iris and normalise calendar as required by ISIMIP3BASD.
                # Force-compute dask arrays NOW (in the main process) so that worker
                # processes receive plain numpy arrays via the initializer.
                # Without this each of the n_processes workers independently triggers
                # NetCDF reads for every file in the dataset (n_processes × n_files
                # concurrent disk reads), causing heavy I/O contention.
                _obs_cube = _to_proleptic_gregorian(obs_hist[var].to_iris())
                _obs_arr = _obs_cube.core_data()
                _obs_cube.data = np.ma.masked_invalid(
                    _obs_arr.compute() if hasattr(_obs_arr, "compute") else _obs_arr
                )
                obs_hist_cubes.append(_obs_cube)

                _sh_cube = _to_proleptic_gregorian(sim_hist[var].to_iris())
                _sh_arr = _sh_cube.core_data()
                _sh_cube.data = np.ma.masked_invalid(
                    _sh_arr.compute() if hasattr(_sh_arr, "compute") else _sh_arr
                )
                sim_hist_cubes.append(_sh_cube)

                _sf_cube = _to_proleptic_gregorian(sim_fut[var].to_iris())
                _sf_arr = _sf_cube.core_data()
                _sf_cube.data = np.ma.masked_invalid(
                    _sf_arr.compute() if hasattr(_sf_arr, "compute") else _sf_arr
                )
                sim_fut_cubes.append(_sf_cube)
                print(f"   [{var}] cubes materialised (numpy, proleptic_gregorian) — "
                      f"obs_hist {obs_hist[var].shape}, "
                      f"sim_hist {sim_hist[var].shape}, "
                      f"sim_fut  {sim_fut[var].shape}")

                # Output path: per-variable when multiple vars, single path when one var
                if output_path is None:
                    ba_path = (tmpdir / f"sim_fut_ba_{var}.nc").resolve()
                elif len(self.variables) == 1:
                    ba_path = Path(output_path).resolve()
                else:
                    p = Path(output_path)
                    ba_path = (p.parent / f"{p.stem}_{var}{p.suffix}").resolve()

                sim_fut_ba_paths.append(ba_path)

            # ── Run ISIMIP3BASD once for ALL variables jointly ────────────────
            # Passing multiple cubes + n_iterations > 0 activates MBCn copula
            # adjustment, preserving inter-variable dependence structure.
            print(f"   Running ISIMIP3BASD bias adjustment (n_iterations={self.n_iterations}) ...")
            print(f"   (This may take a while for large datasets)")

            ba.adjust_bias(
                obs_hist=obs_hist_cubes,
                sim_hist=sim_hist_cubes,
                sim_fut=sim_fut_cubes,
                sim_fut_ba_path=[str(p) for p in sim_fut_ba_paths],
                n_processes=self.n_processes,
                n_quantiles=self.n_quantiles,
                n_iterations=self.n_iterations,
                randomization_seed=self.randomization_seed,
                distribution=          [c['distribution']                   for c in self._var_configs],
                trend_preservation=    [c['trend_preservation']              for c in self._var_configs],
                detrend=               [c['detrend']                         for c in self._var_configs],
                adjust_p_values=       [c['adjust_p_values']                 for c in self._var_configs],
                halfwin_upper_bound_climatology=[c['halfwin_upper_bound_climatology'] for c in self._var_configs],
                lower_bound=           [c['lower_bound']      for c in self._var_configs],
                lower_threshold=       [c['lower_threshold']  for c in self._var_configs],
                upper_bound=           [c['upper_bound']       for c in self._var_configs],
                upper_threshold=       [c['upper_threshold']   for c in self._var_configs],
                if_all_invalid_use=           [np.nan  for _ in self._var_configs],
                unconditional_ccs_transfer=   [False   for _ in self._var_configs],
                trendless_bound_frequency=    [False   for _ in self._var_configs],
                **kwargs,
                **self.kwargs
            )

            print("   ✅ Bias correction complete!")
            print("   📦 Collecting results from npy_stack...")

            # ── Collect per-variable results and merge ────────────────────────
            result_datasets = []
            for var, ba_path, sim_fut_cube in zip(self.variables, sim_fut_ba_paths, sim_fut_cubes):
                npy_stack_path = uf.npy_stack_dir(str(ba_path))
                d = da.from_npy_stack(npy_stack_path, mmap_mode=None).reshape(sim_fut_cube.shape)

                sim_fut_ba_cube = sim_fut_cube.copy()
                sim_fut_ba_cube.data = np.ma.masked_array(d.compute())

                print(f"   💾 Saving {var} → {ba_path}")
                iris.save(sim_fut_ba_cube, str(ba_path),
                          saver=iris.fileformats.netcdf.save,
                          unlimited_dimensions=['time'],
                          zlib=True, complevel=1)

                result_cube = iris.load_cube(str(ba_path))
                result_cube.data  # force load
                ds_var = xr.DataArray.from_iris(result_cube).to_dataset(name=var)
                ds_var.load()
                result_datasets.append(ds_var)

            # Merge all variables into one Dataset
            result = xr.merge(result_datasets)
            result.attrs['bias_correction_method'] = 'ISIMIP3BASD-MBCn' if self.n_iterations > 0 else 'ISIMIP3BASD'
            result.attrs['n_iterations'] = self.n_iterations
            result.attrs['variables'] = str(self.variables)

            if output_path is not None and len(self.variables) == 1:
                result.to_netcdf(output_path)
                print(f"   💾 Saved to: {output_path}")
            elif output_path is not None:
                # Multi-var: save combined file too
                result.to_netcdf(output_path)
                print(f"   💾 Combined output saved to: {output_path}")

        finally:
            import shutil
            if tmpdir.exists():
                shutil.rmtree(tmpdir, ignore_errors=True)

        return result


class StatisticalDownscaling:
    """
    Statistical downscaling using ISIMIP3BASD modified MBCn algorithm.
    
    This class downscales bias-corrected coarse-resolution climate data
    to fine resolution using high-resolution observations.
    
    Parameters
    ----------
    variable : str
        Climate variable name
    downscaling_factor : tuple of int, optional
        Downscaling factors (lat_factor, lon_factor). If None, automatically computed
    n_processes : int, optional
        Number of parallel processes
    n_iterations : int, optional
        Number of MBCn iterations (default: 20)
    lower_bound : float, optional
        Physical lower bound
    lower_threshold : float, optional
        Lower threshold for censoring
    upper_bound : float, optional
        Physical upper bound
    upper_threshold : float, optional
        Upper threshold for censoring
    
    Examples
    --------
    >>> # Downscale from 1° to 0.25°
    >>> sd = StatisticalDownscaling(
    ...     variable='tas',
    ...     n_iterations=20
    ... )
    >>> tas_fine = sd.downscale(
    ...     obs_fine=obs_fine_ds,
    ...     sim_coarse=gcm_coarse_corrected_ds
    ... )
    """
    
    DEFAULT_CONFIGS = BiasCorrection.DEFAULT_CONFIGS
    
    def __init__(
        self,
        variable: str,
        downscaling_factor: Optional[tuple] = None,
        n_processes: int = 1,
        n_iterations: int = 20,
        randomization_seed: Optional[int] = None,
        lower_bound: Optional[float] = None,
        lower_threshold: Optional[float] = None,
        upper_bound: Optional[float] = None,
        upper_threshold: Optional[float] = None,
        **kwargs
    ):
        self.variable = variable
        self.downscaling_factor = downscaling_factor
        self.n_processes = n_processes
        self.n_iterations = n_iterations
        self.randomization_seed = randomization_seed
        
        # Get default config
        default_config = self.DEFAULT_CONFIGS.get(variable, {})
        self.lower_bound = lower_bound if lower_bound is not None else default_config.get('lower_bound')
        self.lower_threshold = lower_threshold if lower_threshold is not None else default_config.get('lower_threshold')
        self.upper_bound = upper_bound if upper_bound is not None else default_config.get('upper_bound')
        self.upper_threshold = upper_threshold if upper_threshold is not None else default_config.get('upper_threshold')
        
        self.kwargs = kwargs
        
        print(f"🔧 StatisticalDownscaling initialized for {variable}")
        print(f"   n_iterations: {n_iterations}")
        print(f"   n_processes : {n_processes}")
        print(f"   randomization_seed: {self.randomization_seed}")
    
    def downscale(
        self,
        obs_fine: xr.Dataset,
        sim_coarse: xr.Dataset,
        output_path: Optional[Union[str, Path]] = None,
        **kwargs
    ) -> xr.Dataset:
        """
        Apply statistical downscaling to coarse-resolution data.
        
        Parameters
        ----------
        obs_fine : xr.Dataset
            Fine-resolution observations
        sim_coarse : xr.Dataset
            Coarse-resolution bias-corrected simulations
        output_path : str or Path, optional
            Path to save downscaled output
        
        Returns
        -------
        xr.Dataset
            Downscaled simulations at fine resolution
        """
        try:
            from climdata._vendor.isimip3basd import statistical_downscaling as sd
        except ImportError:
            raise ImportError(
                "ISIMIP3BASD code not found. Please download from "
                "https://github.com/ISI-MIP/isimip3basd"
            )
        
        if not IRIS_AVAILABLE:
            raise ImportError(
                "iris is required for statistical downscaling. Install it with: pip install scitools-iris"
            )
        
        print(f"\n🔄 Starting statistical downscaling for {self.variable}...")
        
        # Create temporary directory in current working directory
        tmpdir = Path.cwd() / "tmp_bcsd_downscale"
        tmpdir.mkdir(parents=True, exist_ok=True)
        
        try:
            import dask.array as da

            # Convert xarray → iris DIRECTLY in memory (no disk round-trip).
            print("   Converting xarray datasets to iris cubes (in-memory)...")

            # .to_iris() preserves dask lazy arrays — no materialisation until needed.
            # Normalise calendar to proleptic_gregorian as required by ISIMIP3BASD.
            _of_cube = _to_proleptic_gregorian(obs_fine[self.variable].to_iris())
            _of_cube.data = da.ma.masked_invalid(_of_cube.core_data())
            obs_fine_cube = _of_cube

            # ── Reproject to regular lat/lon if obs_fine is curvilinear ──────
            # iris.analysis.Linear (RectilinearRegridder) requires the TARGET
            # cube to have 1D DimCoords.  Rotated-pole grids (e.g. HYRAS) only
            # have 2D AuxCoords, so we reproject to a regular grid first.
            # This must happen BEFORE the crop so the crop operates on the
            # already-regular (1D DimCoord) cube whose shape is well-defined
            # in terms of lat/lon axes.
            if _is_curvilinear(obs_fine_cube):
                print("   ⚠  obs_fine has curvilinear (rotated-pole) coords — "
                      "reprojecting to regular lat/lon grid via scipy griddata...")
                obs_fine_cube = _regrid_curvilinear_to_rectilinear(obs_fine_cube)
                obs_fine_cube = _to_proleptic_gregorian(obs_fine_cube)
                print(f"   obs_fine after reprojection: {obs_fine_cube.shape}")

            # ── Crop obs_fine so its spatial shape is an exact integer multiple
            # of sim_coarse's spatial shape.  ISIMIP3BASD's get_downscaling_factors
            # asserts shape_fine % shape_coarse == 0.  Crop is done after
            # reprojection so we work with the regular 1D-coord cube.
            lat_coarse = sim_coarse['lat'].values if 'lat' in sim_coarse.coords else sim_coarse['latitude'].values
            lon_coarse = sim_coarse['lon'].values if 'lon' in sim_coarse.coords else sim_coarse['longitude'].values
            nlat_coarse, nlon_coarse = len(lat_coarse), len(lon_coarse)

            nlat_fine, nlon_fine = obs_fine_cube.shape[1], obs_fine_cube.shape[2]
            nlat_crop = (nlat_fine // nlat_coarse) * nlat_coarse
            nlon_crop = (nlon_fine // nlon_coarse) * nlon_coarse

            if nlat_crop != nlat_fine or nlon_crop != nlon_fine:
                lat_trim = (nlat_fine - nlat_crop) // 2
                lon_trim = (nlon_fine - nlon_crop) // 2
                obs_fine_cube = obs_fine_cube[
                    :,
                    lat_trim : lat_trim + nlat_crop,
                    lon_trim : lon_trim + nlon_crop,
                ]
                print(f"   ⚠  obs_fine cropped to ({nlat_crop}×{nlon_crop}) "
                      f"[was {nlat_fine}×{nlon_fine}] so shape is divisible by "
                      f"coarse ({nlat_coarse}×{nlon_coarse})")

            _sc_cube = _to_proleptic_gregorian(sim_coarse[self.variable].to_iris())
            _sc_cube.data = da.ma.masked_invalid(_sc_cube.core_data())
            sim_coarse_cube = _sc_cube

            # ── Stamp identical GeogCS on both cubes ──────────────────────────
            # iris.analysis.Linear (RectilinearRegridder) requires either:
            #   (a) both cubes have the same coordinate system, OR
            #   (b) neither has a CRS AND all coordinate metadata strings match.
            # xarray→iris conversion produces 'degrees_north'/'degrees_east'
            # while _regrid_curvilinear_to_rectilinear produces 'degrees'.
            # Stamping GeogCS on both eliminates the string-metadata comparison.
            _stamp_geogcs(obs_fine_cube)
            _stamp_geogcs(sim_coarse_cube)
            print(f"   obs_fine  {obs_fine_cube.shape}, sim_coarse {sim_coarse_cube.shape}")

            # Bilinear interpolation of coarse → fine grid.
            # obs_fine may have a rotated-pole grid (2D lat/lon on y/x dims),
            # so xr.interp cannot be used.  Use the vendor remapbil() instead:
            # it calls iris.analysis.Linear which handles any CRS, with a
            # iris.analysis.Nearest fallback for masked edge cells.
            #
            # Process in yearly chunks in parallel (ThreadPoolExecutor).
            # iris.analysis.Linear is CPU-bound but releases the GIL during
            # the scipy RegularGridInterpolator calls, so threading gives a
            # meaningful speed-up without multiprocessing complexity.
            print("   Creating bilinearly interpolated intermediate data "
                  "(parallel by year)...")
            from climdata._vendor.isimip3basd import utility_functions as uf
            import iris.cube
            from concurrent.futures import ThreadPoolExecutor, as_completed
            try:
                from tqdm import tqdm as _tqdm
            except ImportError:
                _tqdm = None

            # get year labels from the coarse cube time axis
            _time_coord = sim_coarse_cube.coord('time')
            _datetimes  = _time_coord.units.num2date(_time_coord.points)
            _years      = np.array([d.year for d in _datetimes])
            _unique_yrs = np.unique(_years)

            def _remapbil_year(yr):
                _idx = np.where(_years == yr)[0]
                # obs_fine_cube is used only as a spatial grid target by remapbil
                # (iris.analysis.Linear just needs one cube with the right grid).
                # Use a single time slice — indexing obs_fine by sim_coarse's
                # time positions would go out of bounds (obs_fine is 1951-1980,
                # sim_coarse is 2015-2100).
                return yr, uf.remapbil(sim_coarse_cube[_idx], obs_fine_cube[0])

            # n_processes workers — matches the SD parallelism already requested
            _n_remap = max(1, self.n_processes)
            _results_map = {}
            with ThreadPoolExecutor(max_workers=_n_remap) as _pool:
                _futures = {_pool.submit(_remapbil_year, yr): yr
                            for yr in _unique_yrs}
                _iter = as_completed(_futures)
                if _tqdm is not None:
                    _iter = _tqdm(_iter, total=len(_unique_yrs),
                                  desc="   remapbil", unit="yr",
                                  dynamic_ncols=True)
                for _fut in _iter:
                    _yr, _slice = _fut.result()
                    _results_map[_yr] = _slice

            _remap_slices = [_results_map[yr] for yr in _unique_yrs]
            _remap_cube = iris.cube.CubeList(_remap_slices).concatenate_cube()
            _remap_cube = _to_proleptic_gregorian(_remap_cube)
            _remap_cube.data = da.ma.masked_invalid(
                da.from_array(np.ma.asarray(_remap_cube.data))
            )
            sim_coarse_remapbil_cube = _remap_cube
            
            # Prepare output path (always use absolute path)
            if output_path is None:
                sim_fine_path = (tmpdir / "sim_fine.nc").resolve()
            else:
                sim_fine_path = Path(output_path).resolve()
            
            # Create npy_stack directory for ISIMIP3BASD output
            from climdata._vendor.isimip3basd import utility_functions as uf
            npy_stack_path = uf.npy_stack_dir(str(sim_fine_path))
            Path(npy_stack_path).mkdir(parents=True, exist_ok=True)
            
            # Setup npy_stack with proper metadata (crucial for loading results)
            uf.setup_npy_stack(str(sim_fine_path), obs_fine_cube.shape)
            print(f"   Created npy_stack directory: {npy_stack_path}")
            
            print("   Running ISIMIP3BASD statistical downscaling (year-by-year loop)...")
            
            # ── Year-by-year loop to manage memory and reduce bottlenecks ─────
            # Extract years from sim_coarse (future data)
            _time_coord_sc = sim_coarse_cube.coord('time')
            _datetimes_sc  = _time_coord_sc.units.num2date(_time_coord_sc.points)
            _years_sc      = np.array([d.year for d in _datetimes_sc])
            _unique_yrs_sd = np.unique(_years_sc)
            
            print(f"   Processing {len(_unique_yrs_sd)} years: "
                  f"{_unique_yrs_sd[0]}—{_unique_yrs_sd[-1]}")
            
            for _yr_idx, _yr in enumerate(_unique_yrs_sd, start=1):
                _idx_yr = np.where(_years_sc == _yr)[0]
                
                # Slice cubes to this year only
                _sc_yr = sim_coarse_cube[_idx_yr]
                _remap_yr = sim_coarse_remapbil_cube[_idx_yr]
                
                print(f"   [{_yr_idx}/{len(_unique_yrs_sd)}] Downscaling year {_yr} "
                      f"({_sc_yr.shape[0]} timesteps)...")
                _t0_yr = __import__('time').time()
                
                # Downscale this year's data
                sd.downscale(
                    obs_fine=obs_fine_cube,
                    sim_coarse=_sc_yr,
                    sim_coarse_remapbil=_remap_yr,
                    sim_fine_path=str(sim_fine_path),
                    n_processes=self.n_processes,
                    n_iterations=self.n_iterations,
                    randomization_seed=self.randomization_seed,
                    lower_bound=self.lower_bound,
                    lower_threshold=self.lower_threshold,
                    upper_bound=self.upper_bound,
                    upper_threshold=self.upper_threshold,
                    **kwargs,
                    **self.kwargs
                )
                
                _elapsed_yr = __import__('time').time() - _t0_yr
                print(f"      ✓ Year {_yr} complete in {_elapsed_yr:.0f}s")
            
            print("   ✅ Statistical downscaling complete!")
            print("   📦 Collecting results from npy_stack...")
            
            # Collect output from npy_stack (similar to ISIMIP3BASD's main())
            import dask.array as da
            from climdata._vendor.isimip3basd import utility_functions as uf
            
            # Start with obs_fine structure
            sim_fine_cube = obs_fine_cube.copy()
            
            # Load npy_stack and reshape to match obs_fine shape
            npy_stack_path = uf.npy_stack_dir(str(sim_fine_path))
            d = da.from_npy_stack(npy_stack_path, mmap_mode=None).reshape(sim_fine_cube.shape)
            
            # Set the collected data
            sim_fine_cube.data = np.ma.masked_array(d.compute())
            
            # Save the result cube to file
            print(f"   💾 Saving result to: {sim_fine_path}")
            iris.save(sim_fine_cube, str(sim_fine_path),
                     saver=iris.fileformats.netcdf.save,
                     unlimited_dimensions=['time'],
                     zlib=True, complevel=1)
            
            print(f"   📂 File saved: {sim_fine_path.exists()}")
            
            # Load result and convert to xarray
            result_cube = iris.load_cube(str(sim_fine_path))
            # Force load data into memory
            result_cube.data
            result = xr.DataArray.from_iris(result_cube).to_dataset(name=self.variable)
            result.load()
            result.attrs['downscaling_method'] = 'ISIMIP3BASD modified MBCn'
            result.attrs['n_iterations'] = self.n_iterations
            
            if output_path is not None:
                result.to_netcdf(output_path)
                print(f"   💾 Saved to: {output_path}")
            
        finally:
            # Clean up temporary files
            import shutil
            if tmpdir.exists():
                shutil.rmtree(tmpdir, ignore_errors=True)
        
        return result


class BCSD:
    """
    Complete Bias Correction and Statistical Downscaling workflow.
    
    This class combines bias correction and statistical downscaling in a
    single workflow, following the ISIMIP3BASD methodology.
    
    Parameters
    ----------
    variable : str
        Climate variable name (tas, pr, rsds, hurs, etc.)
    bias_correction_kwargs : dict, optional
        Arguments for BiasCorrection class
    downscaling_kwargs : dict, optional
        Arguments for StatisticalDownscaling class
    regridding_method : str, optional
        Method for regridding obs_fine to coarse resolution
    regridding_tool : str, optional
        Tool for regridding: 'xesmf' or 'cdo'
    cdo_method : str, optional
        CDO regridding method if using CDO
    cdo_env : str, optional
        Conda environment with CDO
    weights_dir : str, optional
        Directory to save regridding weights
    
    Examples
    --------
    >>> # Complete BCSD workflow for temperature
    >>> bcsd = BCSD(
    ...     variable='tas',
    ...     bias_correction_kwargs={
    ...         'distribution': 'normal',
    ...         'trend_preservation': 'additive',
    ...         'detrend': True
    ...     },
    ...     downscaling_kwargs={
    ...         'n_iterations': 20
    ...     }
    ... )
    >>> 
    >>> # Run complete workflow
    >>> result = bcsd.run(
    ...     obs_fine=obs_025deg,
    ...     sim_hist_coarse=gcm_hist_1deg,
    ...     sim_fut_coarse=gcm_fut_1deg
    ... )
    """
    
    def __init__(
        self,
        variable: str,
        bias_correction_kwargs: Optional[Dict] = None,
        downscaling_kwargs: Optional[Dict] = None,
        regridding_method: str = 'conservative',
        regridding_tool: str = 'xesmf',
        cdo_method: str = 'remapcon',
        cdo_env: str = 'cdo_stable',
        weights_dir: Optional[str] = None
    ):
        self.variable = variable
        self.regridding_method = regridding_method
        self.regridding_tool = regridding_tool
        self.cdo_method = cdo_method
        self.cdo_env = cdo_env
        self.weights_dir = weights_dir
        
        bc_kwargs = bias_correction_kwargs or {}
        sd_kwargs = downscaling_kwargs or {}
        
        self.bias_correction = BiasCorrection(variable=variable, **bc_kwargs)
        self.downscaling = StatisticalDownscaling(variable=variable, **sd_kwargs)
        
        print(f"\n{'='*60}")
        print(f"BCSD Pipeline initialized for {variable}")
        print(f"Regridding: {regridding_tool} ({regridding_method})")
        print(f"{'='*60}")
    
    def run(
        self,
        obs_fine: xr.Dataset,
        sim_hist_coarse: xr.Dataset,
        sim_fut_coarse: xr.Dataset,
        obs_hist_coarse: Optional[xr.Dataset] = None,
        output_path: Optional[Union[str, Path]] = None,
        save_intermediate: bool = False
    ) -> xr.Dataset:
        """
        Run complete BCSD workflow.
        
        Parameters
        ----------
        obs_fine : xr.Dataset
            Historical observations at fine resolution
        sim_hist_coarse : xr.Dataset
            Historical GCM simulations at coarse resolution
        sim_fut_coarse : xr.Dataset
            Future GCM simulations at coarse resolution
        obs_hist_coarse : xr.Dataset, optional
            Historical observations at coarse (GCM) resolution.
            If None, automatically derived from obs_fine by regridding.
        output_path : str or Path, optional
            Path to save final output
        save_intermediate : bool, optional
            Save intermediate bias-corrected data
        
        Returns
        -------
        xr.Dataset
            Bias-corrected and downscaled future simulations at fine resolution
        
        Workflow
        --------
        0. (Optional) Derive coarse observations from fine observations:
           - Regrid obs_fine to match sim_hist_coarse grid
        1. Bias correction at coarse resolution:
           - Correct sim_fut_coarse using obs_hist_coarse and sim_hist_coarse
        2. Statistical downscaling:
           - Downscale corrected coarse data to fine resolution using obs_fine
        """
        print(f"\n{'='*60}")
        print("Starting BCSD Workflow")
        print(f"{'='*60}\n")
        
        # Step 0: Derive coarse observations if not provided
        if obs_hist_coarse is None:
            print("STEP 0: DERIVING COARSE OBSERVATIONS FROM FINE")
            print("-" * 60)
            print("No coarse observations provided. Regridding obs_fine to coarse resolution...")
            
            obs_hist_coarse = regrid_to_coarse(
                obs_fine,
                sim_hist_coarse,
                method=self.regridding_method,
                regridding_tool=self.regridding_tool,
                cdo_method=self.cdo_method,
                cdo_env=self.cdo_env,
                weights_dir=self.weights_dir
            )
            
            if save_intermediate and output_path:
                coarse_obs_path = Path(output_path).parent / f"{Path(output_path).stem}_obs_coarse.nc"
                obs_hist_coarse.to_netcdf(coarse_obs_path)
                print(f"   💾 Saved coarse observations to: {coarse_obs_path}")
            
            print("\n")
        else:
            print("Using provided coarse observations\n")
        
        # Step 1: Bias Correction at coarse resolution
        print("STEP 1: BIAS CORRECTION")
        print("-" * 60)
        
        # Determine output path for bias-corrected data
        # ISIMIP3BASD requires a path for temporary files, so provide one even if not saving
        bc_output = None
        if output_path:
            if save_intermediate:
                bc_output = Path(output_path).parent / f"{Path(output_path).stem}_bias_corrected.nc"
        
        sim_fut_corrected = self.bias_correction.correct(
            obs_hist=obs_hist_coarse,
            sim_hist=sim_hist_coarse,
            sim_fut=sim_fut_coarse,
            output_path=bc_output
        )
        
        print("\n")
        print("STEP 2: STATISTICAL DOWNSCALING")
        print("-" * 60)
        
        # Step 2: Statistical Downscaling to fine resolution
        result = self.downscaling.downscale(
            obs_fine=obs_fine,
            sim_coarse=sim_fut_corrected,
            output_path=output_path
        )
        
        print(f"\n{'='*60}")
        print("✅ BCSD Workflow Complete!")
        print(f"{'='*60}\n")
        
        return result
