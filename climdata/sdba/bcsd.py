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
            filename=weight_file,
            reuse_weights=weight_file and os.path.exists(weight_file)
        )
        
        # Regrid each variable
        result_vars = {}
        for var in fine_data.data_vars:
            print(f"   Regridding variable: {var}")
            regridded = regridder(fine_data[var])
            # Preserve attributes
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
        variable: str,
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
        **kwargs
    ):
        self.variable = variable
        self.method = method
        self.n_processes = n_processes
        self.n_quantiles = n_quantiles
        
        # Get default config for this variable
        default_config = self.DEFAULT_CONFIGS.get(variable, {})
        
        # Use provided values or fall back to defaults
        self.distribution = distribution if distribution is not None else default_config.get('distribution')
        self.trend_preservation = trend_preservation if trend_preservation is not None else default_config.get('trend_preservation', 'additive')
        self.detrend = detrend if detrend is not False else default_config.get('detrend', False)
        self.adjust_p_values = adjust_p_values if adjust_p_values is not False else default_config.get('adjust_p_values', False)
        self.lower_bound = lower_bound if lower_bound is not None else default_config.get('lower_bound')
        self.lower_threshold = lower_threshold if lower_threshold is not None else default_config.get('lower_threshold')
        self.upper_bound = upper_bound if upper_bound is not None else default_config.get('upper_bound')
        self.upper_threshold = upper_threshold if upper_threshold is not None else default_config.get('upper_threshold')
        self.halfwin_upper_bound_climatology = halfwin_upper_bound_climatology if halfwin_upper_bound_climatology != 0 else default_config.get('halfwin_upper_bound_climatology', 0)
        
        # Store additional kwargs
        self.kwargs = kwargs
        
        print(f"🔧 BiasCorrection initialized for {variable}")
        print(f"   Distribution: {self.distribution}")
        print(f"   Trend preservation: {self.trend_preservation}")
        print(f"   Detrend: {self.detrend}")
    
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
        
        print(f"\n🔄 Starting bias correction for {self.variable}...")
        print(f"   Obs hist period: {obs_hist.time.values[0]} to {obs_hist.time.values[-1]}")
        print(f"   Sim hist period: {sim_hist.time.values[0]} to {sim_hist.time.values[-1]}")
        print(f"   Sim fut period: {sim_fut.time.values[0]} to {sim_fut.time.values[-1]}")
        
        # Create temporary directory in current working directory
        tmpdir = Path.cwd() / "tmp_bcsd"
        tmpdir.mkdir(parents=True, exist_ok=True)
        
        try:
            # Convert xarray to NetCDF, then load as iris cubes
            print("   Converting xarray datasets to iris cubes...")
            obs_hist_path = tmpdir / "obs_hist.nc"
            sim_hist_path = tmpdir / "sim_hist.nc"
            sim_fut_path = tmpdir / "sim_fut.nc"
            
            # Save as NetCDF
            obs_hist[self.variable].to_netcdf(obs_hist_path)
            sim_hist[self.variable].to_netcdf(sim_hist_path)
            sim_fut[self.variable].to_netcdf(sim_fut_path)
            
            # Load as iris cubes
            obs_hist_cube = iris.load_cube(str(obs_hist_path))
            sim_hist_cube = iris.load_cube(str(sim_hist_path))
            sim_fut_cube = iris.load_cube(str(sim_fut_path))
            
            # Prepare output path for iris cube (always use absolute path)
            if output_path is None:
                sim_fut_ba_path = (tmpdir / "sim_fut_ba.nc").resolve()
            else:
                sim_fut_ba_path = Path(output_path).resolve()
            
            print("   Running ISIMIP3BASD bias adjustment...")
            print(f"   (This may take a while for large datasets)")
            
            # Run bias adjustment with iris cubes (ISIMIP3BASD expects lists)
            ba.adjust_bias(
                obs_hist=[obs_hist_cube],
                sim_hist=[sim_hist_cube],
                sim_fut=[sim_fut_cube],
                sim_fut_ba_path=[str(sim_fut_ba_path)],
                n_processes=self.n_processes,
                n_quantiles=self.n_quantiles,
                distribution=[self.distribution] if self.distribution else [None],
                trend_preservation=[self.trend_preservation],
                detrend=[self.detrend],
                adjust_p_values=[self.adjust_p_values],
                halfwin_upper_bound_climatology=[self.halfwin_upper_bound_climatology],
                lower_bound=[self.lower_bound] if self.lower_bound is not None else [None],
                lower_threshold=[self.lower_threshold] if self.lower_threshold is not None else [None],
                upper_bound=[self.upper_bound] if self.upper_bound is not None else [None],
                upper_threshold=[self.upper_threshold] if self.upper_threshold is not None else [None],
                **kwargs,
                **self.kwargs
            )
            
            print("   ✅ Bias correction complete!")
            print("   📦 Collecting results from npy_stack...")
            
            # Collect output from npy_stack (similar to ISIMIP3BASD's main())
            import dask.array as da
            from climdata._vendor.isimip3basd import utility_functions as uf
            
            # Load npy_stack and reshape to match sim_fut shape
            npy_stack_path = uf.npy_stack_dir(str(sim_fut_ba_path))
            d = da.from_npy_stack(npy_stack_path, mmap_mode=None).reshape(sim_fut_cube.shape)
            
            # Create result cube with collected data
            sim_fut_ba_cube = sim_fut_cube.copy()
            sim_fut_ba_cube.data = np.ma.masked_array(d.compute())
            
            # Save the result cube to file
            print(f"   💾 Saving result to: {sim_fut_ba_path}")
            iris.save(sim_fut_ba_cube, str(sim_fut_ba_path),
                     saver=iris.fileformats.netcdf.save,
                     unlimited_dimensions=['time'],
                     zlib=True, complevel=1)
            
            print(f"   📂 File saved: {sim_fut_ba_path.exists()}")
            
            # Load result and convert to xarray
            result_cube = iris.load_cube(str(sim_fut_ba_path))
            # Force load data into memory to avoid lazy loading issues after cleanup
            result_cube.data  # This triggers the lazy data to be loaded
            result = xr.DataArray.from_iris(result_cube).to_dataset(name=self.variable)
            # Ensure data is loaded (not lazy)
            result.load()
            result.attrs['bias_correction_method'] = 'ISIMIP3BASD'
            result.attrs['distribution'] = str(self.distribution)
            result.attrs['trend_preservation'] = self.trend_preservation
            
            if output_path is not None:
                result.to_netcdf(output_path)
                print(f"   💾 Saved to: {output_path}")
            
        finally:
            # Clean up temporary files
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
        
        # Get default config
        default_config = self.DEFAULT_CONFIGS.get(variable, {})
        self.lower_bound = lower_bound if lower_bound is not None else default_config.get('lower_bound')
        self.lower_threshold = lower_threshold if lower_threshold is not None else default_config.get('lower_threshold')
        self.upper_bound = upper_bound if upper_bound is not None else default_config.get('upper_bound')
        self.upper_threshold = upper_threshold if upper_threshold is not None else default_config.get('upper_threshold')
        
        self.kwargs = kwargs
        
        print(f"🔧 StatisticalDownscaling initialized for {variable}")
        print(f"   Iterations: {n_iterations}")
    
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
            # Convert xarray to NetCDF, then load as iris cubes
            print("   Converting xarray datasets to iris cubes...")
            obs_fine_path = tmpdir / "obs_fine.nc"
            sim_coarse_path = tmpdir / "sim_coarse.nc"
            
            obs_fine[self.variable].to_netcdf(obs_fine_path)
            sim_coarse[self.variable].to_netcdf(sim_coarse_path)
            
            # Load as iris cubes
            obs_fine_cube = iris.load_cube(str(obs_fine_path))
            sim_coarse_cube = iris.load_cube(str(sim_coarse_path))
            
            # Create remapped coarse data (bilinear interpolation)
            print("   Creating bilinearly interpolated intermediate data...")
            sim_coarse_remapbil = sim_coarse[self.variable].interp(
                lat=obs_fine.lat,
                lon=obs_fine.lon,
                method='linear'
            )
            sim_coarse_remapbil_path = tmpdir / "sim_coarse_remapbil.nc"
            sim_coarse_remapbil.to_netcdf(sim_coarse_remapbil_path)
            sim_coarse_remapbil_cube = iris.load_cube(str(sim_coarse_remapbil_path))
            
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
            
            print("   Running ISIMIP3BASD statistical downscaling...")
            print(f"   (This may take a while for large datasets)")
            
            # Run downscaling with iris cubes (single cubes, not lists)
            sd.downscale(
                obs_fine=obs_fine_cube,
                sim_coarse=sim_coarse_cube,
                sim_coarse_remapbil=sim_coarse_remapbil_cube,
                sim_fine_path=str(sim_fine_path),
                n_processes=self.n_processes,
                n_iterations=self.n_iterations,
                lower_bound=self.lower_bound,
                lower_threshold=self.lower_threshold,
                upper_bound=self.upper_bound,
                upper_threshold=self.upper_threshold,
                **kwargs,
                **self.kwargs
            )
            
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
