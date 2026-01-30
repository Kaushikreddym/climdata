# (C) 2019 Potsdam Institute for Climate Impact Research (PIK)
# 
# This file is part of ISIMIP3BASD.
#
# ISIMIP3BASD is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ISIMIP3BASD is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with ISIMIP3BASD. If not, see <http://www.gnu.org/licenses/>.



"""
Statistical downscaling
=======================

Provides functions for statistical downscaling of climate simulation data using
climate observation data with the same temporal and higher spatial resolution.

The following variable-specific parameter values (variable units in brackets)
were used to produce the results presented in
Lange (2019) <https://doi.org/10.5194/gmd-12-3055-2019>.

hurs (%)
    --lower-bound 0
    --lower-threshold .01
    --upper-bound 100
    --upper-threshold 99.99

pr (mm day-1)
    --lower-bound 0
    --lower-threshold .1

prsnratio (1)
    --lower-bound 0
    --lower-threshold .0001
    --upper-bound 1
    --upper-threshold .9999
    --if-all-invalid-use 0.

psl (Pa)

rlds (W m-2)

rsds (W m-2)
    --lower-bound 0
    --lower-threshold .01

sfcWind (m s-1)
    --lower-bound 0
    --lower-threshold .01

tas (K)

tasrange (K)
    --lower-bound 0
    --lower-threshold .01

tasskew (1)
    --lower-bound 0
    --lower-threshold .0001
    --upper-bound 1
    --upper-threshold .9999

"""



import os
import sys
import dask
import iris
import shutil
import warnings
import numpy as np
import dask.array as da
from . import utility_functions as uf
import iris.coord_categorisation as icc
import multiprocessing as mp
from optparse import OptionParser
from functools import partial



# shared resources
global_lazy_data = {}
global_month_numbers = {}



def initializer(l, m):
    """
    Sets global variables global_lazy_data and global_month_numbers to l and m,
    respectively. These global variables are used as shared resources by all
    processes in the pool for n_processes > 1. Using these shared resources
    drastically reduces the amount of data that needs to be piped to the
    processes at the beginning of every local statistical downscaling.

    Parameters
    ----------
    l : dict of dask arrays
        Dictionary keys are 'obs_fine', 'sim_coarse', 'sim_coarse_remapbil'. 
        Every array represents one climate dataset.
    m : dict of arrays
        Dictionary keys are 'obs_fine', 'sim_coarse', 'sim_coarse_remapbil'.
        Every array represents a month-number time series.

    """
    global global_lazy_data
    global global_month_numbers
    global_lazy_data = l
    global_month_numbers = m



def weighted_sum_preserving_mbcn(
        x_obs, x_sim_coarse, x_sim,
        sum_weights, rotation_matrices=[], n_quantiles=50):
    """
    Applies the core of the modified MBCn algorithm for statistical downscaling
    as described in Lange (2019) <https://doi.org/10.5194/gmd-12-3055-2019>.

    Parameters
    ----------
    x_obs : (M,N) ndarray
        Array of N observed time series of M time steps each at fine spatial
        resolution.
    x_sim_coarse : (M,) array
        Array of simulated time series of M time steps at coarse spatial
        resolution.
    x_sim : (M,N) ndarray
        Array of N simulated time series of M time steps each at fine spatial
        resolution, derived from x_sim_coarse by bilinear interpolation.
    sum_weights : (N,) array
        Array of N grid cell-area weights.
    rotation_matrices : list of (N,N) ndarrays, optional
        List of orthogonal matrices defining a sequence of rotations in the  
        second dimension of x_obs and x_sim.
    n_quantiles : int, optional
        Number of quantile-quantile pairs used for non-parametric quantile
        mapping.

    Returns
    -------
    x_sim : (M,N) ndarray
        Result of application of the modified MBCn algorithm.

    """
    # initialize total rotation matrix
    n_variables = sum_weights.size
    o_total = np.diag(np.ones(n_variables))

    # p-values in percent for non-parametric quantile mapping
    p = np.linspace(0., 1., n_quantiles+1)

    # normalise the sum weights vector to length 1
    sum_weights = sum_weights / np.sqrt(np.sum(np.square(sum_weights)))

    # rescale x_sim_coarse for initial step of algorithm
    x_sim_coarse = x_sim_coarse * np.sum(sum_weights)

    # iterate
    n_loops = len(rotation_matrices) + 2
    for i in range(n_loops):
        if not i:  # rotate to the sum axis
            o = uf.generate_rotation_matrix_fixed_first_axis(sum_weights)
        elif i == n_loops - 1:  # rotate back to original axes for last qm
            o = o_total.T
        else:  # do random rotation
            o = rotation_matrices[i-1]

        # compute total rotation
        o_total = np.dot(o_total, o)

        # rotate data
        x_sim = np.dot(x_sim, o)
        x_obs = np.dot(x_obs, o)
        sum_weights = np.dot(sum_weights, o)

        if not i:
            # restore simulated values at coarse grid scale
            x_sim[:,0] = x_sim_coarse

            # quantile map observations to values at coarse grid scale
            q_sim = uf.percentile1d(x_sim_coarse, p)
            q_obs = uf.percentile1d(x_obs[:,0], p)
            x_obs[:,0] = \
                uf.map_quantiles_non_parametric_with_constant_extrapolation(
                x_obs[:,0], q_obs, q_sim)
        else:
            # do univariate non-parametric quantile mapping for every variable
            x_sim_previous = x_sim.copy()
            for j in range(n_variables):
                q_sim = uf.percentile1d(x_sim[:,j], p)
                q_obs = uf.percentile1d(x_obs[:,j], p)
                x_sim[:,j] = \
                    uf.map_quantiles_non_parametric_with_constant_extrapolation(
                    x_sim[:,j], q_sim, q_obs)

            # preserve weighted sum of original variables
            if i < n_loops - 1:
                x_sim -= np.outer(np.dot(
                   x_sim - x_sim_previous, sum_weights), sum_weights)

    return x_sim



def downscale_one_month(
        data, long_term_mean,
        lower_bound=None, lower_threshold=None,
        upper_bound=None, upper_threshold=None,
        randomization_seed=None, **kwargs):
    """
    1. Replaces invalid values in time series.
    2. Replaces values beyond thresholds by random numbers.
    3. Applies the modified MBCn algorithm for statistical downscaling.
    4. Replaces values beyond thresholds by the respective bound.

    Parameters
    ----------
    data : dict of masked arrays
        The arrays are of shape (M,N), (M,), (M,N) for key 'x_obs_fine', 
        'x_sim_coarse', 'x_sim_coarse_remapbil', respectively.
    long_term_mean : dict
        For keys 'x_obs_fine' and 'x_sim_coarse_remapbil', an array is expected.
        For key 'x_sim_coarse', a scalar is expected.
        Values represents the average of all valid values in the complete time
        series for one climate variable and one location.
    lower_bound : float, optional
        Lower bound of values in data.
    lower_threshold : float, optional
        Lower threshold of values in data. All values below this threshold are
        replaced by random numbers between lower_bound and lower_threshold
        before application of the modified MBCn algorithm.
    upper_bound : float, optional
        Upper bound of values in data.
    upper_threshold : float, optional
        Upper threshold of values in data. All values above this threshold are
        replaced by random numbers between upper_threshold and upper_bound
        before application of the modified MBCn algorithm.
    randomization_seed : int, optional
        Used to seed the random number generator before replacing values beyond
        the specified thresholds.

    Returns
    -------
    x_sim_fine : (M,N) ndarray
        Result of application of the modified MBCn algorithm.

    Other Parameters
    ----------------
    **kwargs : Passed on to weighted_sum_preserving_mbcn.
    
    """
    x = {}
    for key, d in data.items():
        # remove invalid values from masked array and store resulting data array
        x[key] = uf.sample_invalid_values(
            d, randomization_seed, long_term_mean[key])[0]

        # randomize censored values, use high powers to create many values close
        # to the bounds as this keeps weighted sums similar to original values
        x[key] = uf.randomize_censored_values(x[key], 
            lower_bound, lower_threshold, upper_bound, upper_threshold,
            False, False, randomization_seed, 10., 10.)

    # downscale
    x_sim_coarse_remapbil = x['sim_coarse_remapbil'].copy()
    x_sim_fine = weighted_sum_preserving_mbcn(
        x['obs_fine'], x['sim_coarse'], x['sim_coarse_remapbil'], **kwargs)

    # de-randomize censored values
    uf.randomize_censored_values(x_sim_fine, 
        lower_bound, lower_threshold, upper_bound, upper_threshold, True, True)

    # make sure there are no invalid values
    uf.assert_no_infs_or_nans(x_sim_coarse_remapbil, x_sim_fine)

    return x_sim_fine



def downscale_one_location(
        i_loc_coarse, sim_fine_path, downscaling_factors, sum_weights,
        months=[1,2,3,4,5,6,7,8,9,10,11,12], fill_value=1.e20,
        lower_bound=None, lower_threshold=None,
        upper_bound=None, upper_threshold=None,
        if_all_invalid_use=np.nan, resume_job=False, **kwargs):
    """
    Applies the modified MBCn algorithm for statistical downscaling calendar
    month by calendar month to climate data within one coarse grid cell.

    Parameters
    ----------
    i_loc_coarse : tuple
        Coarse location index.
    sim_fine_path : str
        Path used to store result of statistical downscaling.
    downscaling_factors : array of ints
        Downscaling factors for all grid dimensions.
    sum_weights : ndarray
        Array of fine grid cell area weights.
    months : list, optional
        List of ints from {1,...,12} representing calendar months for which 
        results of statistical downscaling are to be returned.
    fill_value : float, optional
        Value used to indicate missing values.
    lower_bound : float, optional
        Lower bound of values in data.
    lower_threshold : float, optional
        Lower threshold of values in data.
    upper_bound : float, optional
        Upper bound of values in data.
    upper_threshold : float, optional
        Upper threshold of values in data.
    if_all_invalid_use : float, optional
        Used to replace invalid values if there are no valid values.
    resume_job : boolean, optional
        Abort if results for this location already exist.

    Returns
    -------
    None.

    Other Parameters
    ----------------
    **kwargs : Passed on to downscale_one_month.

    """
    # abort here if results for this location already exist
    i_locs = uf.get_fine_location_indices(i_loc_coarse, downscaling_factors)
    i_loc_1d = lambda i : np.ravel_multi_index(i, sum_weights.shape)
    i_loc_path = lambda p, i : uf.npy_stack_dir(p) + '%i.npy'%(i_loc_1d(i))
    if resume_job:
        i_loc_done = True
        for i, i_loc in enumerate(i_locs):
            if not os.path.isfile(i_loc_path(sim_fine_path, i_loc)):
                i_loc_done = False
        if i_loc_done:
            print(i_loc_coarse, 'done already')
            sys.stdout.flush()
            return None

    # prevent dask from opening new threads every time lazy data are realized
    # as this results in RuntimeError: can't start new thread
    # see <http://docs.dask.org/en/latest/scheduler-overview.html>
    dask.config.set(scheduler='single-threaded')

    # use shared resources
    lazy_data = global_lazy_data
    month_numbers = global_month_numbers

    # realize local lazy data
    data = {}
    i_loc_fine = tuple([slice(df * i_loc_coarse[i], df * (i_loc_coarse[i] + 1))
        for i, df in enumerate(downscaling_factors)])
    for key, d in lazy_data.items():
        if key == 'sim_coarse':
            data[key] = d[(slice(None, None),) + i_loc_coarse].compute()
        else:
            data[key] = uf.flatten_all_dimensions_but_first(
                        d[(slice(None, None),) + i_loc_fine].compute())

    # abort here if there are only missing values in at least one time series
    # do not abort though if the if_all_invalid_use option has been specified
    if np.isnan(if_all_invalid_use):
        if uf.only_missing_values_in_at_least_one_time_series(data):
            print(i_loc_coarse, 'skipped due to missing data')
            sys.stdout.flush()
            n_times = month_numbers['sim_coarse'].size
            for i, i_loc in enumerate(i_locs):
                np.save(i_loc_path(sim_fine_path, i_loc), np.expand_dims(
                    np.repeat(np.float32(fill_value), n_times), axis=1))
            return None

    # otherwise continue
    print(i_loc_coarse)
    sys.stdout.flush()
    sim_fine = data['sim_coarse_remapbil'].copy()

    # compute mean value over all time steps for invalid value sampling
    long_term_mean = {}
    for key, d in data.items():
        long_term_mean[key] = uf.average_valid_values(d, if_all_invalid_use,
            lower_bound, lower_threshold, upper_bound, upper_threshold)

    # do statistical downscaling calendar month by calendar month
    sum_weights_loc = sum_weights[i_loc_fine].flatten()
    data_this_month = {}
    for month in months:
        # extract data
        for key, d in data.items():
            m = month_numbers[key] == month
            assert np.any(m), f'no data found for month {month} in {key}'
            data_this_month[key] = d[m]

        # do statistical downscaling and store result as a masked array
        sim_fine_this_month = np.ma.array(downscale_one_month(
            data_this_month, long_term_mean,
            lower_bound, lower_threshold, upper_bound, upper_threshold,
            sum_weights=sum_weights_loc, **kwargs), fill_value=fill_value)
    
        # put downscaled data into sim_fine
        m = month_numbers['sim_coarse_remapbil'] == month
        sim_fine[m] = sim_fine_this_month

    # save local result of statistical downscaling
    for i, i_loc in enumerate(i_locs):
        np.save(i_loc_path(sim_fine_path, i_loc), np.expand_dims(
            sim_fine[:,i].data, axis=1))

    return None



def downscale(
        obs_fine, sim_coarse, sim_coarse_remapbil, 
        sim_fine_path, n_processes=1, n_iterations=20,
        randomization_seed=None, **kwargs):
    """
    Applies the modified MBCn algorithm for statistical downscaling calendar
    month by calendar month and coarse grid cell by coarse grid cell.

    Parameters
    ----------
    obs_fine : iris cube
        Cube of observed climate data at fine spatial resolution.
    sim_coarse : iris cube
        Cube of simulated climate data at coarse spatial resolution.
    sim_coarse_remapbil : iris cube
        Cube of simulated climate data at coarse spatial resolution bilinearly
        interpolated to fine spatial resolution.
    sim_fine_path : str
        Path used to store result of statistical downscaling.
    n_processes : int, optional
        Number of processes used for parallel processing.
    n_iterations : int, optional
        Number of iterations used in the modified MBCn algorithm. If not
        specified, then it is set to 2^(3+n/3), where n is the number of fine
        grid cells per coarse grid cell.
    randomization_seed : int, optional
        Used to seed the random number generator before generating random 
        rotation matrices for the modified MBCn algorithm.

    Returns
    -------
    sim_fine : iris cube
        Result of application of the modified MBCn algorithm.

    Other Parameters
    ----------------
    **kwargs : Passed on to downscale_one_location.

    """
    # prepare statistical downscaling location by location and month by month
    cubes = {
    'obs_fine': obs_fine,
    'sim_coarse': sim_coarse,
    'sim_coarse_remapbil': sim_coarse_remapbil,
    }
    lazy_data = {}
    month_numbers = {}
    for key, cube in cubes.items():
        lazy_data[key] = cube.core_data()
        time_coord = cube.coord('time')
        datetimes = time_coord.units.num2date(time_coord.points)
        month_numbers[key] = uf.convert_datetimes(datetimes, 'month_number')

    # get list of rotation matrices to be used for all locations and months
    if randomization_seed is not None: np.random.seed(randomization_seed)
    downscaling_factors = uf.get_downscaling_factors(
        obs_fine.shape[1:], sim_coarse.shape[1:])
    rotation_matrices = [uf.generateCREmatrix(np.prod(downscaling_factors))
                         for i in range(n_iterations)]

    # compute sum weights assuming a regular latitude-longitude grid
    sum_weights = iris.analysis.cartography.cosine_latitude_weights(obs_fine[0])

    # downscale every location individually
    i_locations_coarse = np.ndindex(sim_coarse.shape[1:])
    sdol = partial(downscale_one_location, 
        sim_fine_path=sim_fine_path,
        downscaling_factors=downscaling_factors,
        sum_weights=sum_weights,
        rotation_matrices=rotation_matrices,
        randomization_seed=randomization_seed, **kwargs)
    print('downscaling at coarse location ...')
    if n_processes > 1:
        pool = mp.Pool(n_processes, initializer=initializer, initargs=(
            lazy_data, month_numbers))
        foo = list(pool.imap(sdol, i_locations_coarse))
        pool.close()
        pool.join()
        pool.terminate()
    else:
        initializer(lazy_data, month_numbers)
        foo = list(map(sdol, i_locations_coarse))



def main():
    """
    Prepares and concludes the application of the modified MBCn algorithm for
    statistical downscaling.

    """
    # parse command line options and arguments
    parser = OptionParser()
    parser.add_option('-o', '--obs-fine', action='store',
        type='string', dest='obs_fine', default=None,
        help='path to input netcdf file with observation at fine resolution')
    parser.add_option('-s', '--sim-coarse', action='store',
        type='string', dest='sim_coarse', default=None,
        help='path to input netcdf file with simulation at coarse resolution')
    parser.add_option('-f', '--sim-fine', action='store',
        type='string', dest='sim_fine', default=None,
        help=('path to output netcdf file with simulation statistically '
              'downscaled to fine resolution'))
    parser.add_option('-v', '--variable', action='store',
        type='string', dest='variable', default=None,
        help=('standard name of variable to be downscaled in netcdf files '
              '(has to be the same in all files)'))
    parser.add_option('-m', '--months', action='store',
        type='string', dest='months', default='1,2,3,4,5,6,7,8,9,10,11,12',
        help=('comma-separated list of integers from {1,...,12} representing '
              'calendar months that shall be statistically downscaled'))
    parser.add_option('--n-processes', action='store',
        type='int', dest='n_processes', default=1,
        help='number of processes used for multiprocessing (default: 1)')
    parser.add_option('--n-iterations', action='store',
        type='int', dest='n_iterations', default=20,
        help=('number of iterations used for statistical downscaling (default: '
              '20)'))
    parser.add_option('--o-time-range', action='store',
        type='string', dest='obs_fine_tr', default=None,
        help=('time constraint for data extraction from input netcdf file with '
              'observation of format %Y%m%dT%H%M%S-%Y%m%dT%H%M%S '
              '(if not specified then no time constraint is applied)'))
    parser.add_option('--s-time-range', action='store',
        type='string', dest='sim_coarse_tr', default=None,
        help=('time constraint for data extraction from input netcdf file with '
              'simulation of format %Y%m%dT%H%M%S-%Y%m%dT%H%M%S '
              '(if not specified then no time constraint is applied)'))
    parser.add_option('--lower-bound', action='store',
        type='float', dest='lower_bound', default=None,
        help=('lower bound of variable that has to be respected during '
              'statistical downscaling (default: not specified)'))
    parser.add_option('--lower-threshold', action='store',
        type='float', dest='lower_threshold', default=None,
        help=('lower threshold of variable that has to be respected during '
              'statistical downscaling (default: not specified)'))
    parser.add_option('--upper-bound', action='store',
        type='float', dest='upper_bound', default=None,
        help=('upper bound of variable that has to be respected during '
              'statistical downscaling (default: not specified)'))
    parser.add_option('--upper-threshold', action='store',
        type='float', dest='upper_threshold', default=None,
        help=('upper threshold of variable that has to be respected during '
              'statistical downscaling (default: not specified)'))
    parser.add_option('--randomization-seed', action='store',
        type='int', dest='randomization_seed', default=None,
        help=('seed used during randomization to generate reproducible results '
              '(default: not specified)'))
    parser.add_option('-q', '--n-quantiles', action='store',
        type='int', dest='n_quantiles', default=50,
        help=('number of quantiles used for non-parametric quantile mapping '
              '(default: 50)'))
    parser.add_option('--if-all-invalid-use', action='store',
        type='float', dest='if_all_invalid_use', default=np.nan,
        help=('replace missing values, infs and nans by this value before '
              'statistical downscaling if there are no other values available '
              'in a time series (default: not specified)'))
    parser.add_option('--fill-value', action='store',
        type='float', dest='fill_value', default=1.e20,
        help=('fill value used for missing values in output netcdf file '
              '(default: 1.e20)'))
    parser.add_option('--repeat-warnings', action='store_true',
        dest='repeat_warnings', default=False,
        help='repeat warnings for the same source location (default: do not)')
    parser.add_option('--limit-time-dimension', action='store_true',
        dest='limit_time_dimension', default=False,
        help=('save output netcdf file with a limited time dimension; data '
              'array in output netcdf file are chunked in space (default: '
              'save output netcdf file with an unlimited time dimension; data '
              'array in output netcdf file are chunked in time)'))
    parser.add_option('--keep-npy-stack', action='store_true',
        dest='keep_npy_stack', default=False,
        help=('local results of statistical downscaling are stored in a stack '
              'of npy files before these are collected and saved in one netcdf '
              'file; this flag prevents the eventual removal of the npy stack '
              '(default: remove npy stack)'))
    parser.add_option('--resume-job', action='store_true',
        dest='resume_job', default=False,
        help=('use local statistical downscaling results where they exist '
              '(default: overwrite existing local statistical downscaling '
              'results)'))
    (options, args) = parser.parse_args()
    if options.repeat_warnings: warnings.simplefilter('always', UserWarning)

    # do some preliminary checks
    assert options.n_iterations > 0, 'invalid number of iterations'
    months = list(np.sort(np.unique(np.array(
        options.months.split(','), dtype=int))))
    uf.assert_validity_of_months(months)
    uf.assert_consistency_of_bounds_and_thresholds(
        options.lower_bound, options.lower_threshold,
        options.upper_bound, options.upper_threshold)

    # load input data
    print('loading input ...')
    obs_fine = uf.load_cube(
        options.obs_fine, options.variable, options.obs_fine_tr)
    sim_coarse = uf.load_cube(
        options.sim_coarse, options.variable, options.sim_coarse_tr)

    # make sure the proleptic gregorian calendar is used in all input cubes
    uf.assert_calendar(obs_fine, 'proleptic_gregorian')
    uf.assert_calendar(sim_coarse, 'proleptic_gregorian')

    # make sure that time is the leading coordinate
    uf.assert_coord_axis(obs_fine, 'time', 0)
    uf.assert_coord_axis(sim_coarse, 'time', 0)

    # bilinearly interpolate sim_coarse to grid of obs_fine
    print('interpolating to fine grid ...')
    sim_coarse_remapbil = uf.remapbil(sim_coarse, obs_fine)

    # prepare loading local bias adjustment results using da.from_npy_stack
    uf.setup_npy_stack(options.sim_fine, sim_coarse_remapbil.shape)

    # turn realized data back into lazy data to reduce memory need
    chunksizes = uf.output_chunksizes(sim_coarse_remapbil.shape)
    iris.save(sim_coarse_remapbil, options.sim_fine,
        saver=iris.fileformats.netcdf.save, fill_value=options.fill_value,
        zlib=True, complevel=1, chunksizes=chunksizes)
    del sim_coarse_remapbil
    sim_coarse_remapbil = uf.load_cube(
        options.sim_fine, options.variable, None)
    del sim_coarse
    sim_coarse = uf.load_cube(
        options.sim_coarse, options.variable, options.sim_coarse_tr)

    # do statistical downscaling
    downscale(
        obs_fine, sim_coarse, sim_coarse_remapbil,
        options.sim_fine, options.n_processes,
        options.n_iterations, options.randomization_seed,
        fill_value=options.fill_value,
        months=months,
        lower_bound=options.lower_bound,
        lower_threshold=options.lower_threshold,
        upper_bound=options.upper_bound,
        upper_threshold=options.upper_threshold,
        n_quantiles=options.n_quantiles,
        if_all_invalid_use=options.if_all_invalid_use,
        resume_job=options.resume_job)

    # collect local results of statistical downscaling
    print('saving output ...')
    sim_fine = sim_coarse_remapbil
    npy_stack_dir = uf.npy_stack_dir(options.sim_fine)
    d = da.from_npy_stack(npy_stack_dir, mmap_mode=None).reshape(sim_fine.shape)
    sim_fine.data = np.ma.masked_array(d, fill_value=options.fill_value)

    # write statistical downscaling parameters into attributes of sim_fine
    uf.add_basd_attributes(sim_fine, options, 'sd_')

    # save output data
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        iris.save(sim_fine, options.sim_fine, 
            saver=iris.fileformats.netcdf.save,
            unlimited_dimensions=None
            if options.limit_time_dimension else ['time'], 
            fill_value=options.fill_value, zlib=True, complevel=1,
            chunksizes=chunksizes if options.limit_time_dimension else None)

    # remove local results of statistical downscaling
    if not options.keep_npy_stack:
        shutil.rmtree(npy_stack_dir)



if __name__ == '__main__':
    main()
