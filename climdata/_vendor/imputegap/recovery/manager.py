import datetime
import os
import platform
import time
import numpy as np
import matplotlib
import importlib.resources
from ..tools import utils

import matplotlib.pyplot as plt


def select_backend():
    system = platform.system()
    #headless = os.getenv("DISPLAY") is None or os.getenv("CI") is not None
    if system == "Darwin":
        for backend in ["MacOSX", "Qt5Agg", "TkAgg"]:
            try:
                matplotlib.use(backend)
                return
            except (ImportError, RuntimeError):
                continue
        try:
            matplotlib.use("TkAgg")  # fallback
        except (ImportError, RuntimeError):
            matplotlib.use("Agg")

    # Linux or Windows
    else:
        for backend in ["TkAgg", "QtAgg", "Qt5Agg", "Agg"]:
            try:
                matplotlib.use(backend)
                return
            except (ImportError, RuntimeError):
                continue



class TimeSeries:
    """
    Class for managing and manipulating time series data.

    This class allows importing, normalizing, and visualizing time series datasets. It also provides methods
    to contaminate the datasets with missing values and plot results.

    Methods
    -------
    __init__() :
        Initializes the TimeSeries object.

    import_matrix(data=None) :
        Imports a matrix of time series data.

    load_series(data=None, max_series=None, max_values=None, header=False) :
        Loads time series data from a file or predefined dataset.

    print(limit=10, view_by_series=False) :
        Prints a limited number of time series from the dataset.

    print_results(metrics, algorithm="") :
        Prints the results of the imputation process.

    normalize(normalizer="z_score") :
        Normalizes the time series dataset.

    plot(input_data, incomp_data=None, recov_data=None, max_series=None, max_values=None, size=(16, 8), save_path="", display=True) :
        Plots the time series data, including raw, contaminated, or imputed data.

    Contamination :
        Class containing methods to contaminate time series data with missing values based on different patterns.

    """

    def __init__(self):
        """
        Initialize the TimeSeries object.

        The class works with time series datasets, where each series is separated by space, and values
        are separated by newline characters.

        IMPORT FORMAT : (Values,Series) : series are seperated by "SPACE" et values by "\\n"
        """
        self.data = None
        self.name = "default"
        self.plots = None
        self.algorithms = utils.list_of_algorithms()
        self.patterns = utils.list_of_patterns()
        self.datasets = utils.list_of_datasets()
        self.optimizers = utils.list_of_optimizers()
        self.extractors = utils.list_of_extractors()
        self.forecasting_models = utils.list_of_downstreams()
        self.families = utils.list_of_families()
        self.algorithms_with_families = utils.list_of_algorithms_with_families()
        select_backend()

    def import_matrix(self, data=None):
        """
        Imports a matrix of time series data.

        The data can be provided as a list or a NumPy array. The format is (Series, Values),
        where series are separated by space, and values are separated by newline characters.

        Parameters
        ----------
        data : list or numpy.ndarray, optional
            The matrix of time series data to import.

        Returns
        -------
        TimeSeries
            The TimeSeries object with the imported data.
        """
        if data is not None:
            if isinstance(data, list):
                self.data = np.array(data)

            elif isinstance(data, np.ndarray):
                self.data = data
            else:
                print("\nThe time series have not been loaded, format unknown\n")
                self.data = None
                raise ValueError("Invalid input for import_matrix")

            return self

    def load_series(self, data, nbr_series=None, nbr_val=None, header=False, replace_nan=False, verbose=True):
        """
        Loads time series data from a file or predefined dataset.

        The data is loaded as a matrix of shape (Values, Series). You can limit the number of series
        or values per series for computational efficiency.

        Parameters
        ----------
        data : str
            The file path or name of a predefined dataset (e.g., 'bafu.txt').
        nbr_series : int, optional
            The maximum number of series to load.
        nbr_val : int, optional
            The maximum number of values per series.
        header : bool, optional
            Whether the dataset has a header. Default is False.
        replace_nan : bool, optional
            The Dataset has already NaN values that needs to be replaced by 0 values.
        verbose : bool, optional
            Display information print (default: True).


        Returns
        -------
        TimeSeries
            The TimeSeries object with the loaded data.

        Example
        -------
            >>> ts.load_series(utils.search_path("eeg-alcohol"), nbr_series=50, nbr_val=100)

        """

        if data is not None:
            if isinstance(data, str):
                saved_data = data

                #  update path form inner library datasets
                if data in utils.list_of_datasets(txt=True):
                    self.name = data[:-4]
                    data = importlib.resources.files('..datasets').joinpath(data)

                if not os.path.exists(data):
                    data = ".." + saved_data
                    if not os.path.exists(data):
                        data = data[1:]

                self.data = np.genfromtxt(data, delimiter=' ', max_rows=nbr_val, skip_header=int(header))

                if verbose:
                    print("\n(SYS) The dataset is loaded from " + str(data) + "\n")

                if nbr_series is not None:
                    self.data = self.data[:, :nbr_series]
            else:
                print("\nThe dataset has not been loaded, format unknown\n")
                self.data = None
                raise ValueError("Invalid input for load_series")

            if replace_nan:
                print("\nThe NaN values has been set to zero...\n")
                self.data = np.nan_to_num(self.data)  # Replace NaNs with 0

            self.data = self.data.T

            return self

    def print(self, nbr_val=10, nbr_series=7, view_by_series=False):
        """
        Prints a limited number of time series from the dataset.

        Parameters
        ----------
        nbr_val : int, optional
        The number of timestamps to print. Default is 15. Use -1 for no restriction.
        nbr_series : int, optional
        The number of series to print. Default is 10. Use -1 for no restriction.
        view_by_series : bool, optional
        Whether to view by series (True) or by values (False).

        Returns
        -------
        None
        """
        to_print = self.data
        nbr_tot_series, nbr_tot_values = to_print.shape
        print_col, print_row = "idx", "TS"
        print_col_inc, print_row_inc = 0, 1


        print(f"\nshape of {self.name} : {self.data.shape}\n\tnumber of series = { nbr_tot_series}\n\tnumber of values = {nbr_tot_values}\n")

        if nbr_val == -1:
            nbr_val = to_print.shape[1]
        if nbr_series == -1:
            nbr_series = to_print.shape[0]
        to_print = to_print[:nbr_series, :nbr_val]

        if not view_by_series:
            to_print = to_print.T
            print_col, print_row = "TS", "idx"
            print_col_inc, print_row_inc = 1, 0

        header_format = "{:<15}"  # Fixed size for headers
        value_format = "{:>15.10f}"  # Fixed size for values
        # Print the header
        print(f"{'':<18}", end="")  # Empty space for the row labels
        for i in range(to_print.shape[1]):
            print(header_format.format(f"{print_col}_{i + print_col_inc}"), end="")
        print()

        # Print each limited series with fixed size
        for i, series in enumerate(to_print):
            print(header_format.format(f"{print_row}_{i + print_row_inc}"), end="")
            print("".join([value_format.format(elem) for elem in series]))

        if nbr_series < nbr_tot_series:
            print("...")

    def print_results(self, metrics, algorithm="", text="Results"):
        """
        Prints the results of the imputation process.

        Parameters
        ----------
        metrics : dict
           A dictionary containing the imputation metrics to display.
        algorithm : str, optional
           The name of the algorithm used for imputation.
        algorithm : str, optional
           Output text to help the user.

        Returns
        -------
        None

        Example
        -------
            >>> ts.print_results(imputer.metrics, imputer.algorithm)
        """

        if algorithm != "":
            print(f"\n{text} ({algorithm}) :")
        else:
            print(f"\n{text} :")

        for key, value in metrics.items():
            print(f"{key:<20} = {value}")

    def normalize(self, normalizer="z_score", verbose=True):
        """
        Normalize the time series dataset.

        Supported normalization techniques are "z_score" and "min_max". The method also logs
        the execution time for the normalization process.

        Parameters
        ----------
        normalizer : str, optional
            The normalization technique to use. Options are "z_score" or "min_max". Default is "z_score".
        verbose : bool, optional
        Whether to display the contamination information (default is False).

        Returns
        -------
        numpy.ndarray
            The normalized time series data.

        Example
        -------
            >>> ts.normalize(normalizer="z_score")
        """
        self.data = self.data.T

        if normalizer == "min_max":
            start_time = time.time()  # Record start time

            # Compute the min and max for each series (column-wise), ignoring NaN
            ts_min = np.nanmin(self.data, axis=0)
            ts_max = np.nanmax(self.data, axis=0)

            # Compute the range for each series, and handle cases where the range is 0
            range_ts = ts_max - ts_min
            range_ts[range_ts == 0] = 1  # Prevent division by zero for constant series

            # Apply min-max normalization
            self.data = (self.data - ts_min) / range_ts

            end_time = time.time()
        elif normalizer == "z_lib":
            from scipy.stats import zscore

            start_time = time.time()  # Record start time

            self.data = zscore(self.data, axis=0)

            end_time = time.time()

        elif normalizer == "m_lib":
            from sklearn.preprocessing import MinMaxScaler

            start_time = time.time()  # Record start time

            scaler = MinMaxScaler()
            self.data = scaler.fit_transform(self.data)

            end_time = time.time()
        else:
            start_time = time.time()  # Record start time

            mean = np.mean(self.data, axis=0)
            std_dev = np.std(self.data, axis=0)

            # Avoid division by zero: set std_dev to 1 where it is zero
            std_dev[std_dev == 0] = 1

            # Apply z-score normalization
            self.data = (self.data - mean) / std_dev

            end_time = time.time()

        self.data = self.data.T

        if verbose:
            print(f"> logs: normalization ({normalizer}) of the data - runtime: {(end_time - start_time):.4f} seconds")

    def plot(self, input_data, incomp_data=None, recov_data=None, nbr_series=None, nbr_val=None, series_range=None,
             subplot=False, size=(16, 8), algorithm=None, save_path="./imputegap_assets", cont_rate=None, display=True, verbose=True):
        """
        Plot the time series data, including raw, contaminated, or imputed data.

        Parameters
        ----------
        input_data : numpy.ndarray
            The original time series data without contamination.
        incomp_data : numpy.ndarray, optional
            The contaminated time series data.
        recov_data : numpy.ndarray, optional
            The imputed time series data.
        nbr_series : int, optional
            The maximum number of series to plot.
        nbr_val : int, optional
            The maximum number of values per series to plot.
        series_range : int, optional
            The index of a specific series to plot. If set, only this series will be plotted.
        subplot : bool, optional
            Print one time series by subplot or all in the same plot.
        size : tuple, optional
            Size of the plot in inches. Default is (16, 8).
        algorithm : str, optional
            Name of the algorithm used for imputation.
        save_path : str, optional
            Path to save the plot locally.
        cont_rate : str, optional
            Percentage of contamination in each series to plot.
        display : bool, optional
            Whether to display the plot. Default is True.
        verbose : bool, optional
            Whether to display the plot information. Default is True.

        Returns
        -------
        str or None
            The file path of the saved plot, if applicable.

        Example
        -------
            >>> ts.plot(input_data=ts.data, nbr_series=9, nbr_val=100, save_path="./imputegap_assets") # plain data
            >>> ts.plot(ts.data, ts_m, nbr_series=9, subplot=True, save_path="./imputegap_assets") # contamination
            >>> ts.plot(input_data=ts.data, incomp_data=ts_m, recov_data=imputer.recov_data, nbr_series=9, subplot=True, save_path="./imputegap_assets") # imputation
        """
        select_backend()

        number_of_series = 0
        if algorithm is None:
            algorithm = "imputegap"
            title_imputation = "Imputed Data"
            title_contamination = "Missing Data"
        else:
            title_imputation = algorithm.lower()
            title_contamination = algorithm.lower()

        if nbr_series is None or nbr_series == -1:
            nbr_series = input_data.shape[0]
        if nbr_val is None or nbr_val == -1:
            nbr_val = input_data.shape[1]

        if subplot:
            series_indices = [i for i in range(incomp_data.shape[0]) if np.isnan(incomp_data[i]).any()]
            count_series = [series_range] if series_range is not None else range(min(len(series_indices), nbr_series))
            n_series_to_plot = len(count_series)
        else:
            series_indices = [series_range] if series_range is not None else range(min(input_data.shape[0], nbr_series))
            n_series_to_plot = len(series_indices)

        if n_series_to_plot == 0:
            n_series_to_plot = min(nbr_series, incomp_data.shape[0])

        if subplot:
            n_cols = min(3, n_series_to_plot)
            n_rows = (n_series_to_plot + n_cols - 1) // n_cols

            x_size, y_size = size
            x_size = x_size * n_cols
            y_size = y_size * n_rows

            scale_factor = 0.85
            x_size_screen = (1920 / 100) * scale_factor
            y_size_screen = (1080 / 100) * scale_factor

            if n_rows < 4:
                x_size = x_size_screen
                y_size = y_size_screen

            fig, axes = plt.subplots(n_rows, n_cols, figsize=(x_size, y_size), squeeze=False)
            fig.canvas.manager.set_window_title(algorithm)
            axes = axes.flatten()
        else:
            plt.figure(figsize=size)
            plt.grid(True, linestyle='--', color='#d3d3d3', linewidth=0.6)

        if input_data is not None:
            colors = utils.load_parameters("default", algorithm="colors", verbose=False)

            for idx, i in enumerate(series_indices):

                if subplot:
                    color = colors[0]
                else:
                    color = colors[i % len(colors)]

                timestamps = np.arange(min(input_data.shape[1], nbr_val))

                # Select the current axes if using subplots
                if subplot:
                    ax = axes[idx]
                    ax.grid(True, linestyle='--', color='#d3d3d3', linewidth=0.6)
                else:
                    ax = plt

                if incomp_data is None and recov_data is None:  # plot only raw matrix
                    ax.plot(timestamps, input_data[i, :nbr_val], linewidth=2.5,
                            color=color, linestyle='-', label=f'Series_' + str(i+1))

                if incomp_data is not None and recov_data is None:  # plot infected matrix
                    if np.isnan(incomp_data[i, :]).any():
                        ax.plot(timestamps, input_data[i, :nbr_val], linewidth=2,
                                color="blue", linestyle='--', label=title_contamination)

                    if np.isnan(incomp_data[i, :]).any() or not subplot:
                        ax.plot(np.arange(min(incomp_data.shape[1], nbr_val)), incomp_data[i, :nbr_val],
                                color=color, linewidth=7, linestyle='-', label=f'Series')

                if recov_data is not None:  # plot imputed matrix
                    if np.isnan(incomp_data[i, :]).any():
                        ax.plot(np.arange(min(recov_data.shape[1], nbr_val)), recov_data[i, :nbr_val],
                                linestyle='-', color="r", label=title_imputation)

                        ax.plot(timestamps, input_data[i, :nbr_val], linewidth=1.5,
                                linestyle='--', color=color, label=f'Missing Data')

                    if np.isnan(incomp_data[i, :]).any() or not subplot:
                        ax.plot(np.arange(min(incomp_data.shape[1], nbr_val)), incomp_data[i, :nbr_val],
                                color=color, linewidth=2.5, linestyle='-', label=f'Series')

                # Label and legend for subplot
                if subplot:
                    ax.set_title('Series ' + str(i+1), fontsize=9)
                    #ax.plot([], [], ' ', label='Series ' + str(i + 1))  # invisible line with label
                    ax.set_xlabel('Timestamp', fontsize=7)
                    ax.set_ylabel('Values', fontsize=7)
                    ax.legend(loc='upper left', fontsize=6, frameon=True, fancybox=True, framealpha=0.8)
                    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
                    fig.subplots_adjust(top=0.96, hspace=0.4)
                else:
                    plt.tight_layout(rect=[0.01, 0.03, 0.88, 0.95])

                number_of_series += 1
                if number_of_series == nbr_series:
                    break

        if subplot:
            for idx in range(len(series_indices), len(axes)):
                axes[idx].axis('off')

        if not subplot:
            plt.xlabel('Timestamp')
            plt.ylabel('Values')
            plt.legend(
                loc='upper left',
                fontsize=10,
                frameon=True,
                fancybox=True,
                shadow=True,
                borderpad=1.5,
                bbox_to_anchor=(1.02, 1),  # Adjusted to keep the legend inside the window
            )

        file_path = None

        if save_path:
            os.makedirs(save_path, exist_ok=True)

            now = datetime.datetime.now()
            current_time = now.strftime("%y_%m_%d_%H_%M_%S")

            if cont_rate is None:
                file_path = os.path.join(save_path + "/" + current_time + "_" + algorithm + "_plot.jpg")
            else:
                file_path = os.path.join(save_path + "/" + cont_rate + "_" + algorithm + "_plot.jpg")

            plt.savefig(file_path, bbox_inches='tight')

            if verbose:
                print("\nplots saved in:", file_path)

        if display:
            plt.show()

        self.plots = plt

        return file_path

    class Contamination:
        """
        Inner class to apply contamination patterns to selected series.

        Methods
        -------
        mcar(ts, series_rate=0.2, missing_rate=0.2, block_size=10, offset=0.1, seed=True, explainer=False, verbose=True) :
            Apply Missing Completely at Random (MCAR) contamination to selected series.

        aligned(ts, series_rate=0.2, missing_rate=0.2, offset=0.1) :
            Apply missing percentage contamination to selected series.

        blackout(ts, missing_rate=0.2, offset=0.1) :
            Apply blackout contamination to selected series.

        gaussian(input_data, series_rate=0.2, missing_rate=0.2, std_dev=0.2, offset=0.1, seed=True, verbose=True):
            Apply Gaussian contamination to selected series.

        distribution(input_data, rate_dataset=0.2, rate_series=0.2, probabilities=None, offset=0.1, seed=True, verbose=True):
            Apply any distribution contamination to the time series data based on their probabilities.

        disjoint(input_data, missing_rate=0.1, limit=1, offset=0.1, verbose=True):
            Apply Disjoint contamination to selected series.

        overlap(input_data, missing_rate=0.2, limit=1, shift=0.05, offset=0.1, verbose=True):
            Apply Overlapping contamination to selected series.

        References
        ----------
            https://imputegap.readthedocs.io/en/latest/patterns.html

        """

        def mcar(input_data, rate_dataset=0.2, rate_series=0.2, block_size=10, offset=0.1, seed=True, explainer=False, verbose=True):
            """
            Apply Missing Completely at Random (MCAR) contamination to selected series.

            Parameters
            ----------
            input_data : numpy.ndarray
                The time series dataset to contaminate.
            rate_dataset : float, optional
                Percentage of series to contaminate (default is 0.2).
            rate_series : float, optional
                Percentage of missing values per series (default is 0.2).
            block_size : int, optional
                Size of the block of missing data (default is 10).
            offset : float, optional
                Size of the uncontaminated section at the beginning of the series (default is 0.1).
            seed : bool, optional
                Whether to use a seed for reproducibility (default is True).
            explainer : bool, optional
                Only used within the Explainer Module to contaminate one series at a time (default: False).
            verbose : bool, optional
                Whether to display the contamination information (default is True).

            Returns
            -------
            numpy.ndarray
                The contaminated time series data.

            Example
            -------
                >>> ts_m = ts.Contamination.mcar(ts.data, rate_dataset=0.2, rate_series=0.4, block_size=10):

            """

            if seed:
                seed_value = 42
                np.random.seed(seed_value)
            else:
                seed_value = -1

            ts_contaminated = input_data.copy()
            M, NS = ts_contaminated.shape

            if not explainer:  # use random series
                rate_series = utils.verification_limitation(rate_series)
                rate_dataset = utils.verification_limitation(rate_dataset)
                offset = utils.verification_limitation(offset)

                nbr_series_impacted = int(np.ceil(M * rate_dataset))
                series_selected = [str(idx) for idx in np.random.choice(M, nbr_series_impacted, replace=False)]

            else:  # use fix series
                series_selected = [str(rate_dataset)]

            offset_nbr = int(offset * NS)
            values_nbr = int(NS * rate_series)

            if not explainer and verbose:
                print(f"\n(CONT) missigness pattern: MCAR"
                      f"\n\tselected series: {', '.join(str(int(n)+1) for n in sorted(series_selected, key=int))}"
                      f"\n\tpercentage of contaminated series: {rate_dataset * 100}%"
                      f"\n\trate of missing data per series: {rate_series * 100}%"
                      f"\n\tblock size: {block_size}"
                      f"\n\tsecurity offset: [0-{offset_nbr}]"
                      f"\n\tseed value: {seed_value}")


            if offset_nbr + values_nbr > NS:
                raise ValueError(
                    f"\n\tError: The sum of offset ({offset_nbr}) and missing values ({values_nbr}) exceeds the limit of of the series."
                    f" ({offset_nbr+values_nbr} must be smaller than {NS}).")


            for series in series_selected:
                S = int(series)
                N = len(ts_contaminated[S])  # number of values in the series
                P = int(N * offset)  # values to protect in the beginning of the series
                W = int(N * rate_series)  # number of data to remove
                B = int(W / block_size)  # number of block to remove

                if B <= 0:
                    raise ValueError("The number of block to remove must be greater than 0. "
                                     "The dataset or the number of blocks may not be appropriate."
                                     "One series has", str(N), "population is ", str((N - P)), "the number to remove",
                                     str(W), "and block site", str(block_size), "")

                data_to_remove = np.random.choice(range(P, N), B, replace=False)

                for start_point in data_to_remove:
                    for jump in range(block_size):  # remove the block size for each random position
                        position = start_point + jump

                        if position >= N:  # If block exceeds the series length
                            position = P + (position - N)  # Wrap around to the start after protection

                        while np.isnan(ts_contaminated[S, position]):
                            position = position + 1

                            if position >= N:  # If block exceeds the series length
                                position = P + (position - N)  # Wrap around to the start after protection

                        ts_contaminated[S, position] = np.nan

            return ts_contaminated

        def aligned(input_data, rate_dataset=0.2, rate_series=0.2, offset=0.1, explainer=False, verbose=True):
            """
            Create aligned missing blocks across the selected series.

            Parameters
            ----------
            input_data : numpy.ndarray
                The time series dataset to contaminate.
            rate_dataset : float, optional
                Percentage of series to contaminate (default is 0.2).
            rate_series : float, optional
                Percentage of missing values per series (default is 0.2).
            offset : float, optional
                Size of the uncontaminated section at the beginning of the series (default is 0.1).
            explainer : bool, optional
                Only used within the Explainer Module to contaminate one series at a time (default: False).
            verbose : bool, optional
                Whether to display the contamination information (default is True).

            Returns
            -------
            numpy.ndarray
                The contaminated time series data.

            Example
            -------
                >>> ts_m = ts.Contamination.aligned(ts.data, rate_dataset=0.2, rate_series=0.4, offset=0.1):

            """

            ts_contaminated = input_data.copy()
            M, NS = ts_contaminated.shape
            default_init = 0

            offset_nbr = int(offset * NS)
            values_nbr = int(NS * rate_series)

            if not explainer:  # use random series
                rate_series = utils.verification_limitation(rate_series)
                rate_dataset = utils.verification_limitation(rate_dataset)
                offset = utils.verification_limitation(offset)
                nbr_series_impacted = int(np.ceil(M * rate_dataset))
            else:  # use fix series
                nbr_series_impacted = int(rate_dataset)
                default_init = nbr_series_impacted
                nbr_series_impacted = nbr_series_impacted + 1

            if not explainer and verbose:
                print(f"\n(CONT) missigness pattern: ALIGNED"
                      f"\n\tpercentage of contaminated series: {rate_dataset * 100}%"
                      f"\n\trate of missing data per series: {rate_series * 100}%"
                      f"\n\tsecurity offset: [0-{offset_nbr}]"
                      f"\n\tindex impacted : {offset_nbr} -> {offset_nbr + values_nbr}")

            if offset_nbr + values_nbr > NS:
                raise ValueError(
                    f"\n\tError: The sum of offset ({offset_nbr}) and missing values ({values_nbr}) exceeds the limit of of the series."
                    f" ({offset_nbr+values_nbr} must be smaller than {NS}).")


            for series in range(default_init, nbr_series_impacted):
                S = int(series)
                N = len(ts_contaminated[S])  # number of values in the series
                P = int(N * offset)  # values to protect in the beginning of the series
                W = int(N * rate_series)  # number of data to remove

                for to_remove in range(0, W):
                    index = P + to_remove
                    ts_contaminated[S, index] = np.nan

            return ts_contaminated

        def scattered(input_data, rate_dataset=0.2, rate_series=0.2, offset=0.1, seed=True, explainer=False, verbose=True):
            """
            Apply percentage shift contamination with random starting position to selected series.

            Parameters
            ----------
            input_data : numpy.ndarray
                The time series dataset to contaminate.
            rate_dataset : float, optional
                Percentage of series to contaminate (default is 0.2).
            rate_series : float, optional
                Percentage of missing values per series (default is 0.2).
            offset : float, optional
                Size of the uncontaminated section at the beginning of the series (default is 0.1).
            seed : bool, optional
                Whether to use a seed for reproducibility (default is True).
            explainer : bool, optional
                Only used within the Explainer Module to contaminate one series at a time (default: False).
            verbose : bool, optional
                Whether to display the contamination information (default is True).

            Returns
            -------
            numpy.ndarray
                The contaminated time series data.

            Example
            -------
                >>> ts_m = ts.Contamination.scattered(ts.data, rate_dataset=0.2, rate_series=0.4, offset=0.1)

            """

            if seed:
                seed_value = 42
                np.random.seed(seed_value)

            ts_contaminated = input_data.copy()
            M, NS = ts_contaminated.shape
            default_init = 0

            offset_nbr = int(offset * NS)
            values_nbr = int(NS * rate_series)

            if not explainer:  # use random series
                rate_series = utils.verification_limitation(rate_series)
                rate_dataset = utils.verification_limitation(rate_dataset)
                offset = utils.verification_limitation(offset)
                nbr_series_impacted = int(np.ceil(M * rate_dataset))
            else:  # use fix series
                nbr_series_impacted = int(rate_dataset)
                default_init = nbr_series_impacted
                nbr_series_impacted = nbr_series_impacted + 1

            if not explainer and verbose:
                print(f"\n(CONT) missigness pattern: SCATTER"
                      f"\n\tpercentage of contaminated series: {rate_dataset * 100}%"
                      f"\n\trate of missing data per series: {rate_series * 100}%"
                      f"\n\tsecurity offset: [0-{offset_nbr}]"
                      f"\n\tindex impacted : {offset_nbr} -> {offset_nbr + values_nbr}")


            if offset_nbr + values_nbr > NS:
                raise ValueError(
                    f"\n\tError: The sum of offset ({offset_nbr}) and missing values ({values_nbr}) exceeds the limit of of the series."
                    f" ({offset_nbr+values_nbr} must be smaller than {NS}).")


            for series in range(default_init, nbr_series_impacted):
                S = int(series)
                N = len(ts_contaminated[S])  # number of values in the series
                P = int(N * offset)  # values to protect in the beginning of the series
                W = int(N * rate_series)  # number of data to remove
                L = (N - W - P) +1

                start_index = np.random.randint(0, L)  # Random start position

                for to_remove in range(0, W):
                    index = P + start_index + to_remove
                    ts_contaminated[S, index] = np.nan

            return ts_contaminated

        def blackout(input_data, series_rate=0.2, offset=0.1, verbose=True):
            """
            Apply blackout contamination to selected series

            Parameters
            ----------
            input_data : numpy.ndarray
                The time series dataset to contaminate.
            series_rate : float, optional
                Percentage of missing values per series (default is 0.2).
            offset : float, optional
                Size of the uncontaminated section at the beginning of the series (default is 0.1).
            verbose : bool, optional
                Whether to display the contamination information (default is True).

            Returns
            -------
            numpy.ndarray
                The contaminated time series data.

            Example
            -------
                >>> ts_m = ts.Contamination.blackout(ts.data, series_rate=0.2)

            """
            return TimeSeries.Contamination.aligned(input_data, rate_dataset=1, rate_series=series_rate, offset=offset, verbose=verbose)

        def gaussian(input_data, rate_dataset=0.2, rate_series=0.2, std_dev=0.2, offset=0.1, seed=True, explainer=False, verbose=True):
            """
            Apply contamination with a Gaussian distribution to selected series

            Parameters
            ----------
            input_data : numpy.ndarray
                The time series dataset to contaminate.
            rate_dataset : float, optional
                Percentage of series to contaminate (default is 0.2).
            rate_series : float, optional
                Percentage of missing values per series (default is 0.2).
            std_dev : float, optional
                Standard deviation of the Gaussian distribution for missing values (default is 0.4).
            offset : float, optional
                Size of the uncontaminated section at the beginning of the series (default is 0.1).
            seed : bool, optional
                Whether to use a seed for reproducibility (default is True).
            explainer : bool, optional
                Only used within the Explainer Module to contaminate one series at a time (default: False).
            verbose : bool, optional
                Whether to display the contamination information (default is True).

            Returns
            -------
            numpy.ndarray
                The contaminated time series data.

            Example
            -------
                >>> ts_m = ts.Contamination.gaussian(ts.data, rate_series=0.2, std_dev=0.4, offset=0.1):

            """
            from scipy.stats import norm

            ts_contaminated = input_data.copy()
            M, NS = ts_contaminated.shape
            default_init = 0

            if seed:
                seed_value = 42
                np.random.seed(seed_value)

            offset_nbr = int(offset * NS)
            values_nbr = int(NS * rate_series)

            if not explainer:  # use random series
                # Validation and limitation of input parameters
                rate_series = utils.verification_limitation(rate_series)
                rate_dataset = utils.verification_limitation(rate_dataset)
                offset = utils.verification_limitation(offset)
                nbr_series_impacted = int(np.ceil(M * rate_dataset))
            else:  # use fix series
                nbr_series_impacted = int(rate_dataset)
                default_init = nbr_series_impacted
                nbr_series_impacted = nbr_series_impacted + 1

            if not explainer and verbose:
                print(f"\n(CONT) missigness pattern: GAUSSIAN"
                      f"\n\tpercentage of contaminated series: {rate_dataset * 100}%"
                      f"\n\trate of missing data per series: {rate_series * 100}%"
                      f"\n\tsecurity offset: [0-{offset_nbr}]"
                      f"\n\tseed value: {seed_value}"
                      f"\n\tstandard deviation : {std_dev}")

            if offset_nbr + values_nbr > NS:
                raise ValueError(
                    f"\n\tError: The sum of offset ({offset_nbr}) and missing values ({values_nbr}) exceeds the limit of of the series.")


            for series in range(default_init, nbr_series_impacted):
                S = int(series)
                N = len(ts_contaminated[S])  # number of values in the series
                P = int(N * offset)  # values to protect in the beginning of the series
                W = int(N * rate_series)  # number of data points to remove
                R = np.arange(P, N)

                # probability density function
                mean = np.mean(ts_contaminated[S])
                mean = max(min(mean, 1), -1)

                probabilities = norm.pdf(R, loc=P + mean * (N - P), scale=std_dev * (N - P))

                # normalizes the probabilities so that their sum equals 1
                probabilities /= probabilities.sum()

                # select the values based on the probability
                missing_indices = np.random.choice(R, size=W, replace=False, p=probabilities)

                # apply missing values
                ts_contaminated[S, missing_indices] = np.nan

            return ts_contaminated

        def distribution(input_data, rate_dataset=0.2, rate_series=0.2, probabilities_list=None, offset=0.1, seed=True, explainer=False, verbose=True):
            """
            Apply any distribution contamination to the time series data based on their probabilities.

            Parameters
            ----------
            input_data : numpy.ndarray
                The time series dataset to contaminate.
            rate_dataset : float, optional
                Percentage of series to contaminate (default is 0.2).
            rate_series : float, optional
                Percentage of missing values per series (default is 0.2).
            probabilities_list : 2-D array-like, optional
                The probabilities of being contaminated associated with each values of a series.
                Most match the shape of input data without the offset : (e.g. [[0.1, 0, 0.3, 0], [0.2, 0.1, 0.2, 0.9]])
            offset : float, optional
                Size of the uncontaminated section at the beginning of the series (default is 0.1).
            seed : bool, optional
                Whether to use a seed for reproducibility (default is True).
            explainer : bool, optional
                Only used within the Explainer Module to contaminate one series at a time (default: False).
            verbose : bool, optional
                Whether to display the contamination information (default is True).

            Returns
            -------
            numpy.ndarray
                The contaminated time series data.

            Example
            -------
                >>> ts_m = ts.Contamination.distribution(ts.data, rate_dataset=0.2, rate_series=0.2, probabilities_list=probabilities_list, offset=0.1)

            """

            ts_contaminated = input_data.copy()
            M, NS = ts_contaminated.shape
            default_init = 0

            if seed:
                seed_value = 42
                np.random.seed(seed_value)

            offset_nbr = int(offset * NS)
            values_nbr = int(NS * rate_series)

            if not explainer:  # use random series
                # Validation and limitation of input parameters
                rate_series = utils.verification_limitation(rate_series)
                rate_dataset = utils.verification_limitation(rate_dataset)
                offset = utils.verification_limitation(offset)
                nbr_series_impacted = int(np.ceil(M * rate_dataset))
            else:  # use fix series
                nbr_series_impacted = int(rate_dataset)
                default_init = nbr_series_impacted
                nbr_series_impacted = nbr_series_impacted + 1

            if not explainer and verbose:
                print(f"\n(CONT) missigness pattern: DISTRIBUTION"
                      f"\n\tpercentage of contaminated series: {rate_dataset * 100}%"
                      f"\n\trate of missing data per series: {rate_series * 100}%"
                      f"\n\tsecurity offset: [0-{offset_nbr}]"
                      f"\n\tseed value: {seed_value}"
                      f"\n\tprobabilities list : {np.array(probabilities_list).shape}")

            if offset_nbr + values_nbr > NS:
                raise ValueError(
                    f"\n\tError: The sum of offset ({offset_nbr}) and missing values ({values_nbr}) exceeds the limit of of the series.")

            if np.array(probabilities_list).shape != (M, NS - offset_nbr):
                raise ValueError(
                    f"\n\tError: The probability list does not match the matrix in input {np.array(probabilities_list).shape} != ({M},{NS - offset_nbr}).")

            for series in range(default_init, nbr_series_impacted):
                S = int(series)
                N = len(ts_contaminated[S])  # number of values in the series
                P = int(N * offset)  # values to protect in the beginning of the series
                W = int(N * rate_series)  # number of data points to remove
                R = np.arange(P, N)
                D = probabilities_list[S]

                missing_indices = np.random.choice(R, size=W, replace=False, p=D)

                # apply missing values
                ts_contaminated[S, missing_indices] = np.nan

            return ts_contaminated


        def disjoint(input_data, rate_series=0.1, limit=1, offset=0.1, verbose=True):
            """
            Apply disjoint contamination to selected series

            Parameters
            ----------
            input_data : numpy.ndarray
                The time series dataset to contaminate.
            rate_series : float, optional
                Percentage of missing values per series (default is 0.1).
            limit : float, optional
                Percentage expressing the limit index of the end of the contamination (default is 1: all length).
            offset : float, optional
                Size of the uncontaminated section at the beginning of the series (default is 0.1).
            verbose : bool, optional
                Whether to display the contamination information (default is True).

            Returns
            -------
            numpy.ndarray
                The contaminated time series data.

            Example
            -------
                >>> ts_m = ts.Contamination.disjoint(ts.data, rate_series=0.1, limit=1, offset=0.1)

            """
            ts_contaminated = input_data.copy()
            M, NS = ts_contaminated.shape

            rate_series = utils.verification_limitation(rate_series)
            offset = utils.verification_limitation(offset)

            offset_nbr = int(offset * NS)
            values_nbr = int(NS * rate_series)

            if verbose:
                print(f"\n(CONT) missigness pattern: DISJOINT"
                      f"\n\tpercentage of contaminated series: {rate_series * 100}%"
                      f"\n\trate of missing data per series: {rate_series * 100}%"
                      f"\n\tsecurity offset: [0-{offset_nbr}]"
                      f"\n\tlimit: {limit}")

            if offset_nbr + values_nbr > NS:
                raise ValueError(
                    f"\n\tError: The sum of offset ({offset_nbr}) and missing values ({values_nbr}) exceeds the limit of of the series.")

            S = 0
            X = 0
            final_limit = int(NS*limit)-1

            while S < M:
                N = len(ts_contaminated[S])  # number of values in the series
                P = int(N * offset)  # values to protect in the beginning of the series
                W = int(N * rate_series)  # number of data to remove
                L = X + W  # new limit

                for to_remove in range(X, L):
                    index = P + to_remove
                    ts_contaminated[S, index] = np.nan

                    if index >= final_limit:  # reach the limitation
                        return ts_contaminated

                X = L
                S = S + 1

            return ts_contaminated

        def overlap(input_data, rate_series=0.2, limit=1, shift=0.05, offset=0.1, verbose=True):
            """
            Apply overlap contamination to selected series

            Parameters
            ----------
            input_data : numpy.ndarray
                The time series dataset to contaminate.
            rate_series : float, optional
                Percentage of missing values per series (default is 0.2).
            limit : float, optional
                Percentage expressing the limit index of the end of the contamination (default is 1: all length).
            shift : float, optional
                Percentage of shift inside each the last disjoint contamination.
            offset : float, optional
                Size of the uncontaminated section at the beginning of the series (default is 0.1).
            verbose : bool, optional
                Whether to display the contamination information (default is True).

            Returns
            -------
            numpy.ndarray
                The contaminated time series data.

            Example
            -------
                >>> ts_m = ts.Contamination.overlap(ts.data, rate_series=0.1, limit=1, shift=0.05, offset=0.1)

            """
            ts_contaminated = input_data.copy()
            M, NS = ts_contaminated.shape

            rate_series = utils.verification_limitation(rate_series)
            offset = utils.verification_limitation(offset)

            offset_nbr = int(offset * NS)
            values_nbr = int(NS * rate_series)

            if verbose:
                print(f"\n(CONT) missigness pattern: OVERLAP"
                      f"\n\tpercentage of contaminated series: {rate_series * 100}%"
                      f"\n\trate of missing data per series: {rate_series * 100}%"
                      f"\n\tsecurity offset: [0-{offset_nbr}]"
                      f"\n\tshift: {shift * 100} %"
                      f"\n\tlimit: {limit}")


            if offset_nbr + values_nbr > NS:
                raise ValueError(
                    f"\n\tError: The sum of offset ({offset_nbr}) and missing values ({values_nbr}) exceeds the limit of of the series.")

            if int(NS*shift) > int(NS*offset):
                raise ValueError(f"Shift too big for this dataset and offset: shift ({int(NS*shift)}), offset ({int(NS*offset)}).")

            S, X = 0, 0
            final_limit = int(NS * limit) - 1

            while S < M:
                N = len(ts_contaminated[S])  # number of values in the series
                P = int(N * offset)  # values to protect in the beginning of the series
                W = int(N * rate_series)  # number of data to remove

                if X != 0:
                    X = X - int(N * shift)

                L = X + W  # new limit

                for to_remove in range(X, L):
                    index = P + to_remove
                    ts_contaminated[S, index] = np.nan

                    if index >= final_limit:  # reach the limitation
                        return ts_contaminated

                X = L
                S = S + 1

            return ts_contaminated