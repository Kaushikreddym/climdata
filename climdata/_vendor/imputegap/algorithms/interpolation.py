import time
import numpy as np
from scipy.interpolate import interp1d

def interpolation(incomp_data, method="linear", poly_order=2, logs=True, verbose=True):
    """
    Perform imputation using the interpolation algorithm, methods to estimate missing values by looking at
    the known values in a dataset.

    Parameters
    ----------
    incomp_data : numpy.ndarray
        The input matrix with contamination (missing values represented as NaNs).
    method : str, optional
        Interpolation method ("linear", "polynomial", "spline", "nearest").
    poly_order : int, optional
        Polynomial degree for "polynomial" and "spline" methods.
    logs : bool, optional
        Whether to log execution time (default: True).
    verbose : bool, optional
        Whether to display the contamination information (default is True).

    Returns
    -------
    numpy.ndarray
        The imputed matrix with missing values recovered.

    Example
    -------
        >>> recov_data = interpolation(incomp_data, method="linear", poly_order=2)
        >>> print(recov_data)

    """

    if verbose:
        print(f"(IMPUTATION) interpolation\n\tMatrix: {incomp_data.shape[0]}, {incomp_data.shape[1]}\n\tmethod: {method}\n\tpolynomial order: {poly_order}\n")

    start_time = time.time()  # Record start time

    recov_data = np.copy(incomp_data)  # Copy data to avoid modifying original
    num_rows, num_cols = recov_data.shape

    for i in range(num_rows):  # Process row-wise (like pyts)
        row_data = recov_data[i, :]
        mask = ~np.isnan(row_data)  # Identify known values

        if np.all(np.isnan(row_data)):  # Skip row if all values are NaN
            continue

        x_known = np.where(mask)[0]  # Indices of known values
        y_known = row_data[mask]  # Known values
        x_missing = np.where(~mask)[0]  # Indices of missing values

        if method in ["nearest", "linear"]:
            interp_func = interp1d(x_known, y_known, kind=method, fill_value="extrapolate", bounds_error=False)
        elif method == "polynomial":
            if len(x_known) > poly_order:
                interp_func = np.poly1d(np.polyfit(x_known, y_known, poly_order))  # Polynomial fit
            else:
                print(f"Not enough polynomial degree for method {method}, fall into linear")
                interp_func = interp1d(x_known, y_known, kind="linear", fill_value="extrapolate", bounds_error=False)
        elif method == "spline":
            if len(x_known) > poly_order:
                interp_func = interp1d(x_known, y_known, kind="cubic", fill_value="extrapolate", bounds_error=False)
            else:
                print(f"Not enough polynomial degree for method {method}, fall into linear")
                interp_func = interp1d(x_known, y_known, kind="linear", fill_value="extrapolate", bounds_error=False)
        else:
            raise ValueError("Invalid interpolation method. Choose from 'linear', 'polynomial', 'spline', or 'nearest'")

        recov_data[i, x_missing] = interp_func(x_missing)  # Fill missing values

    end_time = time.time()
    if logs and verbose:
        print(f"\n> logs: imputation with interpolation - Execution Time: {(end_time - start_time):.4f} seconds\n")

    return recov_data
