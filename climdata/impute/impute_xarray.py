"""Gap filling for gridded datasets, backed by the vendored imputegap library.

Each ``(variable, lat, lon)`` cell is treated as an independent time series: the
dataset is flattened to a 2-D ``(n_series, seq_len)`` matrix, handed to one of
imputegap's algorithms, and folded back into its original shape. Nothing exploits
spatial correlation between neighbouring cells — a cell with no valid values of
its own cannot be recovered from its neighbours.

The deep-learning methods (BRITS, GRIN, BayOTIDE, DeepMVI) need PyTorch, which
climdata does not install by default; the matrix-completion and tree methods
(SoftImpute, CDRec, XGBOOST) do not.

Example:
    >>> imputer = Imputer(ds, method="SoftImpute")     # doctest: +SKIP
    >>> filled = imputer.impute()                      # doctest: +SKIP
    >>> imputer.summary()                              # doctest: +SKIP
"""

import warnings
import numpy as np
import xarray as xr
import logging

try:
    from .._vendor.imputegap.recovery.imputation import Imputation
    from .._vendor.imputegap.recovery.manager import TimeSeries
    _IMPUTEGAP_AVAILABLE = True
except ImportError:
    _IMPUTEGAP_AVAILABLE = False
    Imputation = None
    TimeSeries = None

try:
    import pycatch22  # noqa: F401
    _PYCATCH22_AVAILABLE = True
except ImportError:
    _PYCATCH22_AVAILABLE = False
    warnings.warn(
        "[climdata.impute] pycatch22 is not installed. "
        "Some imputation feature-extraction methods will be unavailable. "
        "Install it with: conda install -c conda-forge pycatch22",
        ImportWarning,
        stacklevel=2,
    )

try:
    from sklearn.metrics import mean_squared_error, mean_absolute_error
    _SKLEARN_AVAILABLE = True
except ImportError:
    _SKLEARN_AVAILABLE = False
    warnings.warn(
        "[climdata.impute] scikit-learn is not installed. "
        "Install it with: pip install scikit-learn",
        ImportWarning,
        stacklevel=2,
    )

logger = logging.getLogger(__name__)

# ANSI color codes (works on Linux terminals)
_COLOR_YELLOW = "\033[93m"
_COLOR_RESET = "\033[0m"

class Imputer:
    """Fill missing values in a ``(time, lat, lon)`` dataset, cell by cell.

    Reshapes the dataset into one time series per variable and grid cell, runs an
    imputegap algorithm over the resulting matrix, and reassembles the output.
    Optional z-score normalisation is applied before imputation and inverted
    afterwards, which matters for the neural methods and is harmless otherwise.

    Attributes:
        ds (xr.Dataset): The input dataset.
        method (str): Algorithm name — ``"BRITS"``, ``"SoftImpute"``,
            ``"CDRec"`` or ``"XGBOOST"``.
        normalize (bool): Whether series are z-scored before imputation.
        variables (list[str]): Data variables that will be imputed.
        recovered_ds (xr.Dataset | None): Result of the last :meth:`impute`.

    Example:
        >>> imputer = Imputer(ds, method="CDRec")          # doctest: +SKIP
        >>> imputer.missing_fraction()["global"]           # doctest: +SKIP
        0.043
        >>> filled = imputer.impute()                      # doctest: +SKIP
    """

    # Methods that require PyTorch (BRITS, GRIN, BayOTIDE, DeepMVI)
    _TORCH_METHODS = {"BRITS", "GRIN", "BayOTIDE", "DeepMVI"}

    def __init__(
        self,
        ds: xr.Dataset,
        time_dim: str = "time",
        lat_dim: str = "lat",
        lon_dim: str = "lon",
        method: str = "BRITS",
        normalize: bool = True,
    ):
        """Bind a dataset and check that the chosen method can actually run.

        The PyTorch check happens here rather than at :meth:`impute`, so a
        missing backend is reported before any reshaping work is done.

        Args:
            ds (xr.Dataset): Dataset with time, latitude and longitude dimensions.
            time_dim (str): Name of the time dimension. Defaults to ``"time"``.
            lat_dim (str): Name of the latitude dimension. Defaults to ``"lat"``.
            lon_dim (str): Name of the longitude dimension. Defaults to ``"lon"``.
            method (str): Algorithm to use. Defaults to ``"BRITS"``.
            normalize (bool): Z-score each series before imputing, and invert
                afterwards. Defaults to ``True``.

        Raises:
            ImportError: If the vendored imputegap is unavailable, or if
                ``method`` needs PyTorch and PyTorch is not installed.
        """
        if not _IMPUTEGAP_AVAILABLE:
            raise ImportError(
                "[climdata.impute] imputegap is not available. "
                "It is bundled with climdata — please reinstall the package."
            )
        if method in self._TORCH_METHODS:
            try:
                import torch  # noqa: F401
            except ImportError:
                raise ImportError(
                    f"[climdata.impute] Method '{method}' requires PyTorch, which is not installed. "
                    f"Install it with:\n"
                    f"  pip install torch torchvision torchaudio "
                    f"--index-url https://download.pytorch.org/whl/cpu\n"
                    f"See docs/installation.md for GPU and full dependency instructions."
                )
        self.ds = ds
        self.time_dim = time_dim
        self.lat_dim = lat_dim
        self.lon_dim = lon_dim
        self.method = method
        self.normalize = normalize

        self.variables = list(ds.data_vars)
        self.coords = ds.coords
        self.attrs = ds.attrs

        self.ts = None
        self.mask = None
        self.shape_info = None
        self.recovered_ds = None
    def missing_fraction(self):
        """Report how much data is missing, per variable and overall.

        Returns:
            dict[str, float]: One entry per data variable plus a ``"global"``
            key. The global figure is the unweighted mean of the per-variable
            fractions, so variables of different sizes count equally.
        """
        frac = {
            v: float(self.ds[v].isnull().mean())
            for v in self.variables
        }
        frac["global"] = float(
            np.mean(list(frac.values()))
        )
        return frac
    def _to_timeseries(self):
        """Flatten the dataset into the ``(n_series, seq_len)`` matrix imputegap wants.

        Series are ordered ``(variable, lat, lon)``, and the shape and coordinate
        values are recorded on ``self.shape_info`` so :meth:`_from_timeseries`
        can invert the operation. The NaN mask is captured before any
        normalisation.

        Returns:
            None: Sets ``self.ts``, ``self.mask`` and ``self.shape_info``.
        """
        da = self.ds[self.variables].to_array("variable")
        # dims: (variable, time, lat, lon)

        da = da.transpose(
            "variable",
            self.time_dim,
            self.lat_dim,
            self.lon_dim,
        )

        data = da.values
        n_var, t, n_lat, n_lon = data.shape

        self.shape_info = {
            "variables": np.array(self.variables, dtype=object),
            self.time_dim: self.ds[self.time_dim].values,
            self.lat_dim: self.ds[self.lat_dim].values,
            self.lon_dim: self.ds[self.lon_dim].values,
            "n_var": n_var,
            "n_lat": n_lat,
            "n_lon": n_lon,
        }

        data_2d = data.reshape(
            n_var * n_lat * n_lon,
            t
        )

        self.mask = np.isnan(data_2d)

        ts = TimeSeries()
        ts.data = data_2d
        ts.n_series, ts.seq_len = ts.data.shape

        if self.normalize:
            ts.normalize(normalizer="z_score")

        self.ts = ts
    def _from_timeseries(self, data_2d):
        """Rebuild a dataset from the imputed matrix.

        Coordinates are restored from ``self.shape_info``. Any coordinate that
        does not match the expected length is replaced by a plain integer range
        rather than raising, so a shape surprise costs coordinate labels rather
        than the whole result.

        Args:
            data_2d (numpy.ndarray): Imputed matrix of shape
                ``(n_series, seq_len)``.

        Returns:
            xr.Dataset: One variable per input variable, dimensioned
            ``(time, lat, lon)``, carrying the original global attributes.
        """
        n_var = self.shape_info["n_var"]
        n_lat = self.shape_info["n_lat"]
        n_lon = self.shape_info["n_lon"]

        # rebuild 4D array to (variable, time, lat, lon)
        data_4d = data_2d.reshape(
            n_var,
            n_lat,
            n_lon,
            -1
        ).transpose(0, 3, 1, 2)
        # (variable, time, lat, lon)

        # --- sanitize coords (ensure 1-D arrays and correct lengths) ---
        # helper to coerce to 1D and validate length
        def _safe_coord(key, expected_len, fallback_range=True):
            val = self.shape_info.get(key, None)
            val = np.asarray(val) if val is not None else None
            if val is None:
                return np.arange(expected_len) if fallback_range else None
            # flatten to 1D if possible
            if val.ndim != 1 or len(val) != expected_len:
                # try ravel() then trim / pad if needed
                val = val.ravel()
                if len(val) >= expected_len:
                    return val[:expected_len]
                else:
                    # fallback to integer index
                    return np.arange(expected_len)
            return val

        variable_coord = _safe_coord("variables", n_var)
        time_coord = _safe_coord(self.time_dim, data_4d.shape[1])
        lat_coord = _safe_coord(self.lat_dim, n_lat)
        lon_coord = _safe_coord(self.lon_dim, n_lon)

        da = xr.DataArray(
            data_4d,
            dims=("variable", self.time_dim, self.lat_dim, self.lon_dim),
            coords=(
                ("variable", np.asarray(variable_coord).astype(str)),
                (self.time_dim, time_coord),
                (self.lat_dim, lat_coord),
                (self.lon_dim, lon_coord),
            ),
            attrs=self.attrs,
        )

        return da.to_dataset(dim="variable")
    def impute(self, epochs: int = 300):
        """Fill the missing values and return the completed dataset.

        Returns a copy unchanged when there is nothing missing, so calling this
        unconditionally in a workflow is safe and cheap.

        Note that every value is reconstructed, not just the gaps: the returned
        dataset carries the algorithm's estimate everywhere. Use :meth:`metrics`
        to see how closely it reproduces the values that were present.

        Args:
            epochs (int): Training epochs, used only by ``BRITS``. Defaults to ``300``.

        Returns:
            xr.Dataset: The completed dataset, also stored on
            :attr:`recovered_ds`.

        Raises:
            ValueError: If :attr:`method` is not one of the supported algorithms.
        """
        if self.missing_fraction()["global"] == 0.0:
            logger.info(f"{_COLOR_YELLOW}No missing data found. Imputation not required.{_COLOR_RESET}")
            self.recovered_ds = self.ds.copy(deep=True)
            return self.recovered_ds

        self._to_timeseries()

        data = self.ts.data
        # print(data.shape)
        if self.method == "BRITS":
            imputer = Imputation.DeepLearning.BRITS(data)
            imputer.epochs = epochs
        elif self.method == "SoftImpute":
            imputer = Imputation.MatrixCompletion.SoftImpute(data)
        elif self.method == "CDRec":
            imputer = Imputation.MatrixCompletion.CDRec(data)
        elif self.method == "XGBOOST":
            imputer = Imputation.MachineLearning.XGBOOST(data)
        else:
            raise ValueError(f"Unknown method: {self.method}")

        imputer.impute()
        rec = imputer.recov_data

        if self.normalize:
            mean = self.ts.data.mean(axis=1, keepdims=True)
            std = self.ts.data.std(axis=1, keepdims=True)
            rec = rec * std + mean

        self.recovered_ds = self._from_timeseries(rec)
        return self.recovered_ds
    def metrics(self):
        """Score the reconstruction against the values that were actually present.

        Compares imputed against original on the cells that were *not* missing —
        the only cells where ground truth exists. This measures how faithfully
        the algorithm reproduces known data, which is an upper bound on, not a
        measure of, its accuracy in the gaps. For a real estimate, mask out known
        values yourself and score the reconstruction there.

        Variables with no missing values are omitted entirely.

        Returns:
            dict[str, dict]: Per variable, ``rmse``, ``mae`` and
            ``missing_fraction``.

        Raises:
            RuntimeError: If :meth:`impute` has not been called.
        """
        if self.recovered_ds is None:
            raise RuntimeError("Call impute() first")

        scores = {}

        for v in self.variables:
            orig = self.ds[v].values
            rec = self.recovered_ds[v].values
            mask = np.isnan(orig)

            if mask.sum() == 0:
                continue

            scores[v] = {
                "rmse": mean_squared_error(
                    orig[~mask],
                    rec[~mask],
                    squared=False,
                ),
                "mae": mean_absolute_error(
                    orig[~mask],
                    rec[~mask],
                ),
                "missing_fraction": float(mask.mean()),
            }

        return scores
    def summary(self):
        """Describe the imputation setup and the gaps it faces.

        Returns:
            dict: ``method``, ``normalize``, ``variables``, ``missing_fraction``
            and ``dims``. Cheap to call before :meth:`impute` to check what will
            be run.
        """
        return {
            "method": self.method,
            "normalize": self.normalize,
            "variables": self.variables,
            "missing_fraction": self.missing_fraction(),
            "dims": dict(self.ds.dims),
        }
