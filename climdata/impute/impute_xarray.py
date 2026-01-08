import numpy as np
import xarray as xr
from .._vendor.imputegap.recovery.imputation import Imputation
from .._vendor.imputegap.recovery.manager import TimeSeries
from sklearn.metrics import mean_squared_error, mean_absolute_error
import logging

logger = logging.getLogger(__name__)

# ANSI color codes (works on Linux terminals)
_COLOR_YELLOW = "\033[93m"
_COLOR_RESET = "\033[0m"

class Imputer:
    """
    Impute missing values in an xarray.Dataset with dims (time, lat, lon)
    using imputegap (BRITS / matrix completion).

    Each (variable, lat, lon) grid cell is treated as an independent time series.
    """

    def __init__(
        self,
        ds: xr.Dataset,
        time_dim: str = "time",
        lat_dim: str = "lat",
        lon_dim: str = "lon",
        method: str = "BRITS",
        normalize: bool = True,
    ):
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
        """
        Missing fraction per variable and globally.
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
        """
        Convert Dataset → 2D array (n_series, seq_len)
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
        """
        Convert (n_series, seq_len) → xarray.Dataset
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
        """
        Run imputation unless missing fraction is zero.
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
        """
        RMSE and MAE on originally missing values only.
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
        return {
            "method": self.method,
            "normalize": self.normalize,
            "variables": self.variables,
            "missing_fraction": self.missing_fraction(),
            "dims": dict(self.ds.dims),
        }
