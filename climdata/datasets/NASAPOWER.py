"""NASA POWER point access.

The POWER (Prediction Of Worldwide Energy Resources) daily API serves
agro-climatology time series for a single coordinate. It is the only provider in
climdata that is inherently point-based: there are no files to download and no
grid to subset, so a request returns exactly one location.

Example:
    >>> power = POWER(cfg)          # doctest: +SKIP
    >>> power.fetch()               # doctest: +SKIP
    >>> power.load()                # doctest: +SKIP
    >>> power.ds                    # doctest: +SKIP
"""

import requests
import pandas as pd
import xarray as xr


class POWER:
    """Fetch and assemble a NASA POWER daily time series for one point.

    Follows the same ``fetch`` → ``load`` → ``extract`` lifecycle as the gridded
    providers, so :class:`~climdata.utils.wrapper_workflow.ClimateExtractor` can
    drive it identically. Because the API is point-based, ``extract`` only
    narrows the time axis; there is no spatial subsetting to do.

    Variable names are translated through ``cfg.dsinfo.POWER.variables``: CF names
    (``tasmax``, ``pr``) on the climdata side, POWER ``api_id`` codes
    (``T2M_MAX``, ``PRECTOTCORR``) on the wire.

    Attributes:
        cfg (DictConfig): Hydra configuration supplying ``lat``, ``lon``,
            ``variables``, ``time_range`` and ``dsinfo.POWER``.
        raw (dict | None): The decoded JSON response, set by :meth:`fetch`.
        ds (xr.Dataset | None): The assembled dataset, set by :meth:`load`.
    """

    def __init__(self, cfg):
        """Bind a configuration.

        No network call is made here.

        Args:
            cfg (DictConfig): Configuration with ``lat``, ``lon``, ``variables``,
                ``time_range.start_date`` / ``.end_date``, and a
                ``dsinfo.POWER.variables`` mapping.
        """
        self.cfg = cfg
        self.raw = None
        self.ds = None

    # ---------------------------------------------------------
    # Step 1: Fetch data from NASA POWER
    # ---------------------------------------------------------
    def fetch(self):
        """Request the daily point series from the POWER API.

        Builds a single query for all configured variables over the full
        configured date range and stores the decoded JSON on :attr:`raw`.

        Returns:
            None: The response is stored on :attr:`raw`.

        Raises:
            requests.HTTPError: If the API returns a non-2xx status — most often
                an unsupported variable code or an out-of-range date.
            KeyError: If a configured variable has no ``api_id`` in
                ``cfg.dsinfo.POWER.variables``.
        """
        lat = self.cfg.lat
        lon = self.cfg.lon

        params = self.cfg.variables
        api_vars =  ",".join([
            self.cfg.dsinfo.POWER.variables[v].api_id
            for v in params
        ])
        start = self.cfg.time_range.start_date.replace("-", "")
        end   = self.cfg.time_range.end_date.replace("-", "")

        url = (
            "https://power.larc.nasa.gov/api/temporal/daily/point"
            f"?parameters={api_vars}"
            f"&community=AG"
            f"&latitude={lat}"
            f"&longitude={lon}"
            f"&start={start}"
            f"&end={end}"
            f"&format=JSON"
        )

        r = requests.get(url)
        r.raise_for_status()
        self.raw = r.json()

    # ---------------------------------------------------------
    # Step 2: Load JSON into xarray Dataset
    # ---------------------------------------------------------
    def load(self):
        """Convert the fetched JSON into an :class:`xarray.Dataset`.

        Renames POWER ``api_id`` columns back to CF names, parses the
        ``YYYYMMDD`` keys into a ``time`` coordinate, attaches the request's
        latitude and longitude as scalar coordinates, and copies ``long_name``
        and ``units`` from the configuration onto each variable.

        POWER encodes gaps as the sentinel ``-999``; those values are passed
        through unchanged, so screen for them before computing statistics.

        Returns:
            None: The dataset is stored on :attr:`ds`.

        Raises:
            TypeError: If :meth:`fetch` has not been called (``raw`` is ``None``).
            KeyError: If the response carries no ``properties.parameter`` block,
                which is how the API reports a rejected request.
        """
        data = self.raw["properties"]["parameter"]

        df = pd.DataFrame(data)
        df.index = pd.to_datetime(df.index, format="%Y%m%d")
        var_map = {
            v.api_id: cmip_id
            for cmip_id, v in self.cfg.dsinfo.POWER.variables.items()
        }
        df = df.rename(columns=var_map)
        self.ds = xr.Dataset.from_dataframe(df)

        # Add coords
        self.ds = self.ds.assign_coords(
            latitude=self.cfg.lat,
            longitude=self.cfg.lon
        )
        self.ds = self.ds.rename({"index":"time"})
        self.ds.latitude.attrs["units"] = "degrees_north"
        self.ds.longitude.attrs["units"] = "degrees_east"

        for cmip_id in self.ds.data_vars:
            if cmip_id in self.cfg.dsinfo.POWER.variables:
                vinfo = self.cfg.dsinfo.POWER.variables[cmip_id]

                self.ds[cmip_id].attrs.update({
                    "long_name": vinfo.long_name,
                    "units": vinfo.units,
                    "source": "NASA POWER",
                    "api_id": vinfo.api_id,
                })

    # ---------------------------------------------------------
    # Step 3: Extract (temporal subsetting etc.)
    # ---------------------------------------------------------
    def extract(self, start=None, end=None):
        """Narrow the loaded dataset to a time window.

        There is no spatial argument: the API already returned exactly one point.
        Calling with neither bound is a no-op, which is what lets the generic
        workflow invoke ``extract`` unconditionally.

        Args:
            start (str | datetime, optional): Inclusive start of the window.
                Open-ended if ``None``.
            end (str | datetime, optional): Inclusive end of the window.
                Open-ended if ``None``.

        Returns:
            None: :attr:`ds` is replaced in place.
        """
        if start or end:
            self.ds = self.ds.sel(time=slice(start, end))

    # ---------------------------------------------------------
    # Step 4: Save methods (same API as CMIP)
    # ---------------------------------------------------------
    def save_netcdf(self, filename):
        """Write the loaded dataset to NetCDF.

        Args:
            filename (str): Destination path. The parent directory must exist.

        Returns:
            None
        """
        self.ds.to_netcdf(filename)

    def save_csv(self, filename):
        """Write the loaded dataset to CSV, one row per timestep.

        Args:
            filename (str): Destination path. The parent directory must exist.

        Returns:
            None
        """
        self.ds.to_dataframe().to_csv(filename)
