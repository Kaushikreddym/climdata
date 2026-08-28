try:
    import intake
    import intake_esm  # noqa: F401 — registers the esm driver with intake
    _INTAKE_AVAILABLE = True
except ImportError:
    _INTAKE_AVAILABLE = False

import xarray as xr
import pandas as pd
from omegaconf import DictConfig
import logging
from datetime import datetime

logger = logging.getLogger(__name__)


class CMIPCloud:
    """Stream CMIP6 model output from the Pangeo cloud catalogue.

    Reads analysis-ready Zarr stores from Google Cloud through ``intake-esm``,
    so nothing is downloaded up front — only the requested subset is transferred
    when the data is finally computed.

    The lifecycle is ``fetch`` → ``load`` → ``extract``, and unlike the
    file-based providers the order is strict: :meth:`extract` needs a loaded
    dataset and raises without one.

    The catalogue lookups (:meth:`open_cmip6_catalog`, :meth:`get_experiment_ids`,
    :meth:`get_source_ids`) are classmethods, so a caller browsing what is
    available — a model picker, say — can query them without building an
    extraction config first.

    Attributes:
        experiment_id (str): CMIP6 experiment, e.g. ``"historical"``, ``"ssp585"``.
        source_id (str): Model name, e.g. ``"GFDL-ESM4"``.
        table_id (str): MIP table setting the frequency, e.g. ``"day"``.
        variables (list[str]): Requested CF variable names.
        ds (xr.Dataset | None): The loaded dataset, set by :meth:`load`.
        col_subsets (list): Per-variable catalogue subsets, set by :meth:`fetch`.

    Example:
        >>> cmip = CMIPCloud(cfg)                             # doctest: +SKIP
        >>> cmip.fetch()                                      # doctest: +SKIP
        >>> cmip.load()                                       # doctest: +SKIP
        >>> ds = cmip.extract(point=(13.4, 52.5))             # doctest: +SKIP
    """

    # Catalogue lookups are classmethods so callers that only want to browse
    # the Pangeo catalogue (e.g. the GUI's model/experiment pickers) can query
    # them without building a full extraction config.
    _catalog = None   # process-wide cache; the catalogue is a ~10 s download

    @classmethod
    def open_cmip6_catalog(cls, refresh=False):
        """Open the Pangeo CMIP6 ESM datastore, caching it for the process.

        The catalogue is a multi-megabyte download taking roughly ten seconds,
        so the first call pays for it and every later call reuses the result.

        Args:
            refresh (bool): Re-download even if a cached catalogue exists. Use
                this to pick up newly published models within a long session.

        Returns:
            intake_esm.esm_datastore: The catalogue, whose ``.df`` is a pandas
            frame with one row per published Zarr store.
        """
        if cls._catalog is None or refresh:
            cls._catalog = intake.open_esm_datastore(
                "https://storage.googleapis.com/cmip6/pangeo-cmip6.json"
            )
        return cls._catalog

    @classmethod
    def get_experiment_ids(cls):
        """List the experiments climdata supports, from the live catalogue.

        Filtered deliberately to ``historical`` plus the canonical ``sspNNN``
        scenarios. The raw catalogue carries hundreds of DECK, CFMIP and
        diagnostic experiments whose time axes and forcings do not fit the
        historical/projection workflow this class is built for.

        Returns:
            list[str]: Sorted experiment IDs, e.g.
            ``["historical", "ssp119", "ssp126", ...]``.
        """
        import re
        col = cls.open_cmip6_catalog()

        experiments = sorted(col.df["experiment_id"].unique())

        pattern = re.compile(r"^ssp\d{3}$")  # ssp + exactly 3 digits

        experiments = [
            e for e in experiments
            if e == "historical" or pattern.match(e)
        ]

        return experiments

    @classmethod
    def get_source_ids(cls, experiment_id, table_id=None, variables=None):
        """List the CMIP6 models available for *experiment_id*.

        Args:
            experiment_id: CMIP6 experiment, e.g. ``"historical"`` or ``"ssp585"``.
            table_id: Restrict to a MIP table, e.g. ``"day"``. ``None`` keeps
                every frequency.
            variables: Restrict to models that publish **all** of these
                ``variable_id`` values (within *table_id* when given).
                ``None`` keeps every model.

        Returns:
            list[str]: Sorted model names.

        Notes:
            Filtering matters for callers that go on to extract: a model listed
            for an experiment does not necessarily publish the requested
            variables at the requested frequency, and ``fetch()`` would then
            fail. Passing the ``table_id``/``variables`` that will be used for
            the extraction yields only models that can actually serve it.
        """
        col = cls.open_cmip6_catalog()

        df = col.df[col.df["experiment_id"] == experiment_id]
        if len(df) == 0:
            raise ValueError(f"No data found for experiment_id={experiment_id}")

        if table_id:
            df = df[df["table_id"] == table_id]

        if variables:
            # Keep only models that carry every requested variable — the
            # intersection, not the union, since extract() needs them all.
            per_var = [
                set(df[df["variable_id"] == v]["source_id"].unique())
                for v in variables
            ]
            usable = set.intersection(*per_var) if per_var else set()
            df = df[df["source_id"].isin(usable)]

        sources = sorted(df["source_id"].unique())
        if not sources:
            raise ValueError(
                f"No models found for experiment_id={experiment_id}"
                + (f", table_id={table_id}" if table_id else "")
                + (f", variables={list(variables)}" if variables else "")
            )
        logger.info(f"{len(sources)} models found for experiment '{experiment_id}'")
        return sources


    

    def get_variables(self, *, experiment_id, source_id, table_id="day"):
        """List the climdata variables common to every requested model and experiment.

        Returns the *intersection*, not the union: only variables that every
        combination of ``experiment_id`` and ``source_id`` publishes. That is
        what a multi-model extraction needs, since one model missing one variable
        breaks the merge. The result is further narrowed to the six variables
        climdata's downstream index and bias-correction code understands
        (``tas``, ``tasmin``, ``tasmax``, ``pr``, ``hurs``, ``sfcWind``).

        Args:
            experiment_id (str | list[str]): One experiment or several.
            source_id (str | list[str]): One model or several.
            table_id (str | None): MIP table, e.g. ``"day"``. ``None`` searches
                every frequency. Defaults to ``"day"``.

        Returns:
            list[str]: Sorted CF variable names available across all combinations.

        Raises:
            ValueError: If no variable is common to every combination — including
                the case where a model/experiment pair publishes nothing at all.

        Example:
            >>> CMIPCloud(cfg).get_variables(                      # doctest: +SKIP
            ...     experiment_id="historical", source_id="GFDL-ESM4")
            ['hurs', 'pr', 'sfcWind', 'tas', 'tasmax', 'tasmin']
        """
        TARGET_VARS = {"tas", "tasmin", "tasmax", "pr", "hurs", "sfcWind"}
        
        col = self.open_cmip6_catalog()
        
        # Normalize to lists
        if isinstance(experiment_id, str):
            experiment_ids = [experiment_id]
        else:
            experiment_ids = list(experiment_id)
            
        if isinstance(source_id, str):
            source_ids = [source_id]
        else:
            source_ids = list(source_id)
        
        common_vars = None
        
        for exp in experiment_ids:
            for src in source_ids:
                query = dict(
                    experiment_id=[exp],
                    source_id=[src],
                )
                
                if table_id is not None:
                    query["table_id"] = [table_id]
                
                subset = col.search(**query)
                
                if len(subset.df) == 0:
                    continue
                
                available = set(subset.df["variable_id"].unique())
                selected = available & TARGET_VARS
                
                if common_vars is None:
                    common_vars = selected
                else:
                    common_vars = common_vars & selected
        
        if not common_vars:
            raise ValueError(
                f"No common variables found for "
                f"experiment_id={experiment_ids}, "
                f"source_id={source_ids}, "
                f"table_id={table_id}"
            )
        
        return sorted(common_vars)


    def __init__(self, cfg: DictConfig):
        """Bind a configuration and validate the requested period.

        Args:
            cfg (DictConfig): Configuration with ``experiment_id``, ``source_id``,
                ``table_id``, ``variables`` and ``time_range.start_date`` /
                ``.end_date``.

        Raises:
            ImportError: If ``intake`` and ``intake-esm`` are not installed.
            ValueError: If the time range cannot belong to the experiment — see
                :meth:`_validate_time_range`.
        """
        if not _INTAKE_AVAILABLE:
            raise ImportError(
                "CMIPCloud requires intake and intake-esm. "
                "Install with: pip install intake intake-esm"
            )
        # Directly read from flat config
        self.experiment_id = cfg.experiment_id
        self.source_id = cfg.source_id
        self.table_id = cfg.table_id
        self.variables = cfg.variables
        self.start_date = cfg.time_range.start_date
        self.end_date = cfg.time_range.end_date
        self.cfg = cfg
        self.col_subsets = []
        self.ds = None
        self.col = None
        self._validate_time_range()
    def _validate_time_range(self):
        """Check the requested period against the experiment's simulated period.

        CMIP6 historical runs cover 1850-2014 and SSP scenarios 2015-2100, so
        asking an SSP for 1990 returns nothing — a failure worth catching at
        construction rather than after a catalogue query. A range that merely
        overruns the period at one end warns and proceeds; one that misses it
        entirely raises. ``picontrol`` and unrecognised experiments are skipped.

        Returns:
            None

        Raises:
            ValueError: If the requested range lies wholly outside the
                experiment's period.
        """
        start_date = datetime.fromisoformat(self.cfg.time_range.start_date)
        end_date = datetime.fromisoformat(self.cfg.time_range.end_date)
        
        start_year = start_date.year
        end_year = end_date.year
        
        # Define valid periods for each experiment type
        if self.experiment_id == 'historical':
            valid_start = 1850
            valid_end = 2014
            period_name = "Historical"
        elif self.experiment_id.startswith('ssp'):
            valid_start = 2015
            valid_end = 2100
            period_name = f"SSP scenario ({self.experiment_id})"
        elif self.experiment_id == 'picontrol':
            # Pre-industrial control - typically long runs, less strict
            return
        else:
            # Unknown experiment, skip validation
            return
        
        # Check if requested period is outside valid range
        if end_year < valid_start or start_year > valid_end:
            raise ValueError(
                f"❌ Time range mismatch for experiment '{self.experiment_id}'!\n"
                f"   Requested: {start_year}-{end_year}\n"
                f"   Valid period for {period_name}: {valid_start}-{valid_end}\n"
                f"   \n"
                f"   Hint: Use 'historical' for years 1850-2014, and SSP scenarios (ssp126, ssp370, ssp585) for 2015-2100."
            )
        
        # Warn if requested period extends beyond valid range
        if start_year < valid_start or end_year > valid_end:
            print(f"⚠️  Warning: Requested time range {start_year}-{end_year} extends beyond")
            print(f"   the typical {period_name} period ({valid_start}-{valid_end}).")
            print(f"   Data availability may be limited.")
    def fetch(self):
        """Collect intake catalog subsets for each variable.

        Raises:
            ValueError: If the model/experiment/table combination does not
                publish every requested variable. Partial coverage is treated
                as an error rather than silently dropped, because downstream
                extraction indexes the dataset by every name in
                ``cfg.variables`` and would otherwise fail with a bare
                ``KeyError`` far from the cause.
        """
        col = self.open_cmip6_catalog()
        self.col_subsets = []
        missing = []

        for var in self.variables:
            query = dict(
                experiment_id=[self.experiment_id],
                source_id=self.source_id,
                table_id=self.table_id,
                variable_id=var,
            )
            col_subset = col.search(require_all_on=["source_id"], **query)

            if len(col_subset.df) == 0:
                missing.append(var)
                continue

            self.col_subsets.append(col_subset)
            self.col = col

        if missing:
            # Report what this model *does* offer at this frequency so the
            # caller can pick a workable variable set or a different model.
            available = sorted(
                col.df[
                    (col.df["experiment_id"] == self.experiment_id)
                    & (col.df["source_id"] == self.source_id)
                    & (col.df["table_id"] == self.table_id)
                ]["variable_id"].unique()
            )
            raise ValueError(
                f"CMIP6 model '{self.source_id}' has no "
                f"{self.table_id}-frequency data for {missing} "
                f"under experiment '{self.experiment_id}'.\n"
                f"   Requested variables : {list(self.variables)}\n"
                f"   Available for this model/table: {available or 'none'}\n"
                f"   Hint: choose a different model, or reduce the requested "
                f"variables to the available ones."
            )

        if not self.col_subsets:
            raise ValueError(
                f"No matching CMIP6 data found for: "
                f"experiment_id={self.experiment_id}, "
                f"source_id={self.source_id}, "
                f"table_id={self.table_id}, "
                f"variables={self.variables}"
            )

        return self.col_subsets
    def convert_to_noleap(self, ds):
        """Rewrite a CMIP time axis as day-floored pandas timestamps.

        CMIP6 models use whichever calendar their modelling centre chose —
        ``noleap``, ``360_day``, ``proleptic_gregorian`` — and stamp daily data
        at 12:00. Both properties block merging across models and joining against
        observations. Converting to ``pandas.Timestamp`` at midnight makes the
        axes comparable.

        This is lossy where the source calendar is not Gregorian: 360-day
        timestamps become real dates, which can produce duplicate or missing days
        across a long record. It is applied at the end of :meth:`extract`, once
        the data has been subset.

        Args:
            ds (xr.Dataset): Dataset with a CMIP time coordinate.

        Returns:
            xr.Dataset: The dataset with a ``pandas``-backed daily time axis, or
            unchanged if it has no ``time`` coordinate.
        """
        if "time" not in ds.coords:
            return ds
        
        t = ds.indexes["time"]  # can be pandas.DatetimeIndex or CFTimeIndex
        new_times = []

        for ti in t:
            year, month, day = ti.year, ti.month, ti.day
            
            # Floor to day: ignore hour, minute, second
            # Convert all time types to pandas Timestamp for compatibility
            new_times.append(pd.Timestamp(year=year, month=month, day=day))
        
        ds = ds.assign_coords(time=("time", new_times))
        return ds

    def load(self):
        """Open the Zarr store behind each catalogue subset and merge them.

        Stores are opened lazily over HTTPS, so this is fast and transfers almost
        nothing; the data moves when the result is computed. The merge uses
        ``compat="override"``, taking coordinate values from the first dataset
        where models disagree in the last bits of their grid coordinates.

        Only the first store per variable is opened, so a variable split across
        several ensemble members contributes one member.

        Returns:
            xr.Dataset | None: The merged dataset, also stored on :attr:`ds`.
            ``None`` if :meth:`fetch` collected nothing.
        """
        datasets = []
        for col_subset in self.col_subsets:
            zstore_path = col_subset.df.zstore.values[0].replace(
                "gs:/", "https://storage.googleapis.com"
            )
            ds_var = xr.open_zarr(zstore_path)
            datasets.append(ds_var)
        if datasets:
            self.ds = xr.merge(datasets,compat='override')
        else:
            self.ds = None

        return self.ds

    def extract(self, *, point=None, box=None, shapefile=None, buffer_km=0.0):
        """Subset the loaded dataset in time and space.

        Always clips to the configured time range first, then applies exactly one
        spatial selector. Afterwards the model name is added as a ``source_id``
        dimension — so several models can be concatenated — and the time axis is
        normalised by :meth:`convert_to_noleap`.

        CMIP6 longitudes run 0-360, so ``point`` and ``box`` coordinates must be
        given in that convention, not as negative western longitudes.

        Args:
            point (tuple[float, float], optional): ``(lon, lat)`` in degrees,
                longitude first.
            box (dict, optional): Bounding box with keys ``lon_min``, ``lon_max``,
                ``lat_min``, ``lat_max``.
            shapefile (str | geopandas.GeoDataFrame, optional): Polygon(s) to clip
                to, as a path or an in-memory frame.
            buffer_km (float): For ``point``, half-width of a box selected around
                it instead of the nearest cell; for ``shapefile``, a dilation
                applied to each geometry. Converted at a flat 111 km per degree.
                Defaults to ``0.0``.

        Returns:
            xr.Dataset: The subset, also stored on :attr:`ds`.

        Raises:
            ValueError: If :meth:`load` has not run, or if none of ``point``,
                ``box`` or ``shapefile`` was given.
        """
        import geopandas as gpd
        from shapely.geometry import mapping

        if self.ds is None:
            raise ValueError("No dataset loaded. Call `load()` first.")
        
        self._subset_time(self.start_date, self.end_date) 
        
        ds = self.ds
        if point is not None:
            lon, lat = point
            if buffer_km > 0:
                buffer_deg = buffer_km / 111
                ds_subset = ds.sel(
                    lon=slice(lon - buffer_deg, lon + buffer_deg),
                    lat=slice(lat - buffer_deg, lat + buffer_deg),
                )
            else:
                ds_subset = ds.sel(lon=lon, lat=lat, method="nearest")

        elif box is not None:
            ds_subset = ds.sel(
                lon=slice(box["lon_min"], box["lon_max"]),
                lat=slice(box["lat_min"], box["lat_max"]),
            )

        elif shapefile is not None:
            if isinstance(shapefile, str):
                gdf = gpd.read_file(shapefile)
            else:
                gdf = shapefile
            if buffer_km > 0:
                gdf = gdf.to_crs(epsg=3857)
                gdf["geometry"] = gdf.buffer(buffer_km * 1000)
                gdf = gdf.to_crs(epsg=4326)
            geom = [mapping(g) for g in gdf.geometry]
            import rioxarray  # noqa: F401 — registers the .rio accessor

            ds = ds.rio.write_crs("EPSG:4326", inplace=False)
            ds_subset = ds.rio.clip(geom, gdf.crs, drop=True)

        else:
            raise ValueError("Must provide either point, box, or shapefile.")
        self.ds = ds_subset
        self.ds = self.ds.assign_coords(source_id=self.source_id)
        self.ds = self.ds.expand_dims("source_id")
        self.ds = self.convert_to_noleap(self.ds)
        return ds_subset

    def _subset_time(self, start_date, end_date):
        """Clip :attr:`ds` to a time range, in place.

        Args:
            start_date (str): Inclusive start, ISO format.
            end_date (str): Inclusive end, ISO format.

        Returns:
            xr.Dataset | None: The clipped dataset, or ``None`` if nothing is loaded.
        """
        if self.ds is None:
            return None
        ds_time = self.ds.sel(time=slice(start_date, end_date))
        self.ds = ds_time
        return ds_time

    def save_netcdf(self, filename):
        """Write the dataset to NetCDF, doing nothing if none is loaded.

        Clears the time encoding inherited from the Zarr store first: it names a
        calendar that no longer matches the axis rewritten by
        :meth:`convert_to_noleap`, and the NetCDF writer would otherwise fail on
        the contradiction.

        Args:
            filename (str): Destination path. The parent directory must exist.

        Returns:
            None
        """
        if self.ds is not None:
            if "time" in self.ds.variables:
                self.ds["time"].encoding.clear()
            self.ds.to_netcdf(filename)
            # print(f"Saved NetCDF to {filename}")

    def save_zarr(self, store_path):
        """Write the dataset to a Zarr store, overwriting any existing one.

        Args:
            store_path (str): Destination store path.

        Returns:
            None
        """
        if self.ds is not None:
            self.ds.to_zarr(store_path, mode="w")
            print(f"Saved Zarr to {store_path}")

    def _format(self, df):
        """Reshape a wide frame into climdata's long output form.

        Melts the variable columns into ``variable``/``value`` pairs, attaches
        each variable's units, and stamps the model, experiment and MIP table so
        rows stay identifiable once several models are concatenated.

        Args:
            df (pd.DataFrame): Wide frame from ``ds.to_dataframe().reset_index()``.

        Returns:
            pd.DataFrame: Long frame with columns ``source``, ``experiment``,
            ``table``, ``time``, ``lat``, ``lon``, ``variable``, ``value``,
            ``units`` — those of them that are present.
        """
        value_vars = [v for v in self.variables if v in df.columns]
        id_vars = [c for c in df.columns if c not in value_vars]

        df_long = df.melt(
            id_vars=id_vars,
            value_vars=value_vars,
            var_name="variable",
            value_name="value",
        )

        df_long["units"] = df_long["variable"].map(
            lambda v: self.ds[v].attrs.get("units", "unknown")
            if v in self.ds.data_vars
            else "unknown"
        )

        df_long["source"] = self.source_id
        df_long["experiment"] = self.experiment_id
        df_long["table"] = self.table_id

        cols = [
            "source",
            "experiment",
            "table",
            "time",
            "lat",
            "lon",
            "variable",
            "value",
            "units",
        ]
        df_long = df_long[[c for c in cols if c in df_long.columns]]

        return df_long

    def save_csv(self, filename):
        """Write the dataset to CSV in long form, doing nothing if none is loaded.

        Args:
            filename (str): Destination path. The parent directory must exist.

        Returns:
            None
        """
        if self.ds is not None:
            df = self.ds.to_dataframe().reset_index()
            df = self._format(df)
            df.to_csv(filename, index=False)
            # print(f"Saved CSV to {filename}")