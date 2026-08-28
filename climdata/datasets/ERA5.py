# SPDX-FileCopyrightText: Copyright (c) 2023 - 2024 NVIDIA CORPORATION & AFFILIATES.
# SPDX-FileCopyrightText: All rights reserved.
# SPDX-License-Identifier: Apache-2.0
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os
import tempfile
try:
    import cdsapi
    _CDSAPI_AVAILABLE = True
except ImportError:
    _CDSAPI_AVAILABLE = False
import xarray as xr
import datetime
import json
import dask
import calendar
from dask.diagnostics import ProgressBar
from typing import List, Tuple, Union
import urllib3
import numpy as np
import fsspec

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)


class ERA5Mirror:
    """Build and maintain a local Zarr mirror of ERA5 from the Copernicus CDS.

    Unlike the other providers this is a *mirroring* tool rather than an
    extractor: it downloads whole global months and appends them to one Zarr
    store per variable, to be read afterwards with ``xarray.open_zarr``. It also
    takes a storage path rather than a Hydra config, so it is not driven by
    :class:`~climdata.utils.wrapper_workflow.ClimateExtractor`.

    Progress is journalled in a ``metadata.json`` beside the stores, so an
    interrupted mirror resumes at the first month it has not yet written instead
    of re-downloading. Months are fetched concurrently through Dask.

    Requires a Copernicus CDS account and a ``~/.cdsapirc`` key file.

    Attributes:
        base_path (str): Directory holding the Zarr stores and ``metadata.json``.
        fs (fsspec.AbstractFileSystem): Filesystem the stores live on — local by
            default, but any fsspec target (S3, GCS) works.
        metadata (dict): The journal of downloaded chunks.

    Example:
        >>> import datetime
        >>> mirror = ERA5Mirror("./era5")                        # doctest: +SKIP
        >>> paths = mirror.download(                             # doctest: +SKIP
        ...     ["2m_temperature", ("temperature", 500)],
        ...     (datetime.date(2020, 1, 1), datetime.date(2020, 3, 1)),
        ... )
    """

    def __init__(self, base_path: str, fs: fsspec.AbstractFileSystem = None):
        """Open (or create) a mirror directory and read its journal.

        Args:
            base_path (str): Directory for the Zarr stores. Created if absent.
            fs (fsspec.AbstractFileSystem, optional): Filesystem to use.
                Defaults to the local filesystem.

        Raises:
            ImportError: If ``cdsapi`` is not installed.
        """
        if not _CDSAPI_AVAILABLE:
            raise ImportError(
                "ERA5 requires cdsapi. "
                "Install with: pip install cdsapi"
            )
        # Get parameters
        self.base_path = base_path
        if fs is None:
            fs = fsspec.filesystem("file")
        self.fs = fs

        # Create the base path if it doesn't exist
        if not self.fs.exists(self.base_path):
            self.fs.makedirs(self.base_path)

        # Create metadata that will be used to track which chunks have been downloaded
        self.metadata_file = os.path.join(self.base_path, "metadata.json")
        self.metadata = self.get_metadata()

    def get_metadata(self):
        """Read the download journal from ``metadata.json``.

        Returns:
            dict: ``{"chunks": [...]}``. An absent or corrupt file yields an
            empty journal rather than raising, so a truncated write from an
            interrupted run costs a re-download rather than a crash.
        """
        if self.fs.exists(self.metadata_file):
            with self.fs.open(self.metadata_file, "r") as f:
                try:
                    metadata = json.load(f)
                except json.decoder.JSONDecodeError:
                    metadata = {"chunks": []}
        else:
            metadata = {"chunks": []}
        return metadata

    def save_metadata(self):
        """Write the download journal back to ``metadata.json``.

        Returns:
            None
        """
        with self.fs.open(self.metadata_file, "w") as f:
            json.dump(self.metadata, f)

    def chunk_exists(self, variable, year, month, pressure_level):
        """Test whether a month has already been mirrored.

        Args:
            variable (str): CDS variable name.
            year (int): Calendar year.
            month (int): Calendar month, 1-12.
            pressure_level (int | None): Pressure level in hPa, or ``None`` for
                single-level data.

        Returns:
            bool: ``True`` if the journal records this chunk.
        """
        for chunk in self.metadata["chunks"]:
            if (
                chunk["variable"] == variable
                and chunk["year"] == year
                and chunk["month"] == month
                and chunk["pressure_level"] == pressure_level
            ):
                return True
        return False

    def download_chunk(
        self,
        variable: str,
        year: int,
        month: int,
        pressure_level: int = None,
    ):
        """Download one global month of daily-mean ERA5 data.

        Requests the CDS daily-statistics product, so what arrives is already
        aggregated to daily means from 6-hourly steps rather than raw hourly
        data. The NetCDF lands in a temporary directory and is opened before
        that directory is removed, so the returned dataset must be consumed or
        loaded before it goes out of scope.

        Args:
            variable (str): CDS variable name, e.g. ``"2m_temperature"`` or
                ``"total_precipitation"``.
            year (int): Calendar year.
            month (int): Calendar month, 1-12.
            pressure_level (int, optional): Pressure level in hPa. ``None``, the
                default, requests the single-level product.

        Returns:
            xr.Dataset: The downloaded month.

        Raises:
            Exception: Whatever ``cdsapi`` raises for a rejected request — an
                unknown variable name, or a missing/invalid ``~/.cdsapirc``.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            # Get all days in the month
            days_in_month = calendar.monthrange(year, month)[1]

            # Make tmpfile to store the data
            output_file = os.path.join(
                tmpdir,
                f"{variable}_{year}_{month:02d}_{str(pressure_level)}.nc",
            )

            # start the CDS API client (maybe need to move this outside the loop?)
            c = cdsapi.Client(quiet=True)

            # Setup the request parameters
            request_params = {
                "product_type": "reanalysis",
                "variable": [variable],
                "year": str(year),
                "month": str(month),
                "day": [f"{day:02d}" for day in range(1, days_in_month + 1)],
                "time_zone": "utc+00:00",
                "frequency": "6_hourly",
                "daily_statistic": "daily_mean",
                "data_format": "netcdf"
            }
            if pressure_level:
                request_params["pressure_level"] = [str(pressure_level)]
                dataset_name = "derived-era5-pressure-levels-daily-statistics"
            else:
                dataset_name = "derived-era5-single-levels-daily-statistics"

            # Download the data
            c.retrieve(
                dataset_name,
                request_params,
                output_file,
            )

            # Open the downloaded data
            ds = xr.open_dataset(output_file)
        return ds

    def variable_to_zarr_name(self, variable: str, pressure_level: int = None):
        """Build the Zarr store path for a variable.

        Args:
            variable (str): CDS variable name.
            pressure_level (int, optional): Pressure level in hPa. When given it
                is folded into the name, so each level gets its own store.

        Returns:
            str: Path of the form ``<base_path>/<variable>.zarr`` or
            ``<base_path>/<variable>_pressure_level_<n>.zarr``.
        """
        # create zarr path for variable
        zarr_path = f"{self.base_path}/{variable}"
        if pressure_level:
            zarr_path += f"_pressure_level_{pressure_level}"
        zarr_path += ".zarr"
        return zarr_path

    def download_and_upload_chunk(
        self,
        variable: str,
        year: int,
        month: int,
        pressure_level: int = None,
    ):
        """Download one month and append it to the variable's Zarr store.

        Creates the store on the first month and appends along ``time``
        thereafter. Chunking is one timestep by the full global grid, which suits
        reading a few dates over a wide area — the common ERA5 access pattern —
        rather than long time series at one point.

        The journal is updated only after a successful write, so a crash
        mid-append leaves the month marked as missing and it is retried.

        Args:
            variable (str): CDS variable name.
            year (int): Calendar year.
            month (int): Calendar month, 1-12.
            pressure_level (int, optional): Pressure level in hPa.

        Returns:
            None
        """

        # Download the data
        ds = self.download_chunk(variable, year, month, pressure_level)
        if "valid_time" in ds.dims:
            ds = ds.rename({"valid_time": "time"})

        # Create the Zarr path
        zarr_path = self.variable_to_zarr_name(variable, pressure_level)

        # Specify the chunking options
        chunking = {"time": 1, "latitude": 721, "longitude": 1440}
        if "level" in ds.dims:
            chunking["level"] = 1

        # Re-chunk the dataset
        ds = ds.chunk(chunking)

        # Check if the Zarr dataset exists
        if self.fs.exists(zarr_path):
            mode = "a"
            append_dim = "time"
            create = False
        else:
            mode = "w"
            append_dim = None
            create = True

        # Upload the data to the Zarr dataset
        mapper = self.fs.get_mapper(zarr_path, create=create)
        ds.to_zarr(mapper, mode=mode, consolidated=True, append_dim=append_dim)

        # Update the metadata
        self.metadata["chunks"].append(
            {
                "variable": variable,
                "year": year,
                "month": month,
                "pressure_level": pressure_level,
            }
        )
        self.save_metadata()

    def download(
        self,
        variables: List[Union[str, Tuple[str, int]]],
        date_range: Tuple[datetime.date, datetime.date],
    ):
        """Mirror a set of variables over a date range, month by month.

        Walks the range one calendar month at a time, skipping months the journal
        already records and submitting the rest to Dask so several download
        concurrently. Both endpoints are rounded down to the first of the month,
        so any month the range touches is fetched whole.

        Every resulting store is verified afterwards to have an evenly spaced
        time axis — the symptom of a download that failed partway through.

        Args:
            variables (list[str | tuple[str, int]]): Variables to mirror. A bare
                string requests single-level data; a ``(variable, level)`` tuple
                requests one pressure level.
            date_range (tuple[datetime.date, datetime.date]): Inclusive start and
                end dates.

        Returns:
            list[str]: Path of the Zarr store for each requested variable.

        Raises:
            AssertionError: If a store ends with an unevenly spaced time axis.
                Delete that store and re-run — the message names the path.
        """

        start_date, end_date = date_range

        # Reformat the variables list so all elements are tuples
        reformated_variables = []
        for variable in variables:
            if isinstance(variable, str):
                reformated_variables.append(tuple([variable, None]))
            else:
                reformated_variables.append(variable)

        # Start Downloading
        with ProgressBar():
            # Round dates to months
            current_date = start_date.replace(day=1)
            end_date = end_date.replace(day=1)

            while current_date <= end_date:
                # Create a list of tasks to download the data
                tasks = []
                for variable, pressure_level in reformated_variables:
                    if not self.chunk_exists(
                        variable,
                        current_date.year,
                        current_date.month,
                        pressure_level,
                    ):
                        task = dask.delayed(self.download_and_upload_chunk)(
                            variable,
                            current_date.year,
                            current_date.month,
                            pressure_level,
                        )
                        tasks.append(task)
                    else:
                        print(
                            f"Chunk for {variable} {pressure_level} {current_date.year}-{current_date.month} already exists. Skipping."
                        )

                # Execute the tasks with Dask
                print(f"Downloading data for {current_date.year}-{current_date.month}")
                if tasks:
                    dask.compute(*tasks)

                # Update the metadata
                self.save_metadata()

                # Update the current date
                days_in_month = calendar.monthrange(
                    year=current_date.year, month=current_date.month
                )[1]
                current_date += datetime.timedelta(days=days_in_month)

        # Return the Zarr paths
        zarr_paths = []
        for variable, pressure_level in reformated_variables:
            zarr_path = self.variable_to_zarr_name(variable, pressure_level)
            zarr_paths.append(zarr_path)

        # Check that Zarr arrays have correct dt for time dimension
        for zarr_path in zarr_paths:
            ds = xr.open_zarr(zarr_path)
            time_stamps = ds.time.values
            dt = time_stamps[1:] - time_stamps[:-1]
            assert np.all(dt == dt[0]), (
                f"Zarr array {zarr_path} has incorrect dt for time dimension. An error may have occurred during download. Please delete the Zarr array and try again."
            )

        return zarr_paths