"""
Test all dataset extractors for point, box, and shapefile AOI modes.
Time range: January 2010 (2010-01-01 to 2010-01-31).

Each dataset is parameterized with the three AOI modes. Tests only verify
that extraction completes and returns a non-empty xarray.Dataset. Datasets
that require credentials, large downloads, or are unavailable in the test
environment are skipped via the SKIP_* flags below.

Run a single dataset:
    pytest tests/test_datasets_extraction.py -k "hyras" -v

Run all enabled datasets:
    pytest tests/test_datasets_extraction.py -v
"""

import json
import os
import unittest
import xarray as xr

# Skip the entire module when running in GitHub Actions CI
# (tests require HPC data paths and external API credentials)
if os.getenv("CI"):
    raise unittest.SkipTest(
        "Integration tests require HPC data paths and credentials — skipped in CI"
    )
import geopandas as gpd
from shapely.geometry import box as shapely_box
from climdata import ClimData

# ---------------------------------------------------------------------------
# Skip flags — set to True to skip datasets that need special setup
# ---------------------------------------------------------------------------
SKIP_MSWX      = True   # requires Google Drive service account
SKIP_ERA5      = True   # requires CDS API key
SKIP_CMIP      = False
SKIP_DWD       = False
SKIP_HYRAS     = False
SKIP_HOSTRADA  = False
SKIP_W5E5      = False
SKIP_CMIP_W5E5 = False
SKIP_NEXGDDP   = False
SKIP_POWER     = False
SKIP_AGRI      = False

DATA_DIR = "/beegfs/muduchuru/data"
START    = "2010-01-01"
END      = "2010-01-31"
AOI_MODES = ["point", "box", "shapefile"]

# ---------------------------------------------------------------------------
# AOI definitions
# ---------------------------------------------------------------------------
POINT_AOI = json.dumps({
    "type": "FeatureCollection",
    "features": [{
        "type": "Feature",
        "properties": {},
        "geometry": {
            "type": "Point",
            "coordinates": [10.0, 52.0]   # lon, lat — central Germany
        }
    }]
})

BOX_AOI = json.dumps({
    "type": "FeatureCollection",
    "features": [{
        "type": "Feature",
        "properties": {},
        "geometry": {
            "type": "Polygon",
            "coordinates": [[
                [9.0, 51.0],
                [11.0, 51.0],
                [11.0, 53.0],
                [9.0, 53.0],
                [9.0, 51.0]
            ]]
        }
    }]
})


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture(scope="session")
def shapefile_path(tmp_path_factory):
    """Create a small in-memory shapefile covering the same region as BOX_AOI."""
    tmp = tmp_path_factory.mktemp("shp")
    gdf = gpd.GeoDataFrame(
        {"id": [1]},
        geometry=[shapely_box(9.0, 51.0, 11.0, 53.0)],
        crs="EPSG:4326"
    )
    path = str(tmp / "test_region.shp")
    gdf.to_file(path)
    return path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def aoi_override(mode: str, shp: str = "") -> str:
    if mode == "point":
        return f"aoi='{POINT_AOI}'"
    elif mode == "box":
        return f"aoi='{BOX_AOI}'"
    elif mode == "shapefile":
        return f"shapefile={shp}"
    raise ValueError(f"Unknown mode: {mode}")


def _run(overrides: list) -> xr.Dataset:
    extractor = ClimData(overrides=overrides)
    ds = extractor.extract()
    assert ds is not None, "extract() returned None"
    assert isinstance(ds, xr.Dataset), f"Expected xr.Dataset, got {type(ds)}"
    assert len(ds.data_vars) > 0, "Dataset has no data variables"
    for var in ds.data_vars:
        assert ds[var].size > 0, f"Variable {var} is empty"
    return ds


# ---------------------------------------------------------------------------
# DWD
# ---------------------------------------------------------------------------
@pytest.mark.skipif(SKIP_DWD, reason="DWD skipped")
@pytest.mark.parametrize("mode", AOI_MODES)
def test_dwd(mode, shapefile_path):
    overrides = [
        "dataset=DWD",
        aoi_override(mode, shapefile_path),
        f"time_range.start_date={START}",
        f"time_range.end_date={END}",
        "variables=[tasmin,tasmax,pr]",
        f"data_dir={DATA_DIR}",
    ]
    _run(overrides)


# ---------------------------------------------------------------------------
# MSWX
# ---------------------------------------------------------------------------
@pytest.mark.skipif(SKIP_MSWX, reason="MSWX requires Google Drive credentials")
@pytest.mark.parametrize("mode", AOI_MODES)
def test_mswx(mode, shapefile_path):
    overrides = [
        "dataset=MSWX",
        aoi_override(mode, shapefile_path),
        f"time_range.start_date={START}",
        f"time_range.end_date={END}",
        "variables=[tasmin,tasmax,pr]",
        f"data_dir={DATA_DIR}",
    ]
    _run(overrides)


# ---------------------------------------------------------------------------
# ERA5
# ---------------------------------------------------------------------------
@pytest.mark.skipif(SKIP_ERA5, reason="ERA5 requires CDS API key")
@pytest.mark.parametrize("mode", AOI_MODES)
def test_era5(mode, shapefile_path):
    overrides = [
        "dataset=ERA5",
        aoi_override(mode, shapefile_path),
        f"time_range.start_date={START}",
        f"time_range.end_date={END}",
        "variables=[tasmin,tasmax,pr]",
        f"data_dir={DATA_DIR}",
    ]
    _run(overrides)


# ---------------------------------------------------------------------------
# CMIP (cloud)
# ---------------------------------------------------------------------------
@pytest.mark.skipif(SKIP_CMIP, reason="CMIP skipped")
@pytest.mark.parametrize("mode", AOI_MODES)
def test_cmip(mode, shapefile_path):
    overrides = [
        "dataset=CMIP",
        aoi_override(mode, shapefile_path),
        f"time_range.start_date={START}",
        f"time_range.end_date={END}",
        "variables=[tasmin,tasmax,pr]",
        "experiment_id=historical",
        "source_id=MIROC6",
        f"data_dir={DATA_DIR}",
    ]
    _run(overrides)


# ---------------------------------------------------------------------------
# HYRAS
# ---------------------------------------------------------------------------
@pytest.mark.skipif(SKIP_HYRAS, reason="HYRAS skipped")
@pytest.mark.parametrize("mode", AOI_MODES)
def test_hyras(mode, shapefile_path):
    overrides = [
        "dataset=HYRAS",
        aoi_override(mode, shapefile_path),
        f"time_range.start_date={START}",
        f"time_range.end_date={END}",
        "variables=[tasmin,tasmax,pr]",
        f"data_dir={DATA_DIR}",
    ]
    _run(overrides)


# ---------------------------------------------------------------------------
# HOSTRADA
# ---------------------------------------------------------------------------
@pytest.mark.skipif(SKIP_HOSTRADA, reason="HOSTRADA skipped")
@pytest.mark.parametrize("mode", AOI_MODES)
def test_hostrada(mode, shapefile_path):
    overrides = [
        "dataset=HOSTRADA",
        aoi_override(mode, shapefile_path),
        f"time_range.start_date={START}",
        f"time_range.end_date={END}",
        "variables=[tasmin,tasmax,pr]",
        f"data_dir={DATA_DIR}",
    ]
    _run(overrides)


# ---------------------------------------------------------------------------
# W5E5
# ---------------------------------------------------------------------------
@pytest.mark.skipif(SKIP_W5E5, reason="W5E5 skipped")
@pytest.mark.parametrize("mode", AOI_MODES)
def test_w5e5(mode, shapefile_path):
    overrides = [
        "dataset=W5E5",
        aoi_override(mode, shapefile_path),
        f"time_range.start_date={START}",
        f"time_range.end_date={END}",
        "variables=[tasmin,tasmax,pr]",
        f"data_dir={DATA_DIR}",
    ]
    _run(overrides)


# ---------------------------------------------------------------------------
# CMIP_W5E5
# ---------------------------------------------------------------------------
@pytest.mark.skipif(SKIP_CMIP_W5E5, reason="CMIP_W5E5 skipped")
@pytest.mark.parametrize("mode", AOI_MODES)
def test_cmip_w5e5(mode, shapefile_path):
    overrides = [
        "dataset=CMIP_W5E5",
        aoi_override(mode, shapefile_path),
        f"time_range.start_date={START}",
        f"time_range.end_date={END}",
        "variables=[tasmin,tasmax,pr]",
        "experiment_id=historical",
        "source_id=MIROC6",
        f"data_dir={DATA_DIR}",
    ]
    _run(overrides)


# ---------------------------------------------------------------------------
# NEXGDDP
# ---------------------------------------------------------------------------
@pytest.mark.skipif(SKIP_NEXGDDP, reason="NEXGDDP skipped")
@pytest.mark.parametrize("mode", AOI_MODES)
def test_nexgddp(mode, shapefile_path):
    overrides = [
        "dataset=NEXGDDP",
        aoi_override(mode, shapefile_path),
        f"time_range.start_date={START}",
        f"time_range.end_date={END}",
        "variables=[tasmin,tasmax,pr]",
        "experiment_id=historical",
        "source_id=MIROC6",
        f"data_dir={DATA_DIR}",
    ]
    _run(overrides)


# ---------------------------------------------------------------------------
# NASA POWER
# ---------------------------------------------------------------------------
@pytest.mark.skipif(SKIP_POWER, reason="POWER skipped")
@pytest.mark.parametrize("mode", AOI_MODES)
def test_power(mode, shapefile_path):
    overrides = [
        "dataset=POWER",
        aoi_override(mode, shapefile_path),
        f"time_range.start_date={START}",
        f"time_range.end_date={END}",
        "variables=[tasmin,tasmax,pr]",
        f"data_dir={DATA_DIR}",
    ]
    _run(overrides)


# ---------------------------------------------------------------------------
# AGRI_ISIMIP
# ---------------------------------------------------------------------------
@pytest.mark.skipif(SKIP_AGRI, reason="AGRI_ISIMIP skipped")
@pytest.mark.parametrize("mode", AOI_MODES)
def test_agri_isimip(mode, shapefile_path):
    overrides = [
        "dataset=AGRI_ISIMIP",
        aoi_override(mode, shapefile_path),
        f"time_range.start_date={START}",
        f"time_range.end_date={END}",
        "variables=[tasmin,tasmax,pr]",
        f"data_dir={DATA_DIR}",
    ]
    _run(overrides)
