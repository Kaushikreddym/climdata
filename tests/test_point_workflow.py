import json
from climdata import ClimData
import logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
# -----------------------------
# Step 1: Define the area of interest (AOI)
# -----------------------------
# The AOI is a single point. In GeoJSON format, the coordinates are [longitude, latitude].
geojson = {
  "type": "FeatureCollection",
  "features": [
    {
      "type": "Feature",
      "properties": {},
      "geometry": {
        "coordinates": [
          24.246667038198012,  # longitude
          12.891982026993958   # latitude
        ],
        "type": "Point"
      }
    }
  ]
}


# -----------------------------
# Step 2: Define configuration overrides
# -----------------------------
# Overrides are strings used by Hydra to modify default configurations at runtime.
overrides = [
    "dataset=cmip",  # Choose the MSWX dataset for extraction
    f"aoi='{json.dumps(geojson)}'",  # Set the AOI as the point defined above
    f"time_range.start_date=2014-12-01",  # Start date for data extraction
    f"time_range.end_date=2014-12-31",    # End date for data extraction
    "variables=[tasmin,tasmax,pr]",       # Variables to extract: min/max temp and precipitation
    "data_dir=/beegfs/muduchuru/data",    # Local directory to store raw/intermediate files
    # "dsinfo.mswx.params.google_service_account=./.climdata_conf/service.json",  # optional . required for MSWS data download
    "index=tn10p",  # Climate extreme index to calculate
    "impute=BRITS"
]

# -----------------------------
# Step 3: Define the workflow sequence
# -----------------------------
seq = ["extract", "impute", "calc_index", "to_nc"]

# -----------------------------
# Step 4: Initialize the ClimData extractor
# -----------------------------
extractor = ClimData(overrides=overrides)

# -----------------------------
# Step 5: Run the Multi-Step workflow
# -----------------------------
result = extractor.run_workflow(
    actions=seq,
)