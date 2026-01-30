#!/usr/bin/env python3
"""
Test script to demonstrate auto-discovery of member_id and grid_label
for NEX-GDDP-CMIP6 data
"""

from omegaconf import OmegaConf
from climdata.datasets.NEXGDDP import NEXGDDP

# Create a minimal configuration
cfg = OmegaConf.create({
    "source_id": "GFDL-ESM4",  # Model name
    "experiment_id": "historical",  # Experiment
    "variables": ["tasmax", "tasmin", "pr"],  # Variables to download
    "time_range": {
        "start_date": "2010-01-01",
        "end_date": "2010-12-31"
    },
    "data_dir": "./outputs/nexgddp_test"
})

print("=" * 80)
print("Testing NEX-GDDP Auto-Discovery of Metadata")
print("=" * 80)

# Initialize NEXGDDP - this will auto-discover member_id and grid_label
print(f"\nInitializing NEXGDDP for {cfg.source_id} / {cfg.experiment_id}...")
print("-" * 80)

nexgddp = NEXGDDP(cfg)

print("\n" + "=" * 80)
print("Discovered Metadata:")
print("=" * 80)
print(f"Model (source_id):     {nexgddp.source_id}")
print(f"Experiment:            {nexgddp.experiment_id}")
print(f"Realization (member_id): {nexgddp.member_id}")
print(f"Grid Label:            {nexgddp.grid_label}")

# Get all available member_ids for this model/experiment
print("\n" + "=" * 80)
print("All Available Realizations:")
print("=" * 80)
member_ids = nexgddp.get_member_ids()
if member_ids:
    for mid in member_ids:
        print(f"  - {mid}")
else:
    print("  (Could not retrieve member IDs)")

# Show sample URL construction
print("\n" + "=" * 80)
print("Sample Download URL:")
print("=" * 80)
sample_url, sample_filename = nexgddp._construct_download_url("tasmax", 2010)
print(f"URL:      {sample_url}")
print(f"Filename: {sample_filename}")

print("\n" + "=" * 80)
print("Test Complete!")
print("=" * 80)
