# Installation

## Stable release

It's recommended to create and activate a conda environment first, then install via pip:

```bash
# create & activate conda environment (recommended)
conda create -n climdata python=3.11 -y
conda activate climdata

# install climdata from PyPI
pip install climdata
```

This is the preferred method to install climdata, as it will always install the most recent stable release.

If you don't have [pip](https://pip.pypa.io) installed, this [Python installation guide](http://docs.python-guide.org/en/latest/starting/installation/) can guide you through the process.

## From sources

To install climdata from sources, create/activate a conda environment and then install from the repository:

```bash
# create & activate conda environment (optional)
conda create -n climdata python=3.11 -y
conda activate climdata

# install from GitHub (editable install if desired)
pip install git+https://github.com/Kaushikreddym/climdata
# or for editable development install:
# git clone https://github.com/Kaushikreddym/climdata
# cd climdata
# pip install -e .
```
