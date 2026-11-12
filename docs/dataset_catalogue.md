---
title: Dataset Catalogue
description: Every climate and crop dataset climdata can fetch, with resolution, period, sources, credentials and notes for crop modellers.
---

# Dataset Catalogue

<p class="cat-lead">
Every dataset climdata can fetch, on one page: what it covers, how fine it is,
where it comes from, what you need to access it, and how well it fits a
crop-model workflow such as DSSAT, APSIM, WOFOST, SIMPLACE or LINTUL.
</p>

All of them are reached the same way: set `dataset=` and nothing else changes.
See [Datasets](datasets.md) for the call pattern and region selection.

## Which one do I need?

<div class="grid cards" markdown>

-   :material-sprout:{ .lg .middle } **Calibrating at a field trial**

    ---

    Use station records where you have them, and a fine grid where you don't.

    [`dwd`](#dwd-stations) · [`hyras`](#hyras) · [`power`](#nasa-power)

-   :material-earth:{ .lg .middle } **Gridded historical runs, any region**

    ---

    Bias-corrected global daily weather with the full crop-model variable set.

    [`mswx`](#mswx) · [`w5e5`](#w5e5) · [`era5`](#era5)

-   :material-chart-timeline-variant-shimmer:{ .lg .middle } **Climate-change projections**

    ---

    Bias-adjusted CMIP6 you can feed straight into a crop model.

    [`nexgddp`](#nex-gddp-cmip6) · [`cmip_w5e5`](#cmip6-isimip3b-cmip_w5e5)

-   :material-scale-balance:{ .lg .middle } **Benchmarking against GGCMI**

    ---

    Simulated yields from the ISIMIP crop-model ensemble, for comparison with your own runs.

    [`agri_isimip`](#isimip3b-crop-model-output)

</div>

## At a glance

| Dataset | `dataset=` | Kind | Coverage | Grid | Time step | Period | Access |
|---|---|---|---|---|---|---|---|
| [MSWX](#mswx) | `mswx` | <span class="cat-tag obs">Observation</span> | Global | 0.1° | daily | 1979 – present | <span class="cat-key">:material-key: key</span> |
| [ERA5](#era5) | `era5` | <span class="cat-tag rean">Reanalysis</span> | Global | 0.25° | hourly · daily | 1940 – present | <span class="cat-key">:material-key: key</span> |
| [W5E5](#w5e5) | `w5e5` | <span class="cat-tag obs">Observation</span> | Global | 0.5° | daily | 1979 – 2019 | <span class="cat-open">open</span> |
| [NASA POWER](#nasa-power) | `power` | <span class="cat-tag rean">Reanalysis-derived</span> | Global, point | 0.5° | daily | 1981 – present | <span class="cat-open">open</span> |
| [HYRAS](#hyras) | `hyras` | <span class="cat-tag obs">Gridded obs.</span> | Germany | 1 km | daily | 1951 – present | <span class="cat-open">open</span> |
| [HOSTRADA](#hostrada) | `hostrada` | <span class="cat-tag obs">Gridded obs.</span> | Germany | 1 km | hourly | 1995 – present | <span class="cat-open">open</span> |
| [DWD stations](#dwd-stations) | `dwd` | <span class="cat-tag obs">Station</span> | Germany | point | daily · hourly | per station | <span class="cat-open">open</span> |
| [CMIP6 raw](#cmip6-raw) | `cmip` | <span class="cat-tag proj">Projection</span> | Global | ~1° (model) | daily | 1850 – 2100 | <span class="cat-open">open</span> |
| [NEX-GDDP-CMIP6](#nex-gddp-cmip6) | `nexgddp` | <span class="cat-tag proj">Projection, BC</span> | Global | 0.25° | daily | 1950 – 2100 | <span class="cat-open">open</span> |
| [ISIMIP3b CMIP6](#cmip6-isimip3b-cmip_w5e5) | `cmip_w5e5` | <span class="cat-tag proj">Projection, BC</span> | Global | 0.5° | daily | 1850 – 2100 | <span class="cat-open">open</span> |
| [ISIMIP3b crop output](#isimip3b-crop-model-output) | `agri_isimip` | <span class="cat-tag crop">Crop-model output</span> | Global | 0.5° | annual (season) | varies, to 2100 | <span class="cat-open">open</span> |

<small>**BC** = bias-corrected. **key** = needs an account or credential file (see the card). **open** = anonymous download.</small>

## Crop-model variables by dataset

A daily crop model typically needs minimum and maximum temperature, rainfall,
solar radiation, and (for Penman–Monteith ET or VPD-driven transpiration)
humidity and wind. The table shows which datasets supply each one.

| | T<sub>min</sub> / T<sub>max</sub> | Rain | Radiation | Humidity | Wind | Also |
|---|:-:|:-:|:-:|:-:|:-:|---|
| `mswx` | :material-check: | :material-check: | :material-check: | :material-check: | :material-check: | |
| `era5` | :material-check: | :material-check: | :material-check: | :material-check: | :material-check: | `rlds` |
| `w5e5` | :material-check: | :material-check: | :material-check: | :material-check: | :material-check: | `rlds`, `ps`, `huss` |
| `power` | :material-check: | :material-check: | :material-check: | :material-minus: | :material-check: | |
| `hyras` | :material-check: | :material-check: | :material-check: | :material-check: | :material-close: | `et0`, `evpot`, soil moisture, soil temp. |
| `hostrada` | :material-close: | :material-close: | :material-check: | :material-close: | :material-check: | hourly only |
| `dwd` | :material-check: | :material-check: | :material-check: | :material-check: | :material-check: | station-dependent |
| `cmip` | :material-check: | :material-check: | :material-check: | :material-close: | :material-check: | |
| `nexgddp` | :material-check: | :material-check: | :material-check: | :material-check: | :material-check: | `rlds`, `huss` |
| `cmip_w5e5` | :material-check: | :material-check: | :material-check: | :material-check: | :material-check: | `rlds`, `ps`, `psl` |

<small>:material-minus: = available from the provider but not yet mapped in climdata.</small>

!!! crop "Units before you write a weather file"
    climdata returns `rsds` in **W m⁻²** as a daily mean. Most crop models want
    daily radiation totals in **MJ m⁻² d⁻¹**, so multiply by **0.0864**.
    Temperatures come back in °C and precipitation in mm d⁻¹.

!!! crop "CO₂ is not in any of these"
    None of the climate datasets carry atmospheric CO₂. For SSP runs, take the
    annual concentrations for your scenario from the CMIP6 forcing files (the
    ISIMIP3b protocol publishes them) and set them in the crop model yourself.
    Otherwise the CO₂-fertilisation effect is silently left out.

---

## Observations and reanalysis

### MSWX

<div class="cat-card" markdown>

<p class="cat-sub">Multi-Source Weather · GloH2O · <span class="cat-tag obs">Observation</span></p>

| | |
|---|---|
| **Spatial resolution** | 0.1° (~11 km), global land |
| **Time step** | daily (climdata); the native product is 3-hourly |
| **Period** | 1979 – present, with near-real-time updates |
| **Variables** | `tas` `tasmin` `tasmax` `pr` `rsds` `hurs` `sfcWind` |
| **Source** | Google Drive, one NetCDF per day per variable · [gloh2o.org/mswx](https://www.gloh2o.org/mswx/) |
| **Licence** | CC BY-NC 4.0, non-commercial use only |

!!! crop "For crop modellers"
    One of the best global choices for gridded historical simulations. It has
    every variable a DSSAT/APSIM weather file needs, bias-corrected against
    station observations at 0.1°, which is fine enough to tell apart
    neighbouring agro-ecological zones. The near-real-time tail also makes it
    usable for in-season yield forecasting. Precipitation merges gauge,
    satellite and reanalysis data (it is MSWEP's sibling), so rain-fed yield
    variance is captured better than with ERA5 alone.

!!! key "Credentials"
    1. Request access to the MSWX data from GloH2O (see the link above).
    2. Create a Google Cloud **service account**, download its JSON key, and
       share the MSWX Drive folder with the service-account e-mail.
    3. Pass the key with `dsinfo.MSWX.params.google_service_account=/path/to/key.json`.

    You only need the key for files that aren't already on disk. The
    [MSWX guide](MSWX_guide.md) covers every step.

**Reference:** Beck, H. E. et al. (2022). MSWX: Global 3-hourly 0.1° bias-corrected meteorological data including near-real-time updates and forecast ensembles. *Bulletin of the American Meteorological Society*, 103(3), E710–E732. [doi:10.1175/BAMS-D-21-0145.1](https://doi.org/10.1175/BAMS-D-21-0145.1)

</div>

### ERA5

<div class="cat-card" markdown>

<p class="cat-sub">ECMWF Reanalysis v5 · Copernicus C3S · <span class="cat-tag rean">Reanalysis</span></p>

| | |
|---|---|
| **Spatial resolution** | 0.25° (~31 km), global |
| **Time step** | hourly, or daily statistics |
| **Period** | 1940 – present (about 5 days behind real time) |
| **Variables** | `tas` `tasmin` `tasmax` `pr` `hurs` `sfcWind` `rsds` `rlds`, plus pressure levels |
| **Source** | Copernicus Climate Data Store via `cdsapi`, mirrored to local Zarr |

!!! crop "For crop modellers"
    Use the **hourly** fields when the process you model is sub-daily, for
    example heat stress during anthesis or hourly canopy temperature and VPD.
    ERA5 is not bias-corrected: precipitation shows the usual reanalysis
    drizzle bias (too many light-rain days, too few heavy ones), which skews
    infiltration and runoff in soil-water modules. For daily rain-fed
    simulations, prefer MSWX or W5E5, which are built on ERA5 but corrected.

!!! key "Credentials"
    Needs a free **Copernicus CDS** account. Accept the ERA5 licence on the
    dataset page, then put your API key in `~/.cdsapirc`. ERA5 is driven via
    `climdata.ERA5` rather than `ClimData(...).extract()`; see
    [Datasets › ERA5](datasets.md#per-provider-notes).

**Reference:** Hersbach, H. et al. (2020). The ERA5 global reanalysis. *Quarterly Journal of the Royal Meteorological Society*, 146(730), 1999–2049. [doi:10.1002/qj.3803](https://doi.org/10.1002/qj.3803)

</div>

### W5E5

<div class="cat-card" markdown>

<p class="cat-sub">WFDE5 over land merged with ERA5 over ocean · ISIMIP · <span class="cat-tag obs">Observation</span></p>

| | |
|---|---|
| **Spatial resolution** | 0.5° (~55 km), global |
| **Time step** | daily |
| **Period** | 1979 – 2019 |
| **Variables** | `tas` `tasmin` `tasmax` `pr` `rsds` `rlds` `hurs` `huss` `sfcWind` `ps` |
| **Source** | ISIMIP repository via `isimip-client` · [data.isimip.org](https://data.isimip.org) |

!!! crop "For crop modellers"
    This is the observational climate of **GGCMI phase 3** and the reference
    that `cmip_w5e5` was bias-adjusted to. If your study follows the
    ISIMIP3b protocol, calibrate on W5E5 and then run the `cmip_w5e5`
    projections: the historical and future climate are then consistent by
    construction, and your results compare directly with the GGCMI ensemble
    in `agri_isimip`. The 0.5° grid is too coarse for field-scale calibration,
    and the record stops in 2019.

!!! key "Credentials"
    None. Public data; install `isimip-client`.

**Reference:** Lange, S. et al. (2021). *WFDE5 over land merged with ERA5 over the ocean (W5E5 v2.0)*. ISIMIP Repository. [doi:10.48364/ISIMIP.342217](https://doi.org/10.48364/ISIMIP.342217) · Cucchi, M. et al. (2020). WFDE5: bias-adjusted ERA5 reanalysis data for impact studies. *Earth System Science Data*, 12, 2097–2120. [doi:10.5194/essd-12-2097-2020](https://doi.org/10.5194/essd-12-2097-2020)

</div>

### NASA POWER

<div class="cat-card" markdown>

<p class="cat-sub">Prediction Of Worldwide Energy Resources · Agroclimatology · <span class="cat-tag rean">Reanalysis-derived</span></p>

| | |
|---|---|
| **Spatial resolution** | point query; meteorology from MERRA-2 (0.5° × 0.625°), solar from CERES/GEWEX (1°) |
| **Time step** | daily |
| **Period** | 1981 – present |
| **Variables** | `tas` `tasmin` `tasmax` `pr` `rsds` `sfcWind` |
| **Source** | NASA POWER REST API, `AG` community · [power.larc.nasa.gov](https://power.larc.nasa.gov) |

!!! crop "For crop modellers"
    Built for agroclimatology: the `AG` community endpoint that climdata
    queries is the one DSSAT's weather tools use. It's the fastest way to get
    a complete weather file for a single site anywhere, and its
    satellite-derived radiation is often better than nothing at data-poor
    sites. Being point-only, it fits site calibration rather than regional
    runs. Precipitation is coarse reanalysis, so check it against a local
    gauge before trusting simulated water stress.

!!! key "Credentials"
    None.

**Reference:** NASA Langley Research Center, POWER Project. [power.larc.nasa.gov](https://power.larc.nasa.gov) · Sparks, A. H. (2018). nasapower: a NASA POWER global meteorology, surface solar energy and climatology data client for R. *Journal of Open Source Software*, 3(30), 1035. [doi:10.21105/joss.01035](https://doi.org/10.21105/joss.01035)

</div>

### HYRAS

<div class="cat-card" markdown>

<p class="cat-sub">Hydrologische Rasterdaten · Deutscher Wetterdienst · <span class="cat-tag obs">Gridded observation</span></p>

| | |
|---|---|
| **Spatial resolution** | 1 km, Germany (and the transboundary river basins for some variables) |
| **Time step** | daily |
| **Period** | 1951 – present (start year varies by variable) |
| **Variables** | `tas` `tasmin` `tasmax` `pr` `rsds` `hurs`, plus DWD grids for `et0` (FAO), `evpot`, `soilMoist`, `soilTemp` (5 cm) |
| **Source** | DWD open-data server · [opendata.dwd.de › grids_germany/daily](https://opendata.dwd.de/climate_environment/CDC/grids_germany/daily/) |

!!! crop "For crop modellers"
    The standard choice for crop modelling in Germany. At 1 km a grid cell is
    close to the scale of a field trial, so it can stand in for a missing
    station. The extra DWD grids are unusually useful: **FAO-56 reference
    ET** (`et0`) can go straight into models that accept ET₀, and **5 cm soil
    temperature** helps simulate sowing-to-emergence. HYRAS has **no wind**, so
    for Penman–Monteith inside the model either use `et0` directly or add
    wind from `dwd` stations or `hostrada`.

!!! key "Credentials"
    None. Anonymous HTTPS. Files are national annual NetCDFs, so the first
    download of a year is large even for a single point.

**Reference:** Razafimaharo, C., Krähenmann, S., Höpp, S., Rauthe, M. & Deutschländer, T. (2020). New high-resolution gridded dataset of daily mean, minimum, and maximum temperature and relative humidity for Central Europe (HYRAS). *Theoretical and Applied Climatology*, 142, 1531–1553. [doi:10.1007/s00704-020-03388-w](https://doi.org/10.1007/s00704-020-03388-w) · Rauthe, M. et al. (2013). A Central European precipitation climatology – Part I. *Meteorologische Zeitschrift*, 22(3), 235–256. [doi:10.1127/0941-2948/2013/0436](https://doi.org/10.1127/0941-2948/2013/0436)

</div>

### HOSTRADA

<div class="cat-card" markdown>

<p class="cat-sub">Hourly Station-based and Raster-based Terrestrial Atmospheric Dataset · DWD · <span class="cat-tag obs">Gridded observation</span></p>

| | |
|---|---|
| **Spatial resolution** | 1 km, Germany (regular lat/lon, EPSG:4326) |
| **Time step** | hourly |
| **Period** | 1995 – present |
| **Variables** | `rsds`, `wspd` (wind speed) |
| **Source** | DWD open-data server, monthly NetCDF · [opendata.dwd.de › grids_germany/hourly/hostrada](https://opendata.dwd.de/climate_environment/CDC/grids_germany/hourly/hostrada/) |

!!! crop "For crop modellers"
    A complement to HYRAS, not a replacement: climdata currently maps only
    radiation and wind. Use it to fill HYRAS's missing wind, or to drive
    hourly canopy photosynthesis and light-interception models that need the
    diurnal radiation curve rather than a daily total. Shapefile requests are
    clipped to the polygon's bounding box, not the polygon itself.

!!! key "Credentials"
    None.

**Reference:** Deutscher Wetterdienst, Climate Data Center. *HOSTRADA v1.0: hourly gridded surface climate data for Germany.* [opendata.dwd.de](https://opendata.dwd.de/climate_environment/CDC/grids_germany/hourly/hostrada/)

</div>

### DWD stations

<div class="cat-card" markdown>

<p class="cat-sub">Deutscher Wetterdienst station network · <span class="cat-tag obs">Station observation</span></p>

| | |
|---|---|
| **Spatial resolution** | point; climdata picks the nearest station(s) |
| **Time step** | daily, hourly |
| **Period** | varies by station; many start before 1950 |
| **Variables** | `tas` `tasmin` `tasmax` `pr` `rsds` `hurs` `sfcWind` |
| **Source** | DWD Climate Data Center via the `wetterdienst` library · [opendata.dwd.de/climate_environment/CDC](https://opendata.dwd.de/climate_environment/CDC/) |

!!! crop "For crop modellers"
    The ground truth for calibrating and validating at German experimental
    sites. Not every station measures every variable, radiation in
    particular, and gaps are common, so run the series through
    `climdata.impute` before writing a weather file. A gap of a single day in
    T<sub>max</sub> stops most crop models. The
    [DWD export format](dwd_format_guide.md) page describes the station-file
    layout.

!!! key "Credentials"
    None. Install `wetterdienst`.

**Reference:** Deutscher Wetterdienst, Climate Data Center. [opendata.dwd.de/climate_environment/CDC](https://opendata.dwd.de/climate_environment/CDC/) · Access library: [wetterdienst](https://github.com/earthobservations/wetterdienst)

</div>

---

## Climate projections

!!! crop "Don't feed raw GCM output to a crop model"
    A raw GCM with a 2 °C warm bias in T<sub>max</sub> moves flowering by
    weeks and wipes out the yield signal you are trying to measure. Use one of
    the bias-corrected products below, or correct `cmip` yourself with
    [BCSD](bcsd_guide.md). Run several GCMs (ideally five or more) and report
    the spread: inter-model spread usually dominates crop-model uncertainty
    by mid-century.

### CMIP6 raw

<div class="cat-card" markdown>

<p class="cat-sub">CMIP6 via the Pangeo cloud catalogue · <span class="cat-tag proj">Projection, raw</span></p>

| | |
|---|---|
| **Spatial resolution** | each model's native grid, typically ~1° (~100 km) |
| **Time step** | daily (`table_id=day`) |
| **Period** | historical 1850 – 2014; SSPs 2015 – 2100 |
| **Variables** | `tas` `tasmin` `tasmax` `pr` `rsds` `sfcWind` |
| **Source** | Google Cloud Storage via `intake-esm` (Pangeo CMIP6 catalogue) |

!!! crop "For crop modellers"
    Use it when you need a model, member or scenario the bias-corrected
    products don't carry, or when you want to control the bias correction
    yourself, e.g. correcting to MSWX or HYRAS instead of a global reference.
    Plan on running it through `climdata.sdba` before it goes near a crop
    model.

!!! key "Credentials"
    None. The data is public; install `intake-esm`.

**Reference:** Eyring, V. et al. (2016). Overview of the Coupled Model Intercomparison Project Phase 6 (CMIP6) experimental design and organization. *Geoscientific Model Development*, 9, 1937–1958. [doi:10.5194/gmd-9-1937-2016](https://doi.org/10.5194/gmd-9-1937-2016)

</div>

### NEX-GDDP-CMIP6

<div class="cat-card" markdown>

<p class="cat-sub">NASA Earth Exchange Global Daily Downscaled Projections · <span class="cat-tag proj">Projection, bias-corrected</span></p>

| | |
|---|---|
| **Spatial resolution** | 0.25° (~25 km), global |
| **Time step** | daily |
| **Period** | historical 1950 – 2014; `ssp126` `ssp245` `ssp370` `ssp585` 2015 – 2100 |
| **Variables** | `tas` `tasmin` `tasmax` `pr` `rsds` `rlds` `hurs` `huss` `sfcWind` |
| **Source** | NASA NCCS THREDDS server (HTTP subsetting) · [nccs.nasa.gov › NEX-GDDP-CMIP6](https://www.nccs.nasa.gov/services/data-collections/land-based-products/nex-gddp-cmip6) |

!!! crop "For crop modellers"
    The finest-resolution ready-to-use projections here, with about 35 GCMs
    and the full crop-model variable set. It's a good default for regional
    yield-impact studies. BCSD corrects each model's distribution against
    GMFD observations, but day-to-day sequencing comes from the GCM, so
    dry-spell lengths and heat-wave persistence are only as realistic as the
    parent model. Check them if drought timing drives your results. Not every
    model provides every variable.

!!! key "Credentials"
    None.

**Reference:** Thrasher, B. et al. (2022). NASA Global Daily Downscaled Projections, CMIP6. *Scientific Data*, 9, 262. [doi:10.1038/s41597-022-01393-4](https://doi.org/10.1038/s41597-022-01393-4)

</div>

### CMIP6 ISIMIP3b (`cmip_w5e5`)

<div class="cat-card" markdown>

<p class="cat-sub">CMIP6 bias-adjusted and downscaled to W5E5 · ISIMIP3BASD · <span class="cat-tag proj">Projection, bias-corrected</span></p>

| | |
|---|---|
| **Spatial resolution** | 0.5°, global (the W5E5 grid) |
| **Time step** | daily |
| **Period** | historical 1850 – 2014; `ssp126` `ssp370` `ssp585` 2015 – 2100 |
| **Variables** | `tas` `tasmin` `tasmax` `pr` `rsds` `rlds` `hurs` `sfcWind` `ps` `psl` |
| **Models** | primary set: `gfdl-esm4` `ipsl-cm6a-lr` `mpi-esm1-2-hr` `mri-esm2-0` `ukesm1-0-ll` (lowercase) |
| **Source** | ISIMIP repository via `isimip-client` · [data.isimip.org](https://data.isimip.org) |

!!! crop "For crop modellers"
    The climate forcing behind GGCMI phase 3. Pick it when your results should
    sit alongside the published global crop-model ensemble, or when you need
    the trend-preserving ISIMIP3BASD method, which keeps each model's
    warming signal per quantile. That matters for heat-stress thresholds in
    T<sub>max</sub>. Five models are few for a full uncertainty range, and the
    0.5° grid smooths out topography in mountainous cropping regions.

!!! key "Credentials"
    None.

**Reference:** Lange, S. (2019). Trend-preserving bias adjustment and statistical downscaling with ISIMIP3BASD (v1.0). *Geoscientific Model Development*, 12, 3055–3070. [doi:10.5194/gmd-12-3055-2019](https://doi.org/10.5194/gmd-12-3055-2019)

</div>

---

## Crop-model output

### ISIMIP3b crop-model output

<div class="cat-card" markdown>

<p class="cat-sub">GGCMI phase 3 via ISIMIP · `agri_isimip` · <span class="cat-tag crop">Simulated</span></p>

| | |
|---|---|
| **Spatial resolution** | 0.5°, global |
| **Time step** | annual, one value per growing season |
| **Period** | historical and SSP experiments, up to 2100; varies by model and forcing |
| **Variables** | `yield` (t ha⁻¹) `biom` `lai` `evap` `aet` `plantyear` |
| **Dimensions** | impact model × climate forcing × period × crop × irrigation (`firr` / `noirr`) |
| **Source** | ISIMIP repository, `ISIMIP3b/OutputData/agriculture` · [data.isimip.org](https://data.isimip.org) |

!!! crop "For crop modellers"
    These aren't inputs: they're other models' answers. Use them to benchmark
    your own model against the GGCMI ensemble (e.g. LPJmL, pDSSAT, EPIC-IIASA,
    SIMPLACE-LINTUL5) driven by the same `cmip_w5e5` climate, or to put
    your site results in a global context. Compare fully-irrigated and
    rain-fed runs to separate heat effects from water stress. Which
    model/crop/forcing combinations exist can't be guessed, so call
    `AGRI_ISIMIP(cfg).explore()` first.

!!! key "Credentials"
    None.

**Reference:** Jägermeyr, J. et al. (2021). Climate impacts on global agriculture emerge earlier in new generation of climate and crop models. *Nature Food*, 2, 873–885. [doi:10.1038/s43016-021-00400-y](https://doi.org/10.1038/s43016-021-00400-y)

</div>

---

<small>The resolutions, periods and variable lists on this page come from
`climdata/conf/mappings/parameters.yaml`. `cd.explore("<NAME>")` prints the same
metadata at runtime. Provider periods and model availability change over
time, so check the source before a large download.</small>
