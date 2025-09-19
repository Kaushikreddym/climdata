# Common Concepts in climdata

This page describes common terminology, configuration patterns, and reusable components in the `climdata` package.

## Configuration Files

- All configuration is managed via Hydra and YAML files in the `conf/` directory.
- See [`config.yaml`](../climdata/conf/config.yaml) for the main entry point.

## Standard Variable Names

- Variables follow CF conventions (see [`variables.yaml`](../climdata/conf/mappings/variables.yaml)).
- Example: `tas` for air temperature, `pr` for precipitation.

## Output Schema

All outputs are standardized to the following columns:

| Column      | Description                      |
|-------------|----------------------------------|
| latitude    | Latitude of observation/grid     |
| longitude   | Longitude of observation/grid    |
| time        | Timestamp                        |
| source      | Data source/provider             |
| variable    | Variable name                    |
| value       | Observed or modeled value        |
| units       | Units of measurement             |

## Regions and Bounds

- Regions are defined in `config.yaml` under `bounds`.
- Example: `europe`, `global`.

## Usage Patterns

- Use `climdata.load_config()` to load configuration.
- Use `climdata.DWD(cfg)` or `climdata.MSWX(cfg)` for dataset access.

---

Add more shared concepts as your documentation grows.