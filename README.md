# BART Ridership Dashboard

![App Screenshot](assets/screenshot.PNG)

Interactive Dash + Plotly app with GeoPandas/Shapely maps to explore BART ridership and routes.

## Features
- Interactive station/route maps
- Station-level monthly ridership comparison
- Filters for route service year, ridership period, route, and station pair
- Geospatial visualizations with GeoPandas and Shapely

## Quickstart (Windows)

```powershell
python -m venv .venv
.\.venv\Scripts\Activate
pip install -r requirements.txt
python app.py
```

App runs at http://127.0.0.1:8050/

## Project Layout

- `app.py`: Dash app entrypoint, layout, and callbacks.
- `src/data_loader.py`: processed artifact loading plus preprocessing helpers.
- `src/route_builder.py`: curated BART route geometry corrections.
- `src/station_mapping.py`: station abbreviation to station-name mapping.
- `scripts/prepare_data.py`: optional raw-to-processed data export script.

## Data
Raw source files belong under `data/raw/`, which is ignored by Git because some files are too large for normal GitHub uploads. The app runs from tracked processed files in `data/processed/`; raw files are only needed when regenerating those artifacts. See [DATA_SETUP.md](DATA_SETUP.md) for the required files, the included sample CSV, and the full-data download link.

Processed ridership tables are exported as both CSV and Parquet. Runtime loaders prefer Parquet and fall back to CSV for compatibility.

## Data Quality
Hourly OD files are the canonical ridership source where available. Workbook Total Trips files are validation benchmarks and fallback inputs only. Missing full hourly OD days are imputed during preprocessing, and imputed periods plus major source-definition discrepancies are documented in `data/processed/hourly_od_completeness_audit.csv`.

## Benchmarking
Run `python scripts/benchmark_data_loading.py` to compare raw hourly CSV aggregation against processed CSV and processed Parquet summaries. Results are written to `data/processed/benchmark_data_loading.csv` for resume/portfolio metrics.

## Tech Stack
Python, Dash, Plotly, GeoPandas, Shapely, Pandas.


## Notes
Workbook station abbreviations are normalized into canonical app codes during preprocessing, including legacy labels around Pittsburg/Bay Point, West Dublin/Pleasanton, and Warm Springs/South Fremont.
