# BART Ridership Dashboard

![App Screenshot](assets/screenshot.PNG)

Interactive Dash + Plotly app with GeoPandas/Shapely maps to explore BART ridership and routes.

## Features
- Interactive station/route maps
- Station-level ridership comparison for 2018
- Filters for route and station pair
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
- `src/data_loader.py`: runtime data loading and ridership shaping.
- `src/route_builder.py`: curated BART route geometry corrections.
- `src/station_mapping.py`: station abbreviation to station-name mapping.
- `scripts/prepare_data.py`: optional raw-to-processed data export script.

## Data
Raw source files belong under `data/raw/`, which is ignored by Git because some files are too large for normal GitHub uploads. See [DATA_SETUP.md](DATA_SETUP.md) for the required files, the included sample CSV, and the full-data download link.

## Tech Stack
Python, Dash, Plotly, GeoPandas, Shapely, Pandas.


## Notes
January 2018 workbook lacks a monthly Total Trips sheet, so the pipeline uses the higher-granularity 2018 hourly OD CSV as a source-specific fallback. 
It aggregates hourly OD rows by Origin to match the station-entry metric used for all other monthly workbooks.
