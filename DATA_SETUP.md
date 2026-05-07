# Data Setup

This project expects raw BART data files under `data/raw/`. The raw folder is intentionally ignored by Git because some source files are too large for normal GitHub uploads.

## App Runtime Files

`app.py` should run from processed files, not raw source files:

```text
data/processed/stations.geojson
data/processed/bart_routes.geojson
data/processed/station_ridership_monthly_summary.csv
data/processed/station_ridership_monthly_summary.parquet
```

The hourly-derived validation outputs are also processed artifacts:

```text
data/processed/hourly_od_station_monthly_summary.csv
data/processed/hourly_od_station_monthly_summary.parquet
data/processed/hourly_od_station_monthly/YYYY-MM.csv
data/processed/hourly_od_station_monthly_parquet/YYYY-MM.parquet
data/processed/hourly_od_total_trips_validation.csv
data/processed/hourly_od_total_trips_validation.parquet
data/processed/hourly_od_completeness_audit.csv
data/processed/hourly_od_completeness_audit.parquet
```

Raw station files, route files, ridership workbooks, and hourly OD CSVs are not
required for normal app runtime when these processed files are present.

The runtime contract is:

```text
Dash app: reads data/processed/* and data/reference/*
prepare_data.py: reads data/raw/* and writes data/processed/*
notebooks: audit processed outputs and may inspect raw files during investigation
```

Where both CSV and Parquet versions exist, runtime loaders prefer Parquet and
fall back to CSV for compatibility.

## Raw Files Needed To Regenerate Processed Data

`scripts/prepare_data.py` reads:

```text
data/raw/BART_Stations_2025/doc.kml
data/raw/BART_Stations_2025/bart_geojson.json
data/raw/BART_Routes/p20/shortexercise1.gdb/
data/raw/ridership_2018/Ridership_201801.xlsx
data/raw/ridership_2018/Ridership_201802.xlsx
...
data/raw/ridership_OD_2025/Ridership_202512.xlsx
data/raw/Hourly_OD/date-hour-soo-dest-2018.csv
data/raw/Hourly_OD/date-hour-soo-dest-2019.csv
```

Hourly origin-destination CSVs are aggregated into monthly station-entry totals
by origin station. These hourly-derived summaries are the app's canonical
ridership source where available. Workbook `Total Trips` / `Total Trips OD`
values are validation benchmarks and fallback data for periods without hourly
OD coverage.

If an hourly OD file is missing a full calendar day, `scripts/prepare_data.py`
imputes that day at the hourly station level using same-month hourly station
patterns before writing monthly totals. The imputed periods and major
source-definition discrepancies are documented in:

```text
data/processed/hourly_od_completeness_audit.csv
```

```text
data/raw/Hourly_OD/date-hour-soo-dest-YYYY.csv
```

Workbook station labels are normalized before summaries are written. BART
workbooks use some legacy abbreviations that differ from the app's canonical
codes:

```text
Workbook WP -> canonical PITT -> Pittsburg/Bay Point
Workbook WD -> canonical WP   -> West Dublin/Pleasanton
Workbook WS -> canonical WD   -> Warm Springs/South Fremont
```

The `Coliseum/Airport Connector` point in the station KML is an infrastructure
marker for the airport connector, not a separate ridership station. It is
excluded from app station markers and dropdowns.

## Option 1: Use The 1000-Row Sample CSV

A 1000-row snapshot from the real 2018 CSV is included so the project has a lightweight example file in Git:

```text
data/sample/date-hour-soo-dest-2018-sample.csv
```

If you want to inspect or experiment with the sample locally, copy it into `data/raw/Hourly_OD/`:

```powershell
New-Item -ItemType Directory -Force data/raw/Hourly_OD
Copy-Item data/sample/date-hour-soo-dest-2018-sample.csv data/raw/Hourly_OD/date-hour-soo-dest-2018.csv
```

The sample is useful for local testing and repository portability. It preserves the real CSV structure and is spread across the full 2018 file, but it is still only a small snapshot of the full ridership dataset.

## Optional Processed Data Export

The Dash app no longer writes processed files during startup. To regenerate the processed exports manually, run:

```powershell
python scripts/prepare_data.py
```

That script writes files under:

```text
data/processed/
```

It writes CSV files for inspection/compatibility and Parquet files for faster
runtime loading and future time-lapse/modeling features.

To benchmark the impact of those processed summaries, run:

```powershell
python scripts/benchmark_data_loading.py --year 2025 --month 5 --runs 3
```

The benchmark compares raw hourly CSV aggregation, processed CSV loading, and
processed Parquet loading, then writes:

```text
data/processed/benchmark_data_loading.csv
```

The validation CSV is station-level. `Has Difference = True` means the station's
hourly OD monthly total and workbook monthly total are not exactly equal for
that period.

The completeness audit is period-level. `missing_full_days_imputed` means the
processed monthly ridership includes hourly imputation for missing source days.
`complete_hours_source_discrepancy` means the hourly data appears complete, but
differs materially from the workbook benchmark; the hourly OD source remains
canonical for those periods.

## Option 2: Use The Full Optional CSV

Download the full CSV from Google Drive:

```text
PASTE_GOOGLE_DRIVE_LINK_HERE
```

After downloading it, place it at exactly this path:

```text
data/raw/Hourly_OD/date-hour-soo-dest-2018.csv
```

Keep the filename unchanged if you add future analysis code that reads the full hourly origin-destination CSV.

## Why The Full Data Is Not In Git

GitHub has practical file size limits, and large raw datasets make the repository slower to clone and harder to review. This repo keeps raw source data out of Git and tracks only lightweight project files plus the sample CSV.
