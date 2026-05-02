# Data Setup

This project expects raw BART data files under `data/raw/`. The raw folder is intentionally ignored by Git because some source files are too large for normal GitHub uploads.

## Required Raw Files

`app.py` currently reads these files at startup:

```text
data/raw/BART_Stations_2025/doc.kml
data/raw/BART_Routes/p20/shortexercise1.gdb/
data/raw/ridership_2018/Ridership_201801.xlsx
```

The optional data preparation script also reads:

```text
data/raw/BART_Stations_2025/bart_geojson.json
```

The large hourly origin-destination CSV is not currently required to run the app:

```text
data/raw/date-hour-soo-dest-2018.csv
```

## Option 1: Use The 1000-Row Sample CSV

A 1000-row snapshot from the real 2018 CSV is included so the project has a lightweight example file in Git:

```text
data/sample/date-hour-soo-dest-2018-sample.csv
```

If you want to inspect or experiment with the sample locally, copy it into `data/raw/`:

```powershell
Copy-Item data/sample/date-hour-soo-dest-2018-sample.csv data/raw/date-hour-soo-dest-2018.csv
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

## Option 2: Use The Full Optional CSV

Download the full CSV from Google Drive:

```text
PASTE_GOOGLE_DRIVE_LINK_HERE
```

After downloading it, place it at exactly this path:

```text
data/raw/date-hour-soo-dest-2018.csv
```

Keep the filename unchanged if you add future analysis code that reads the full hourly origin-destination CSV.

## Why The Full Data Is Not In Git

GitHub has practical file size limits, and large raw datasets make the repository slower to clone and harder to review. This repo keeps raw source data out of Git and tracks only lightweight project files plus the sample CSV.
