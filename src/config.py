"""Shared filesystem paths for the BART ridership dashboard.

Keeping paths here avoids hard-coded path construction throughout the app and
makes future data additions easier to locate.
"""

# region Imports
from pathlib import Path
# endregion


# region Project Paths
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RAW = DATA / "raw"
REFERENCE = DATA / "reference"
PROCESSED = DATA / "processed"
SAMPLE = DATA / "sample"
# endregion


# region Reference Inputs
STATIONS_MASTER_CSV = REFERENCE / "stations_master.csv"
STATION_ALIASES_CSV = REFERENCE / "station_aliases.csv"
# endregion


# region Raw And Sample Inputs
BART_STATIONS_DIR = RAW / "BART_Stations_2025"
BART_ROUTES_GDB = RAW / "BART_Routes" / "p20" / "shortexercise1.gdb"
RIDERSHIP_2018_DIR = RAW / "ridership_2018"
HOURLY_OD_DIR = RAW / "Hourly_OD"
HOURLY_OD_2018_CSV = HOURLY_OD_DIR / "date-hour-soo-dest-2018.csv"
HOURLY_OD_2018_SAMPLE_CSV = SAMPLE / "date-hour-soo-dest-2018-sample.csv"
# endregion


# region Processed Runtime Inputs
PROCESSED_STATIONS_GEOJSON = PROCESSED / "stations.geojson"
PROCESSED_RAW_ROUTES_GEOJSON = PROCESSED / "bart_routes.geojson"
PROCESSED_STATION_RIDERSHIP_MONTHLY_SUMMARY_CSV = PROCESSED / "station_ridership_monthly_summary.csv"
PROCESSED_STATION_RIDERSHIP_MONTHLY_SUMMARY_PARQUET = PROCESSED / "station_ridership_monthly_summary.parquet"
PROCESSED_HOURLY_STATION_MONTHLY_DIR = PROCESSED / "hourly_od_station_monthly"
PROCESSED_HOURLY_STATION_MONTHLY_PARQUET_DIR = PROCESSED / "hourly_od_station_monthly_parquet"
PROCESSED_HOURLY_STATION_MONTHLY_SUMMARY = PROCESSED / "hourly_od_station_monthly_summary.csv"
PROCESSED_HOURLY_STATION_MONTHLY_SUMMARY_PARQUET = PROCESSED / "hourly_od_station_monthly_summary.parquet"
PROCESSED_HOURLY_VALIDATION_CSV = PROCESSED / "hourly_od_total_trips_validation.csv"
PROCESSED_HOURLY_VALIDATION_PARQUET = PROCESSED / "hourly_od_total_trips_validation.parquet"
PROCESSED_HOURLY_COMPLETENESS_CSV = PROCESSED / "hourly_od_completeness_audit.csv"
PROCESSED_HOURLY_COMPLETENESS_PARQUET = PROCESSED / "hourly_od_completeness_audit.parquet"
# endregion
