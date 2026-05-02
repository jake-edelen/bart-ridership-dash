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
PROCESSED = DATA / "processed"
SAMPLE = DATA / "sample"
# endregion


# region Raw And Sample Inputs
BART_STATIONS_DIR = RAW / "BART_Stations_2025"
BART_ROUTES_GDB = RAW / "BART_Routes" / "p20" / "shortexercise1.gdb"
RIDERSHIP_2018_DIR = RAW / "ridership_2018"
HOURLY_OD_2018_CSV = RAW / "date-hour-soo-dest-2018.csv"
HOURLY_OD_2018_SAMPLE_CSV = SAMPLE / "date-hour-soo-dest-2018-sample.csv"
# endregion
