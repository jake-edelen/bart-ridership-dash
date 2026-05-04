"""Runtime data loading and shaping for the BART dashboard.

This module prepares the in-memory tables consumed by `app.py`. It deliberately
does not write processed files; disk exports belong in `scripts/prepare_data.py`.
The optional hourly origin-destination loader is kept here as a future hook for
hourly filters, multi-year analysis, and forecasting work.
"""

# region Imports
from dataclasses import dataclass
import re
from pathlib import Path

import geopandas as gpd
import pandas as pd

from .config import (
    BART_ROUTES_GDB,
    BART_STATIONS_DIR,
    HOURLY_OD_2018_CSV,
    HOURLY_OD_2018_SAMPLE_CSV,
    RAW,
    RIDERSHIP_2018_DIR,
    PROCESSED,
)
from .route_builder import build_current_routes
from .station_mapping import STATION_MAPPING, normalize_station_code
# endregion


# region Column Definitions
STATION_METADATA_COLUMNS = [
    "description",
    "timestamp",
    "begin",
    "end",
    "altitudeMode",
    "tessellate",
    "extrude",
    "visibility",
    "drawOrder",
    "icon",
    "snippet",
]

HOURLY_OD_COLUMNS = [
    "Date",
    "Hour",
    "ORIGIN",
    "DESTINATION",
    "Number of Exits",
]
# endregion


# region Data Containers
@dataclass
class AppData:
    """Container for all data objects needed by the Dash app.

    Attributes:
        stations_gdf: Station point GeoDataFrame with lat/lon and ridership.
        raw_routes_gdf: Infrastructure line GeoDataFrame loaded from the GDB.
        route_to_gdf: Curated route GeoDataFrame with route names and colors.
        station_ridership: Aggregated ridership by entry station abbreviation.
        station_mapping: Abbreviation-to-full-station-name lookup.
    """

    stations_gdf: gpd.GeoDataFrame
    raw_routes_gdf: gpd.GeoDataFrame
    route_to_gdf: gpd.GeoDataFrame
    station_ridership: pd.DataFrame
    station_mapping: dict
# endregion


# region App Data Loader
def load_app_data():
    """Load and assemble the complete dataset required by the Dash app.

    Returns:
        AppData containing station geometry, curated routes, station ridership,
        and station name mappings.
    """
    stations_gdf = load_stations()
    raw_routes_gdf = load_raw_routes()
    route_to_gdf = build_current_routes(raw_routes_gdf, stations_gdf)

    df_ridership_long = load_monthly_ridership_long(RIDERSHIP_2018_DIR / "Ridership_201801.xlsx")
    station_ridership = summarize_station_ridership(df_ridership_long)
    stations_gdf = attach_station_ridership(stations_gdf, station_ridership)

    return AppData(
        stations_gdf=stations_gdf,
        raw_routes_gdf=raw_routes_gdf,
        route_to_gdf=route_to_gdf,
        station_ridership=station_ridership,
        station_mapping=STATION_MAPPING,
    )
# endregion


# region Station And Route Loading
def load_stations():
    """Load BART station points and add latitude/longitude helper columns.

    Returns:
        GeoDataFrame of stations with metadata columns removed where present.
    """
    stations_gdf = gpd.read_file(BART_STATIONS_DIR / "doc.kml", driver="KML")
    stations_gdf = stations_gdf.drop(columns=STATION_METADATA_COLUMNS, errors="ignore")
    stations_gdf["lat"] = stations_gdf["geometry"].y
    stations_gdf["lon"] = stations_gdf["geometry"].x
    return stations_gdf


def load_raw_routes():
    """Load the raw route geometries from the file geodatabase.

    The source layer contains geometry and length, but not stable route names or
    colors; semantic route metadata is added in `route_builder.py`.
    """
    return gpd.read_file(BART_ROUTES_GDB, layer="BARTLine")
# endregion


# region Ridership Loading
def load_monthly_ridership_long(file_path):
    """Load a monthly station matrix workbook and convert it to long format.

    Args:
        file_path: Path to one BART monthly ridership Excel workbook.

    Returns:
        DataFrame with `Exit Station`, `Entry Station`, and `Ridership` columns.
    """
    excel_file = pd.ExcelFile(file_path)
    try:
        target_sheet = _select_total_trips_sheet(excel_file.sheet_names)
    finally:
        excel_file.close()

    df_ridership = pd.read_excel(file_path, sheet_name=target_sheet, header=None)
    header_row = _find_station_header_row(df_ridership)
    if header_row is None:
        raise ValueError(f"Could not find station header row in {file_path} [{target_sheet}].")

    raw_entry_stations = [
        _normalize_station_label(name)
        for name in df_ridership.iloc[header_row, 1:].tolist()
    ]
    entry_columns = [
        (column_index, station_code)
        for column_index, station_code in enumerate(raw_entry_stations, start=1)
        if station_code in STATION_MAPPING
    ]
    exit_station_labels = [
        _normalize_station_label(name)
        for name in df_ridership.iloc[header_row + 1:, 0].tolist()
    ]
    valid_exit_rows = [
        station_code in STATION_MAPPING
        for station_code in exit_station_labels
    ]

    entry_stations = [station_code for _, station_code in entry_columns]
    exit_stations = [
        station_code
        for station_code, is_valid_station in zip(exit_station_labels, valid_exit_rows)
        if is_valid_station
    ]
    entry_stations = _make_unique(entry_stations)

    df_ridership_clean = df_ridership.iloc[
        header_row + 1:,
        [column_index for column_index, _ in entry_columns],
    ].copy()
    df_ridership_clean = df_ridership_clean.loc[valid_exit_rows].copy()
    df_ridership_clean.columns = entry_stations
    df_ridership_clean.insert(0, "Exit Station", exit_stations)
    df_ridership_clean.iloc[:, 1:] = df_ridership_clean.iloc[:, 1:].apply(
        pd.to_numeric,
        errors="coerce",
    )

    df_ridership_long = df_ridership_clean.melt(
        id_vars=["Exit Station"],
        var_name="Entry Station",
        value_name="Ridership",
    )
    df_ridership_long["Entry Station"] = pd.Categorical(
        df_ridership_long["Entry Station"],
        categories=entry_stations,
        ordered=True,
    )
    df_ridership_long["Exit Station"] = pd.Categorical(
        df_ridership_long["Exit Station"],
        categories=exit_stations,
        ordered=True,
    )
    return df_ridership_long


def summarize_station_ridership(df_ridership_long):
    """Aggregate long-format ridership by entry station.

    Args:
        df_ridership_long: Output from `load_monthly_ridership_long`.

    Returns:
        DataFrame keyed by station abbreviation with total ridership and mapped
        full station name.
    """
    station_ridership = (
        df_ridership_long.groupby("Entry Station", observed=False)["Ridership"]
        .sum()
        .reset_index()
    )
    station_ridership = station_ridership[
        station_ridership["Entry Station"].isin(STATION_MAPPING)
    ]
    station_ridership["Full Station Name"] = station_ridership["Entry Station"].map(STATION_MAPPING)
    return station_ridership


def build_monthly_station_ridership_summaries(raw_root=RAW):
    """Aggregate every discovered monthly workbook to station-level summaries.

    Args:
        raw_root: Directory containing `ridership_*` raw workbook folders.

    Returns:
        DataFrame with one row per year/month/entry-station combination.
    """
    monthly_summaries = []

    for workbook in discover_ridership_workbooks(raw_root):
        year, month = parse_ridership_period(workbook)
        station_summary = summarize_station_ridership(load_monthly_ridership_long(workbook))
        station_summary.insert(0, "Month", month)
        station_summary.insert(0, "Year", year)
        station_summary["Period"] = f"{year}-{month:02d}"
        station_summary["Source File"] = str(workbook)
        monthly_summaries.append(station_summary)

    if not monthly_summaries:
        return pd.DataFrame(
            columns=[
                "Year",
                "Month",
                "Entry Station",
                "Ridership",
                "Full Station Name",
                "Period",
                "Source File",
            ]
        )

    return pd.concat(monthly_summaries, ignore_index=True)


def write_monthly_station_ridership_summaries(output_path=None, raw_root=RAW):
    """Write processed monthly station ridership summaries to CSV.

    Args:
        output_path: Optional destination. Defaults to
            `data/processed/station_ridership_monthly_summary.csv`.
        raw_root: Directory containing raw ridership workbooks.

    Returns:
        Path to the written summary CSV.
    """
    output_path = Path(output_path) if output_path else monthly_station_summary_path()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    build_monthly_station_ridership_summaries(raw_root).to_csv(output_path, index=False)
    return output_path


def monthly_station_summary_path():
    """Return the standard processed monthly ridership summary path."""
    return PROCESSED / "station_ridership_monthly_summary.csv"


def load_station_ridership_for_period(year, month, processed_path=None, raw_root=RAW):
    """Load station ridership for one selected year/month.

    The app prefers the processed summary CSV when it exists, then falls back to
    parsing the matching raw workbook. This keeps startup light while allowing
    future production runs to precompute all monthly summaries.
    """
    year = int(year)
    month = int(month)
    processed_path = Path(processed_path) if processed_path else monthly_station_summary_path()

    if processed_path.exists():
        monthly_summary = pd.read_csv(processed_path)
        selected_summary = monthly_summary[
            (monthly_summary["Year"].astype(int) == year)
            & (monthly_summary["Month"].astype(int) == month)
        ].copy()
        if not selected_summary.empty:
            return _normalize_station_summary_columns(selected_summary)

    workbook = find_ridership_workbook(year, month, raw_root)
    return summarize_station_ridership(load_monthly_ridership_long(workbook))


def attach_station_ridership(stations_gdf, station_ridership):
    """Merge aggregated ridership onto station geometries.

    Unmatched stations are currently filled with zero so the map can render, but
    those rows should be audited separately before interpreting zero as true
    zero ridership.
    """
    stations_with_ridership = stations_gdf.merge(
        station_ridership,
        left_on="Name",
        right_on="Full Station Name",
        how="left",
    )
    stations_with_ridership["Ridership"] = pd.to_numeric(
        stations_with_ridership["Ridership"],
        errors="coerce",
    ).fillna(0)
    return stations_with_ridership
# endregion


# region Optional Hourly Origin-Destination Loading
def load_hourly_origin_destination(path=None):
    """Load optional hourly origin-destination data for future features.

    Args:
        path: Optional explicit CSV path. When omitted, the full raw CSV is used
            if present; otherwise the tracked 1000-row sample is loaded.

    Returns:
        DataFrame with `Date`, `Hour`, `ORIGIN`, `DESTINATION`, and
        `Number of Exits` columns.
    """
    csv_path = Path(path) if path else HOURLY_OD_2018_CSV
    if not csv_path.exists():
        csv_path = HOURLY_OD_2018_SAMPLE_CSV

    return pd.read_csv(csv_path, header=None, names=HOURLY_OD_COLUMNS)
# endregion


# region Ridership Workbook Discovery
def discover_ridership_workbooks(raw_root=RAW):
    """Return all monthly ridership workbooks discovered under `data/raw`."""
    raw_root = Path(raw_root)
    workbooks = set(raw_root.glob("ridership*/Ridership_*.xlsx"))
    workbooks.update(raw_root.glob("Ridership_*.xlsx"))
    return tuple(sorted(workbooks))


def get_available_ridership_periods(raw_root=RAW):
    """Return discovered ridership periods as `(year, month)` tuples."""
    return tuple(
        sorted(
            {
                parse_ridership_period(workbook)
                for workbook in discover_ridership_workbooks(raw_root)
            }
        )
    )


def find_ridership_workbook(year, month, raw_root=RAW):
    """Find the raw workbook matching one ridership year/month."""
    requested_period = (int(year), int(month))

    for workbook in discover_ridership_workbooks(raw_root):
        if parse_ridership_period(workbook) == requested_period:
            return workbook

    raise FileNotFoundError(f"No ridership workbook found for {year}-{int(month):02d}.")


def parse_ridership_period(path):
    """Extract `(year, month)` from a `Ridership_YYYYMM.xlsx` file path."""
    match = re.search(r"Ridership_(\d{4})(\d{2})\.xlsx$", Path(path).name)
    if not match:
        raise ValueError(f"Could not infer ridership period from {path}.")
    return int(match.group(1)), int(match.group(2))


def extract_ridership_station_codes(file_path, sheet_name=None):
    """Extract canonical station codes from old and new BART workbook layouts.

    Args:
        file_path: Path to a monthly ridership workbook.
        sheet_name: Optional sheet to inspect. When omitted, the total-trips
            sheet is preferred regardless of sheet order.

    Returns:
        Tuple of canonical station codes found in the OD matrix header row.
    """
    file_path = Path(file_path)
    excel_file = pd.ExcelFile(file_path)
    try:
        target_sheet = sheet_name or _select_total_trips_sheet(excel_file.sheet_names)
    finally:
        excel_file.close()

    df = pd.read_excel(file_path, sheet_name=target_sheet, header=None)
    header_row = _find_station_header_row(df)

    if header_row is None:
        raise ValueError(f"Could not find station header row in {file_path} [{target_sheet}].")

    station_codes = []
    for value in df.iloc[header_row, 1:].dropna().tolist():
        code = _normalize_station_label(value)
        if code != "Exits":
            station_codes.append(code)
    return tuple(station_codes)


def get_ridership_station_codes_by_year(raw_root=RAW):
    """Return a year-to-station-code-set mapping for discovered workbooks."""
    codes_by_year = {}

    for workbook in discover_ridership_workbooks(raw_root):
        year = _year_from_ridership_path(workbook)
        codes_by_year.setdefault(year, set()).update(extract_ridership_station_codes(workbook))

    return {
        year: tuple(sorted(codes))
        for year, codes in sorted(codes_by_year.items())
    }


def _select_total_trips_sheet(sheet_names):
    """Select the total-trips sheet, tolerating 2025 files with shuffled order."""
    for preferred_name in ("Total Trips", "Total Trips OD"):
        if preferred_name in sheet_names:
            return preferred_name
    return sheet_names[0]


def _find_station_header_row(df):
    """Find the OD matrix station-code header row in old or new workbooks."""
    for row_index in range(min(12, len(df))):
        row_values = {
            _normalize_station_label(value)
            for value in df.iloc[row_index, 1:].dropna().tolist()
        }
        if {"RM", "EN", "EP"}.issubset(row_values):
            return row_index
    return None


def _year_from_ridership_path(path):
    """Extract a four-digit year from a ridership workbook path."""
    for part in Path(path).parts:
        if part.startswith("ridership_OD_"):
            return part.removeprefix("ridership_OD_")
        if part.startswith("ridership_"):
            return part.removeprefix("ridership_")
    raise ValueError(f"Could not infer ridership year from {path}.")
# endregion


# region Helpers
def _normalize_station_summary_columns(station_summary):
    """Return period-filtered summary rows with runtime-compatible columns."""
    station_summary = station_summary.copy()
    station_summary["Entry Station"] = station_summary["Entry Station"].astype(str)
    station_summary["Ridership"] = pd.to_numeric(
        station_summary["Ridership"],
        errors="coerce",
    ).fillna(0)

    if "Full Station Name" not in station_summary.columns:
        station_summary["Full Station Name"] = station_summary["Entry Station"].map(STATION_MAPPING)

    return station_summary[["Entry Station", "Ridership", "Full Station Name"]]


def _normalize_station_label(name):
    """Normalize station labels from the Excel matrix header cells."""
    if isinstance(name, (int, float)) and not pd.isna(name):
        name = str(int(name))
    return normalize_station_code(str(name))


def _make_unique(names):
    """Return names with `_dup` suffixes added where duplicates occur."""
    seen = set()
    unique_names = []

    for name in names:
        unique_name = name
        while unique_name in seen:
            unique_name += "_dup"
        seen.add(unique_name)
        unique_names.append(unique_name)

    return unique_names
# endregion
