"""Runtime data loading and shaping for the BART dashboard.

This module prepares the in-memory tables consumed by `app.py` and exposes
writer helpers used by `scripts/prepare_data.py`. The Dash app should call the
loaders only; processed artifact generation remains a separate script step.
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
    HOURLY_OD_DIR,
    HOURLY_OD_2018_CSV,
    HOURLY_OD_2018_SAMPLE_CSV,
    PROCESSED_HOURLY_STATION_MONTHLY_DIR,
    PROCESSED_HOURLY_STATION_MONTHLY_SUMMARY,
    PROCESSED_HOURLY_VALIDATION_CSV,
    PROCESSED_RAW_ROUTES_GEOJSON,
    PROCESSED_STATIONS_GEOJSON,
    RAW,
    RIDERSHIP_2018_DIR,
    PROCESSED,
)
from .route_builder import build_current_routes
from .station_mapping import (
    IGNORED_STATION_NAMES,
    STATION_MAPPING,
    normalize_hourly_od_station_code,
    normalize_workbook_station_code,
)
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

    station_ridership = load_station_ridership_for_period(2018, 1)
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
    """Load processed BART station points and add latitude/longitude columns.

    Returns:
        GeoDataFrame of stations with metadata columns removed where present.
    """
    if not PROCESSED_STATIONS_GEOJSON.exists():
        raise FileNotFoundError(
            f"Processed stations file not found: {PROCESSED_STATIONS_GEOJSON}. "
            "Run `python scripts/prepare_data.py` to regenerate it."
        )

    stations_gdf = gpd.read_file(PROCESSED_STATIONS_GEOJSON)
    stations_gdf = stations_gdf.drop(columns=STATION_METADATA_COLUMNS, errors="ignore")
    stations_gdf = stations_gdf[~stations_gdf["Name"].isin(IGNORED_STATION_NAMES)].copy()
    stations_gdf["lat"] = stations_gdf["geometry"].y
    stations_gdf["lon"] = stations_gdf["geometry"].x
    return stations_gdf


def load_raw_routes():
    """Load processed infrastructure route geometries.

    The processed source contains geometry and length, but not stable route
    names or colors; semantic route metadata is added in `route_builder.py`.
    """
    if not PROCESSED_RAW_ROUTES_GEOJSON.exists():
        raise FileNotFoundError(
            f"Processed route file not found: {PROCESSED_RAW_ROUTES_GEOJSON}. "
            "Run `python scripts/prepare_data.py` to regenerate it."
        )

    return gpd.read_file(PROCESSED_RAW_ROUTES_GEOJSON)
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


def hourly_od_path_for_year(year, raw_root=RAW):
    """Return the expected hourly OD CSV path for one year."""
    if Path(raw_root) == RAW:
        return HOURLY_OD_DIR / f"date-hour-soo-dest-{int(year)}.csv"
    return Path(raw_root) / "Hourly_OD" / f"date-hour-soo-dest-{int(year)}.csv"


def load_hourly_station_ridership_for_month(year, month, path=None):
    """Aggregate hourly OD rows to monthly station entries by origin station.

    Args:
        year: Four-digit ridership year.
        month: One-based ridership month.
        path: Optional explicit hourly OD CSV path.

    Returns:
        DataFrame with `Entry Station`, `Ridership`, and `Full Station Name`.
    """
    year = int(year)
    month = int(month)
    csv_path = Path(path) if path else hourly_od_path_for_year(year)
    if not csv_path.exists() and year == 2018:
        csv_path = HOURLY_OD_2018_CSV
    if not csv_path.exists():
        raise FileNotFoundError(f"No hourly OD CSV found for {year}: {csv_path}")

    period_prefix = f"{year}-{month:02d}"
    monthly_chunks = []

    for chunk in pd.read_csv(
        csv_path,
        header=None,
        names=HOURLY_OD_COLUMNS,
        usecols=["Date", "ORIGIN", "Number of Exits"],
        chunksize=500_000,
    ):
        selected_rows = chunk[chunk["Date"].astype(str).str.startswith(period_prefix)].copy()
        if selected_rows.empty:
            continue

        selected_rows["Entry Station"] = selected_rows["ORIGIN"].map(normalize_hourly_od_station_code)
        selected_rows = selected_rows[selected_rows["Entry Station"].isin(STATION_MAPPING)]
        selected_rows["Ridership"] = pd.to_numeric(
            selected_rows["Number of Exits"],
            errors="coerce",
        ).fillna(0)
        monthly_chunks.append(
            selected_rows.groupby("Entry Station", as_index=False)["Ridership"].sum()
        )

    if monthly_chunks:
        station_ridership = (
            pd.concat(monthly_chunks, ignore_index=True)
            .groupby("Entry Station", as_index=False)["Ridership"]
            .sum()
        )
    else:
        station_ridership = pd.DataFrame(columns=["Entry Station", "Ridership"])

    station_ridership["Full Station Name"] = station_ridership["Entry Station"].map(STATION_MAPPING)
    return station_ridership


def discover_hourly_od_files(raw_root=RAW):
    """Return hourly OD CSVs discovered under `data/raw/Hourly_OD`."""
    return tuple(sorted((Path(raw_root) / "Hourly_OD").glob("date-hour-soo-dest-*.csv")))


def parse_hourly_od_year(path):
    """Extract the four-digit year from a `date-hour-soo-dest-YYYY.csv` path."""
    match = re.search(r"date-hour-soo-dest-(\d{4})\.csv$", Path(path).name)
    if not match:
        raise ValueError(f"Could not infer hourly OD year from {path}.")
    return int(match.group(1))


def build_hourly_station_monthly_summaries(raw_root=RAW):
    """Aggregate discovered hourly OD CSVs to monthly station-entry summaries."""
    monthly_summaries = []

    for hourly_file in discover_hourly_od_files(raw_root):
        hourly_summary = _aggregate_hourly_od_file_to_monthly_stations(hourly_file)
        if hourly_summary.empty:
            continue
        monthly_summaries.append(hourly_summary)

    if not monthly_summaries:
        return pd.DataFrame(
            columns=[
                "Year",
                "Month",
                "Entry Station",
                "Ridership",
                "Full Station Name",
                "Period",
                "Source Type",
                "Source File",
            ]
        )

    return pd.concat(monthly_summaries, ignore_index=True)


def write_hourly_station_monthly_summaries(
    output_dir=None,
    summary_path=None,
    raw_root=RAW,
):
    """Write hourly-derived monthly station summaries by month and combined CSV."""
    output_dir = Path(output_dir) if output_dir else PROCESSED_HOURLY_STATION_MONTHLY_DIR
    summary_path = Path(summary_path) if summary_path else PROCESSED_HOURLY_STATION_MONTHLY_SUMMARY
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    hourly_summary = build_hourly_station_monthly_summaries(raw_root)
    hourly_summary.to_csv(summary_path, index=False)

    for (year, month), monthly_summary in hourly_summary.groupby(["Year", "Month"]):
        monthly_path = output_dir / f"{int(year)}-{int(month):02d}.csv"
        monthly_summary.to_csv(monthly_path, index=False)

    return summary_path


def build_hourly_workbook_validation(raw_root=RAW, hourly_summary=None):
    """Compare hourly OD station totals against workbook Total Trips values."""
    if hourly_summary is None:
        hourly_summary = build_hourly_station_monthly_summaries(raw_root)

    validation_rows = []
    for workbook in discover_ridership_workbooks(raw_root):
        year, month = parse_ridership_period(workbook)
        if not _workbook_has_total_trips_sheet(workbook):
            continue

        workbook_summary = summarize_station_ridership(load_monthly_ridership_long(workbook))
        hourly_month = hourly_summary[
            (hourly_summary["Year"].astype(int) == year)
            & (hourly_summary["Month"].astype(int) == month)
        ]
        if hourly_month.empty:
            continue

        comparison = workbook_summary.merge(
            hourly_month[["Entry Station", "Ridership", "Source File"]],
            on="Entry Station",
            how="outer",
            suffixes=(" Workbook", " Hourly"),
        )
        comparison["Year"] = year
        comparison["Month"] = month
        comparison["Period"] = f"{year}-{month:02d}"
        comparison["Entry Station"] = comparison["Entry Station"].astype("object")
        comparison["Full Station Name"] = comparison["Full Station Name"].astype("object")
        comparison["Full Station Name"] = comparison["Full Station Name"].fillna(
            comparison["Entry Station"].map(STATION_MAPPING)
        )
        comparison["Workbook Ridership"] = pd.to_numeric(
            comparison["Ridership Workbook"],
            errors="coerce",
        ).fillna(0)
        comparison["Hourly OD Ridership"] = pd.to_numeric(
            comparison["Ridership Hourly"],
            errors="coerce",
        ).fillna(0)
        comparison["Difference"] = comparison["Hourly OD Ridership"] - comparison["Workbook Ridership"]
        comparison["Absolute Difference"] = comparison["Difference"].abs()
        comparison["Has Difference"] = comparison["Absolute Difference"] > 0
        comparison["Workbook Source File"] = str(workbook)
        comparison["Hourly Source File"] = comparison["Source File"]
        validation_rows.append(
            comparison[
                [
                    "Year",
                    "Month",
                    "Period",
                    "Entry Station",
                    "Full Station Name",
                    "Workbook Ridership",
                    "Hourly OD Ridership",
                    "Difference",
                    "Absolute Difference",
                    "Has Difference",
                    "Workbook Source File",
                    "Hourly Source File",
                ]
            ]
        )

    if not validation_rows:
        return pd.DataFrame(
            columns=[
                "Year",
                "Month",
                "Period",
                "Entry Station",
                "Full Station Name",
                "Workbook Ridership",
                "Hourly OD Ridership",
                "Difference",
                "Absolute Difference",
                "Has Difference",
                "Workbook Source File",
                "Hourly Source File",
            ]
        )

    return pd.concat(validation_rows, ignore_index=True)


def write_hourly_workbook_validation(output_path=None, raw_root=RAW, hourly_summary_path=None):
    """Write station-level hourly OD versus workbook Total Trips validation."""
    output_path = Path(output_path) if output_path else PROCESSED_HOURLY_VALIDATION_CSV
    hourly_summary_path = (
        Path(hourly_summary_path)
        if hourly_summary_path
        else PROCESSED_HOURLY_STATION_MONTHLY_SUMMARY
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)

    hourly_summary = (
        pd.read_csv(hourly_summary_path)
        if hourly_summary_path.exists()
        else build_hourly_station_monthly_summaries(raw_root)
    )
    validation = build_hourly_workbook_validation(raw_root, hourly_summary)
    validation.to_csv(output_path, index=False)
    return output_path


def _aggregate_hourly_od_file_to_monthly_stations(hourly_file):
    """Aggregate one hourly OD file into monthly station-entry ridership."""
    hourly_file = Path(hourly_file)
    monthly_chunks = []

    for chunk in pd.read_csv(
        hourly_file,
        header=None,
        names=HOURLY_OD_COLUMNS,
        usecols=["Date", "ORIGIN", "Number of Exits"],
        chunksize=500_000,
    ):
        chunk["Entry Station"] = chunk["ORIGIN"].map(normalize_hourly_od_station_code)
        chunk = chunk[chunk["Entry Station"].isin(STATION_MAPPING)].copy()
        if chunk.empty:
            continue

        chunk["Date"] = pd.to_datetime(chunk["Date"], errors="coerce")
        chunk = chunk.dropna(subset=["Date"])
        chunk["Year"] = chunk["Date"].dt.year.astype(int)
        chunk["Month"] = chunk["Date"].dt.month.astype(int)
        chunk["Ridership"] = pd.to_numeric(
            chunk["Number of Exits"],
            errors="coerce",
        ).fillna(0)
        monthly_chunks.append(
            chunk.groupby(["Year", "Month", "Entry Station"], as_index=False)["Ridership"].sum()
        )

    if monthly_chunks:
        hourly_summary = (
            pd.concat(monthly_chunks, ignore_index=True)
            .groupby(["Year", "Month", "Entry Station"], as_index=False)["Ridership"]
            .sum()
        )
    else:
        hourly_summary = pd.DataFrame(columns=["Year", "Month", "Entry Station", "Ridership"])

    hourly_summary["Full Station Name"] = hourly_summary["Entry Station"].map(STATION_MAPPING)
    hourly_summary["Period"] = (
        hourly_summary["Year"].astype(str)
        + "-"
        + hourly_summary["Month"].astype(int).astype(str).str.zfill(2)
    )
    hourly_summary["Source Type"] = "hourly_od"
    hourly_summary["Source File"] = str(hourly_file)
    return hourly_summary[
        [
            "Year",
            "Month",
            "Entry Station",
            "Ridership",
            "Full Station Name",
            "Period",
            "Source Type",
            "Source File",
        ]
    ]


def build_monthly_station_ridership_summaries(raw_root=RAW, hourly_summary=None):
    """Build app-ready monthly station summaries from the best available source.

    Hourly OD is the canonical station ridership source where available because
    it supports future hourly/daily views. Monthly workbooks are used only for
    periods without hourly OD coverage.

    Args:
        raw_root: Directory containing `ridership_*` raw workbook folders.
        hourly_summary: Optional precomputed hourly monthly station summary.

    Returns:
        DataFrame with one row per year/month/entry-station combination.
    """
    monthly_summaries = []
    if hourly_summary is None:
        hourly_summary = build_hourly_station_monthly_summaries(raw_root)

    hourly_periods = set()
    if not hourly_summary.empty:
        hourly_summary = hourly_summary.copy()
        hourly_summary["Year"] = hourly_summary["Year"].astype(int)
        hourly_summary["Month"] = hourly_summary["Month"].astype(int)
        hourly_periods = {
            (int(row["Year"]), int(row["Month"]))
            for _, row in hourly_summary[["Year", "Month"]].drop_duplicates().iterrows()
        }
        monthly_summaries.append(
            hourly_summary[
                [
                    "Year",
                    "Month",
                    "Entry Station",
                    "Ridership",
                    "Full Station Name",
                    "Period",
                    "Source Type",
                    "Source File",
                ]
            ]
        )

    for workbook in discover_ridership_workbooks(raw_root):
        year, month = parse_ridership_period(workbook)
        if (year, month) in hourly_periods:
            continue
        if not _workbook_has_total_trips_sheet(workbook):
            continue

        station_summary = summarize_station_ridership(load_monthly_ridership_long(workbook))
        source_type = "workbook_total_trips"
        source_file = str(workbook)

        station_summary.insert(0, "Month", month)
        station_summary.insert(0, "Year", year)
        station_summary["Period"] = f"{year}-{month:02d}"
        station_summary["Source Type"] = source_type
        station_summary["Source File"] = source_file
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
                "Source Type",
                "Source File",
            ]
        )

    return pd.concat(monthly_summaries, ignore_index=True)


def write_monthly_station_ridership_summaries(
    output_path=None,
    raw_root=RAW,
    hourly_summary=None,
    hourly_summary_path=None,
):
    """Write processed monthly station ridership summaries to CSV.

    Args:
        output_path: Optional destination. Defaults to
            `data/processed/station_ridership_monthly_summary.csv`.
        raw_root: Directory containing raw ridership workbooks.
        hourly_summary: Optional precomputed hourly monthly station summary.
        hourly_summary_path: Optional path to a precomputed hourly summary CSV.

    Returns:
        Path to the written summary CSV.
    """
    output_path = Path(output_path) if output_path else monthly_station_summary_path()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if hourly_summary is None and hourly_summary_path is not None:
        hourly_summary = pd.read_csv(hourly_summary_path)
    build_monthly_station_ridership_summaries(raw_root, hourly_summary).to_csv(output_path, index=False)
    return output_path


def monthly_station_summary_path():
    """Return the standard processed monthly ridership summary path."""
    return PROCESSED / "station_ridership_monthly_summary.csv"


def load_monthly_station_ridership_summary(processed_path=None):
    """Load the processed monthly station ridership summary.

    Runtime code should use this processed artifact instead of parsing raw
    workbooks or hourly OD CSVs. Regenerate it with `scripts/prepare_data.py`.
    """
    processed_path = Path(processed_path) if processed_path else monthly_station_summary_path()
    if not processed_path.exists():
        raise FileNotFoundError(
            f"Processed ridership summary not found: {processed_path}. "
            "Run `python scripts/prepare_data.py` to regenerate it."
        )

    monthly_summary = pd.read_csv(processed_path)
    monthly_summary["Year"] = monthly_summary["Year"].astype(int)
    monthly_summary["Month"] = monthly_summary["Month"].astype(int)
    monthly_summary["Ridership"] = pd.to_numeric(
        monthly_summary["Ridership"],
        errors="coerce",
    ).fillna(0)
    return monthly_summary


def load_station_ridership_for_period(year, month, processed_path=None):
    """Load processed station ridership for one selected year/month.

    The app consumes `data/processed/station_ridership_monthly_summary.csv`.
    Raw parsing belongs to `scripts/prepare_data.py`.
    """
    year = int(year)
    month = int(month)
    monthly_summary = load_monthly_station_ridership_summary(processed_path)
    selected_summary = monthly_summary[
        (monthly_summary["Year"] == year)
        & (monthly_summary["Month"] == month)
    ].copy()
    if not selected_summary.empty:
        return _normalize_station_summary_columns(selected_summary)

    raise ValueError(f"Processed ridership summary has no rows for {year}-{month:02d}.")


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


def get_available_ridership_periods(processed_path=None):
    """Return processed ridership periods as `(year, month)` tuples."""
    monthly_summary = load_monthly_station_ridership_summary(processed_path)
    return tuple(
        sorted(
            {
                (int(row["Year"]), int(row["Month"]))
                for _, row in monthly_summary[["Year", "Month"]].drop_duplicates().iterrows()
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


def _workbook_has_total_trips_sheet(file_path):
    """Return whether a workbook has a comparable monthly total trips sheet."""
    excel_file = pd.ExcelFile(file_path)
    try:
        return any(sheet in excel_file.sheet_names for sheet in ("Total Trips", "Total Trips OD"))
    finally:
        excel_file.close()


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
    return normalize_workbook_station_code(str(name))


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
