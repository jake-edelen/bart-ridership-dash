"""Runtime data loading and shaping for the BART dashboard.

This module prepares the in-memory tables consumed by `app.py` and exposes
writer helpers used by `scripts/prepare_data.py`. The Dash app should call the
loaders only; processed artifact generation remains a separate script step.
"""

# region Imports
from dataclasses import dataclass
from calendar import monthrange
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
    PROCESSED_HOURLY_COMPLETENESS_CSV,
    PROCESSED_HOURLY_COMPLETENESS_PARQUET,
    PROCESSED_HOURLY_STATION_MONTHLY_DIR,
    PROCESSED_HOURLY_STATION_MONTHLY_PARQUET_DIR,
    PROCESSED_HOURLY_STATION_MONTHLY_SUMMARY,
    PROCESSED_HOURLY_STATION_MONTHLY_SUMMARY_PARQUET,
    PROCESSED_HOURLY_VALIDATION_CSV,
    PROCESSED_HOURLY_VALIDATION_PARQUET,
    PROCESSED_RAW_ROUTES_GEOJSON,
    PROCESSED_STATIONS_GEOJSON,
    PROCESSED_STATION_RIDERSHIP_MONTHLY_SUMMARY_CSV,
    PROCESSED_STATION_RIDERSHIP_MONTHLY_SUMMARY_PARQUET,
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
from .station_registry import add_station_identity
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

STATION_SUMMARY_COLUMNS = [
    "Year",
    "Month",
    "Period",
    "Entry Station",
    "station_id",
    "Full Station Name",
    "display_name",
    "Ridership",
    "Source Type",
    "Source File",
]

VALIDATION_COLUMNS = [
    "Year",
    "Month",
    "Period",
    "Entry Station",
    "station_id",
    "Full Station Name",
    "display_name",
    "Workbook Ridership",
    "Hourly OD Ridership",
    "Difference",
    "Absolute Difference",
    "Has Difference",
    "Workbook Source File",
    "Hourly Source File",
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
    """Load and assemble processed datasets required by the Dash app.

    Runtime app loading intentionally reads processed/reference artifacts only.
    Raw workbook, hourly OD, KML, and GDB reads belong to `scripts/prepare_data.py`
    and the writer helpers in this module.

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
    return _with_station_identity(station_ridership)


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

    monthly_summary = _aggregate_hourly_od_file_to_monthly_stations(csv_path)
    selected_summary = monthly_summary[
        (monthly_summary["Year"].astype(int) == year)
        & (monthly_summary["Month"].astype(int) == month)
    ].copy()

    if selected_summary.empty:
        station_ridership = pd.DataFrame(columns=["Entry Station", "Ridership"])
        station_ridership["Full Station Name"] = station_ridership["Entry Station"].map(STATION_MAPPING)
        return _with_station_identity(station_ridership)

    return _normalize_station_summary_columns(selected_summary)


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
        return pd.DataFrame(columns=STATION_SUMMARY_COLUMNS)

    return pd.concat(monthly_summaries, ignore_index=True)[STATION_SUMMARY_COLUMNS]


def write_hourly_station_monthly_summaries(
    output_dir=None,
    summary_path=None,
    raw_root=RAW,
):
    """Write hourly-derived monthly station summaries as CSV and Parquet."""
    output_dir = Path(output_dir) if output_dir else PROCESSED_HOURLY_STATION_MONTHLY_DIR
    summary_path = Path(summary_path) if summary_path else PROCESSED_HOURLY_STATION_MONTHLY_SUMMARY
    parquet_dir = PROCESSED_HOURLY_STATION_MONTHLY_PARQUET_DIR
    parquet_summary_path = PROCESSED_HOURLY_STATION_MONTHLY_SUMMARY_PARQUET
    output_dir.mkdir(parents=True, exist_ok=True)
    parquet_dir.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    hourly_summary = build_hourly_station_monthly_summaries(raw_root)
    _write_table_artifacts(hourly_summary, summary_path, parquet_summary_path)

    for (year, month), monthly_summary in hourly_summary.groupby(["Year", "Month"]):
        monthly_path = output_dir / f"{int(year)}-{int(month):02d}.csv"
        monthly_parquet_path = parquet_dir / f"{int(year)}-{int(month):02d}.parquet"
        _write_table_artifacts(monthly_summary, monthly_path, monthly_parquet_path)

    return summary_path


def build_hourly_completeness_audit(raw_root=RAW, validation=None, hourly_summary=None):
    """Build period-level quality flags for hourly OD coverage and imputation.

    Missing full calendar days are imputed at the hourly station level before
    monthly summaries are written. This audit records where that happened and
    where workbook totals disagree with the canonical hourly-derived totals.
    """
    quality_rows = []
    validation_totals = _period_validation_totals(validation)
    hourly_totals = _period_hourly_summary_totals(hourly_summary)

    for hourly_file in discover_hourly_od_files(raw_root):
        daily_hourly = _read_hourly_od_file_to_daily_hourly_stations(hourly_file)
        if daily_hourly.empty:
            continue

        imputed_daily_hourly = _impute_missing_full_days(daily_hourly)
        raw_period_totals = (
            daily_hourly.groupby(["Year", "Month"], as_index=False)["Ridership"]
            .sum()
            .rename(columns={"Ridership": "Raw Hourly Ridership"})
        )
        imputed_period_totals = (
            imputed_daily_hourly.groupby(["Year", "Month"], as_index=False)["Ridership"]
            .sum()
            .rename(columns={"Ridership": "Imputed Hourly Ridership"})
        )
        raw_and_imputed = raw_period_totals.merge(
            imputed_period_totals,
            on=["Year", "Month"],
            how="outer",
        )

        for (year, month), period_rows in daily_hourly.groupby(["Year", "Month"]):
            period = f"{int(year)}-{int(month):02d}"
            daily_hour_counts = period_rows.groupby("Date")["Hour"].nunique()
            daily_first_hours = period_rows.groupby("Date")["Hour"].min()
            daily_last_hours = period_rows.groupby("Date")["Hour"].max()
            expected_dates = _expected_period_dates(year, month)
            observed_dates = set(pd.to_datetime(period_rows["Date"]).dt.normalize())
            missing_dates = [
                date
                for date in expected_dates
                if date not in observed_dates
            ]
            totals = raw_and_imputed[
                (raw_and_imputed["Year"].astype(int) == int(year))
                & (raw_and_imputed["Month"].astype(int) == int(month))
            ].iloc[0]
            workbook_total = validation_totals.get(period, {}).get("Workbook Ridership")
            validation_hourly_total = validation_totals.get(period, {}).get("Hourly OD Ridership")
            canonical_hourly_total = hourly_totals.get(period, validation_hourly_total)
            difference = (
                None
                if workbook_total is None or canonical_hourly_total is None
                else canonical_hourly_total - workbook_total
            )
            difference_pct = (
                None
                if difference is None or not workbook_total
                else difference / workbook_total
            )

            quality_rows.append(
                {
                    "Period": period,
                    "Year": int(year),
                    "Month": int(month),
                    "Source File": str(hourly_file),
                    "Calendar Days": len(expected_dates),
                    "Observed Days": int(daily_hour_counts.shape[0]),
                    "Missing Full Days": len(missing_dates),
                    "Missing Full Day List": ", ".join(date.strftime("%Y-%m-%d") for date in missing_dates),
                    "Min Present Hours": int(daily_hour_counts.min()),
                    "Median Present Hours": float(daily_hour_counts.median()),
                    "Median First Hour": float(daily_first_hours.median()),
                    "Median Last Hour": float(daily_last_hours.median()),
                    "Raw Hourly Ridership": float(totals["Raw Hourly Ridership"]),
                    "Imputed Hourly Ridership": float(totals["Imputed Hourly Ridership"]),
                    "Imputed Ridership Added": float(
                        totals["Imputed Hourly Ridership"] - totals["Raw Hourly Ridership"]
                    ),
                    "Workbook Ridership": workbook_total,
                    "Canonical Hourly Ridership": canonical_hourly_total,
                    "Workbook Difference": difference,
                    "Workbook Difference Pct": difference_pct,
                    "Quality Flag": _hourly_quality_flag(
                        missing_full_days=len(missing_dates),
                        min_present_hours=float(daily_hour_counts.min()),
                        median_present_hours=float(daily_hour_counts.median()),
                        workbook_difference_pct=difference_pct,
                    ),
                }
            )

    if not quality_rows:
        return pd.DataFrame(
            columns=[
                "Period",
                "Year",
                "Month",
                "Quality Flag",
            ]
        )

    return pd.DataFrame(quality_rows).sort_values(["Year", "Month"]).reset_index(drop=True)


def write_hourly_completeness_audit(
    output_path=None,
    raw_root=RAW,
    validation_path=None,
    hourly_summary_path=None,
):
    """Write the hourly OD completeness/imputation quality artifact."""
    output_path = Path(output_path) if output_path else PROCESSED_HOURLY_COMPLETENESS_CSV
    parquet_path = PROCESSED_HOURLY_COMPLETENESS_PARQUET
    validation_path = Path(validation_path) if validation_path else PROCESSED_HOURLY_VALIDATION_CSV
    hourly_summary_path = (
        Path(hourly_summary_path)
        if hourly_summary_path
        else PROCESSED_HOURLY_STATION_MONTHLY_SUMMARY
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)

    validation = (
        _read_processed_table(validation_path, PROCESSED_HOURLY_VALIDATION_PARQUET)
        if validation_path.exists() or PROCESSED_HOURLY_VALIDATION_PARQUET.exists()
        else None
    )
    hourly_summary = (
        _read_processed_table(hourly_summary_path, PROCESSED_HOURLY_STATION_MONTHLY_SUMMARY_PARQUET)
        if hourly_summary_path.exists() or PROCESSED_HOURLY_STATION_MONTHLY_SUMMARY_PARQUET.exists()
        else None
    )
    audit = build_hourly_completeness_audit(
        raw_root=raw_root,
        validation=validation,
        hourly_summary=hourly_summary,
    )
    _write_table_artifacts(audit, output_path, parquet_path)
    return output_path


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

        comparison = workbook_summary.drop(
            columns=["station_id", "display_name"],
            errors="ignore",
        ).merge(
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
        comparison = _with_station_identity(comparison)
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
        validation_rows.append(comparison[VALIDATION_COLUMNS])

    if not validation_rows:
        return pd.DataFrame(columns=VALIDATION_COLUMNS)

    return pd.concat(validation_rows, ignore_index=True)


def write_hourly_workbook_validation(output_path=None, raw_root=RAW, hourly_summary_path=None):
    """Write station-level hourly OD versus workbook Total Trips validation."""
    output_path = Path(output_path) if output_path else PROCESSED_HOURLY_VALIDATION_CSV
    parquet_path = PROCESSED_HOURLY_VALIDATION_PARQUET
    hourly_summary_path = (
        Path(hourly_summary_path)
        if hourly_summary_path
        else PROCESSED_HOURLY_STATION_MONTHLY_SUMMARY
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)

    hourly_summary = (
        _read_processed_table(hourly_summary_path, PROCESSED_HOURLY_STATION_MONTHLY_SUMMARY_PARQUET)
        if hourly_summary_path.exists() or PROCESSED_HOURLY_STATION_MONTHLY_SUMMARY_PARQUET.exists()
        else build_hourly_station_monthly_summaries(raw_root)
    )
    validation = build_hourly_workbook_validation(raw_root, hourly_summary)
    _write_table_artifacts(validation, output_path, parquet_path)
    return output_path


def _aggregate_hourly_od_file_to_monthly_stations(hourly_file):
    """Aggregate one hourly OD file into monthly station-entry ridership."""
    hourly_file = Path(hourly_file)
    daily_hourly = _read_hourly_od_file_to_daily_hourly_stations(hourly_file)
    daily_hourly = _impute_missing_full_days(daily_hourly)

    if not daily_hourly.empty:
        hourly_summary = (
            daily_hourly.groupby(["Year", "Month", "Entry Station"], as_index=False)["Ridership"]
            .sum()
        )
    else:
        hourly_summary = pd.DataFrame(columns=["Year", "Month", "Entry Station", "Ridership"])

    hourly_summary["Full Station Name"] = hourly_summary["Entry Station"].map(STATION_MAPPING)
    hourly_summary = _with_station_identity(hourly_summary)
    hourly_summary["Period"] = (
        hourly_summary["Year"].astype(str)
        + "-"
        + hourly_summary["Month"].astype(int).astype(str).str.zfill(2)
    )
    hourly_summary["Source Type"] = "hourly_od"
    hourly_summary["Source File"] = str(hourly_file)
    return hourly_summary[STATION_SUMMARY_COLUMNS]


def _read_hourly_od_file_to_daily_hourly_stations(hourly_file):
    """Read one hourly OD CSV into station/hour totals for imputation/audit."""
    hourly_file = Path(hourly_file)
    hourly_chunks = []

    for chunk in pd.read_csv(
        hourly_file,
        header=None,
        names=HOURLY_OD_COLUMNS,
        usecols=["Date", "Hour", "ORIGIN", "Number of Exits"],
        chunksize=500_000,
    ):
        chunk["Entry Station"] = chunk["ORIGIN"].map(normalize_hourly_od_station_code)
        chunk = chunk[chunk["Entry Station"].isin(STATION_MAPPING)].copy()
        if chunk.empty:
            continue

        chunk["Date"] = pd.to_datetime(chunk["Date"], errors="coerce").dt.normalize()
        chunk["Hour"] = pd.to_numeric(chunk["Hour"], errors="coerce")
        chunk["Ridership"] = pd.to_numeric(
            chunk["Number of Exits"],
            errors="coerce",
        ).fillna(0)
        chunk = chunk.dropna(subset=["Date", "Hour"])
        chunk["Hour"] = chunk["Hour"].astype(int)
        chunk["Year"] = chunk["Date"].dt.year.astype(int)
        chunk["Month"] = chunk["Date"].dt.month.astype(int)
        hourly_chunks.append(
            chunk.groupby(
                ["Date", "Year", "Month", "Hour", "Entry Station"],
                as_index=False,
            )["Ridership"].sum()
        )

    if not hourly_chunks:
        return pd.DataFrame(
            columns=["Date", "Year", "Month", "Hour", "Entry Station", "Ridership"]
        )

    return (
        pd.concat(hourly_chunks, ignore_index=True)
        .groupby(["Date", "Year", "Month", "Hour", "Entry Station"], as_index=False)["Ridership"]
        .sum()
    )


def _impute_missing_full_days(daily_hourly):
    """Fill missing full calendar days using same-month hourly station means."""
    if daily_hourly.empty:
        return daily_hourly.copy()

    imputed_parts = [daily_hourly.copy()]

    for (year, month), period_rows in daily_hourly.groupby(["Year", "Month"]):
        expected_dates = _expected_period_dates(year, month)
        observed_dates = set(pd.to_datetime(period_rows["Date"]).dt.normalize())
        missing_dates = [
            date
            for date in expected_dates
            if date not in observed_dates
        ]
        if not missing_dates:
            continue

        source_rows = period_rows.copy()
        source_rows["Weekday"] = source_rows["Date"].dt.dayofweek
        source_rows["Day Type"] = source_rows["Weekday"].map(
            lambda weekday: "weekday" if weekday < 5 else "weekend"
        )
        weekday_means = source_rows.groupby(
            ["Weekday", "Hour", "Entry Station"],
            as_index=False,
        )["Ridership"].mean().rename(columns={"Ridership": "Weekday Mean"})
        day_type_means = source_rows.groupby(
            ["Day Type", "Hour", "Entry Station"],
            as_index=False,
        )["Ridership"].mean().rename(columns={"Ridership": "Day Type Mean"})
        hour_means = source_rows.groupby(
            ["Hour", "Entry Station"],
            as_index=False,
        )["Ridership"].mean().rename(columns={"Ridership": "Hour Mean"})
        stations = sorted(source_rows["Entry Station"].dropna().unique())
        hours = list(range(24))
        templates = []

        for missing_date in missing_dates:
            template = pd.MultiIndex.from_product(
                [[missing_date], hours, stations],
                names=["Date", "Hour", "Entry Station"],
            ).to_frame(index=False)
            template["Year"] = int(year)
            template["Month"] = int(month)
            template["Weekday"] = missing_date.dayofweek
            template["Day Type"] = "weekday" if missing_date.dayofweek < 5 else "weekend"
            templates.append(template)

        imputed_rows = pd.concat(templates, ignore_index=True)
        imputed_rows = imputed_rows.merge(
            weekday_means,
            on=["Weekday", "Hour", "Entry Station"],
            how="left",
        ).merge(
            day_type_means,
            on=["Day Type", "Hour", "Entry Station"],
            how="left",
        ).merge(
            hour_means,
            on=["Hour", "Entry Station"],
            how="left",
        )
        imputed_rows["Ridership"] = (
            imputed_rows["Weekday Mean"]
            .combine_first(imputed_rows["Day Type Mean"])
            .combine_first(imputed_rows["Hour Mean"])
            .fillna(0)
        )
        imputed_parts.append(
            imputed_rows[["Date", "Year", "Month", "Hour", "Entry Station", "Ridership"]]
        )

    return pd.concat(imputed_parts, ignore_index=True)


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
        if {"station_id", "display_name"} - set(hourly_summary.columns):
            hourly_summary = _with_station_identity(hourly_summary)
        hourly_summary["Year"] = hourly_summary["Year"].astype(int)
        hourly_summary["Month"] = hourly_summary["Month"].astype(int)
        hourly_periods = {
            (int(row["Year"]), int(row["Month"]))
            for _, row in hourly_summary[["Year", "Month"]].drop_duplicates().iterrows()
        }
        monthly_summaries.append(hourly_summary[STATION_SUMMARY_COLUMNS])

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
        return pd.DataFrame(columns=STATION_SUMMARY_COLUMNS)

    return pd.concat(monthly_summaries, ignore_index=True)[STATION_SUMMARY_COLUMNS]


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
    parquet_path = monthly_station_summary_parquet_path()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if hourly_summary is None and hourly_summary_path is not None:
        hourly_summary = _read_processed_table(
            hourly_summary_path,
            PROCESSED_HOURLY_STATION_MONTHLY_SUMMARY_PARQUET,
        )
    monthly_summary = build_monthly_station_ridership_summaries(raw_root, hourly_summary)
    _write_table_artifacts(monthly_summary, output_path, parquet_path)
    return output_path


def monthly_station_summary_path():
    """Return the standard processed monthly ridership CSV path."""
    return PROCESSED_STATION_RIDERSHIP_MONTHLY_SUMMARY_CSV


def monthly_station_summary_parquet_path():
    """Return the standard processed monthly ridership Parquet path."""
    return PROCESSED_STATION_RIDERSHIP_MONTHLY_SUMMARY_PARQUET


def load_monthly_station_ridership_summary(processed_path=None):
    """Load the processed monthly station ridership summary.

    Runtime code should use this processed artifact instead of parsing raw
    workbooks or hourly OD CSVs. Regenerate it with `scripts/prepare_data.py`.
    """
    processed_path = Path(processed_path) if processed_path else monthly_station_summary_path()
    parquet_path = monthly_station_summary_parquet_path() if processed_path == monthly_station_summary_path() else None
    if not processed_path.exists() and not (parquet_path and parquet_path.exists()):
        raise FileNotFoundError(
            f"Processed ridership summary not found: {processed_path}. "
            "Run `python scripts/prepare_data.py` to regenerate it."
        )

    monthly_summary = _read_processed_table(processed_path, parquet_path)
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
    ridership_name_column = (
        "display_name"
        if "display_name" in station_ridership.columns
        else "Full Station Name"
    )
    stations_with_ridership = stations_gdf.merge(
        station_ridership,
        left_on="Name",
        right_on=ridership_name_column,
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
def _with_station_identity(station_summary):
    """Attach canonical `station_id` and `display_name` to legacy-code rows."""
    station_summary = station_summary.drop(
        columns=["station_id", "display_name"],
        errors="ignore",
    )
    return add_station_identity(
        station_summary,
        code_column="Entry Station",
        source_system="legacy_app_code",
    )


def _write_table_artifacts(dataframe, csv_path, parquet_path):
    """Write one processed table as CSV compatibility plus Parquet runtime data."""
    csv_path = Path(csv_path)
    parquet_path = Path(parquet_path)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    parquet_path.parent.mkdir(parents=True, exist_ok=True)
    dataframe.to_csv(csv_path, index=False)
    dataframe.to_parquet(parquet_path, index=False)


def _read_processed_table(csv_path, parquet_path=None):
    """Read Parquet when present, falling back to CSV for compatibility."""
    csv_path = Path(csv_path)
    parquet_path = Path(parquet_path) if parquet_path else csv_path.with_suffix(".parquet")

    if parquet_path.exists():
        return pd.read_parquet(parquet_path)
    if csv_path.exists():
        return pd.read_csv(csv_path)
    raise FileNotFoundError(f"No processed table found at {parquet_path} or {csv_path}.")


def _expected_period_dates(year, month):
    """Return every calendar date for one year/month period."""
    year = int(year)
    month = int(month)
    return pd.date_range(
        f"{year}-{month:02d}-01",
        periods=monthrange(year, month)[1],
        freq="D",
    )


def _period_validation_totals(validation):
    """Return workbook/hourly validation totals keyed by `YYYY-MM` period."""
    if validation is None or validation.empty:
        return {}

    totals = (
        validation.groupby("Period", as_index=False)[
            ["Workbook Ridership", "Hourly OD Ridership"]
        ]
        .sum()
    )
    return {
        row["Period"]: {
            "Workbook Ridership": float(row["Workbook Ridership"]),
            "Hourly OD Ridership": float(row["Hourly OD Ridership"]),
        }
        for _, row in totals.iterrows()
    }


def _period_hourly_summary_totals(hourly_summary):
    """Return canonical hourly summary totals keyed by `YYYY-MM` period."""
    if hourly_summary is None or hourly_summary.empty:
        return {}

    summary = hourly_summary.copy()
    if "Period" not in summary.columns:
        summary["Period"] = (
            summary["Year"].astype(int).astype(str)
            + "-"
            + summary["Month"].astype(int).astype(str).str.zfill(2)
        )
    totals = summary.groupby("Period", as_index=False)["Ridership"].sum()
    return {
        row["Period"]: float(row["Ridership"])
        for _, row in totals.iterrows()
    }


def _hourly_quality_flag(
    missing_full_days,
    min_present_hours,
    median_present_hours,
    workbook_difference_pct,
):
    """Classify the period-level hourly source quality issue, if any."""
    if missing_full_days:
        return "missing_full_days_imputed"
    if median_present_hours <= 18:
        return "covid_service_hours"
    if min_present_hours < max(1, median_present_hours * 0.6):
        return "partial_low_hour_day"
    if workbook_difference_pct is not None and abs(workbook_difference_pct) > 0.05:
        return "complete_hours_source_discrepancy"
    return "clean_or_minor_difference"


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

    if {"station_id", "display_name"} - set(station_summary.columns):
        station_summary = _with_station_identity(station_summary)

    return station_summary[
        [
            "Entry Station",
            "station_id",
            "Ridership",
            "Full Station Name",
            "display_name",
        ]
    ]


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
