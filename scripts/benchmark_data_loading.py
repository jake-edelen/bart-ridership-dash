"""Benchmark raw CSV aggregation against processed CSV and Parquet summaries.

This script produces local, machine-specific timing evidence for the data
engineering payoff of pre-aggregated summary tables. It does not run during app
startup.
"""

# region Imports And Path Setup
from argparse import ArgumentParser
from dataclasses import dataclass
from pathlib import Path
from statistics import median
from time import perf_counter
import sys

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.config import (  # noqa: E402
    PROCESSED,
    PROCESSED_STATION_RIDERSHIP_MONTHLY_SUMMARY_CSV,
    PROCESSED_STATION_RIDERSHIP_MONTHLY_SUMMARY_PARQUET,
)
from src.data_loader import (  # noqa: E402
    HOURLY_OD_COLUMNS,
    discover_hourly_od_files,
    hourly_od_path_for_year,
)
from src.station_mapping import STATION_MAPPING, normalize_hourly_od_station_code  # noqa: E402
# endregion


# region Constants And Scenarios
BENCHMARK_OUTPUT_CSV = PROCESSED / "benchmark_data_loading.csv"
DEFAULT_YEAR = 2025
DEFAULT_MONTH = 5


@dataclass(frozen=True)
class BenchmarkScenario:
    """One representative data-access query to benchmark."""

    query: str
    period: str
    years: tuple[int, ...]
    months: tuple[int, ...] | None = None
    station_codes: tuple[str, ...] | None = None


def representative_scenarios():
    """Return a compact suite of representative station-ridership queries."""
    return (
        BenchmarkScenario("one_month_all_stations", "2025-05", (2025,), (5,), None),
        BenchmarkScenario("one_full_year_all_stations", "2025", (2025,), None, None),
        BenchmarkScenario("all_available_years", "2018-2025", _available_hourly_years(), None, None),
        BenchmarkScenario("one_station_month", "2025-05", (2025,), (5,), ("EM",)),
        BenchmarkScenario("high_volume_station_pair_month", "2025-05", (2025,), (5,), ("EM", "MT")),
    )


def single_month_scenario(year, month):
    """Return the historical one-month benchmark scenario."""
    return BenchmarkScenario(
        "one_month_all_stations",
        f"{int(year)}-{int(month):02d}",
        (int(year),),
        (int(month),),
        None,
    )
# endregion


# region Benchmark Operations
def benchmark_raw_hourly_od(scenario):
    """Scan raw hourly OD CSVs and aggregate the requested station summary."""
    rows_read = 0
    monthly_chunks = []
    raw_files = _raw_files_for_scenario(scenario)
    start = perf_counter()

    for raw_file in raw_files:
        for chunk in pd.read_csv(
            raw_file,
            header=None,
            names=HOURLY_OD_COLUMNS,
            usecols=["Date", "ORIGIN", "Number of Exits"],
            chunksize=750_000,
        ):
            rows_read += len(chunk)
            selected_rows = _filter_raw_chunk(chunk, scenario)
            if selected_rows.empty:
                continue

            selected_rows["Entry Station"] = selected_rows["ORIGIN"].map(normalize_hourly_od_station_code)
            selected_rows = selected_rows[selected_rows["Entry Station"].isin(STATION_MAPPING)]
            if scenario.station_codes:
                selected_rows = selected_rows[selected_rows["Entry Station"].isin(scenario.station_codes)]
            if selected_rows.empty:
                continue

            selected_rows["Date"] = pd.to_datetime(selected_rows["Date"], errors="coerce")
            selected_rows = selected_rows.dropna(subset=["Date"])
            selected_rows["Year"] = selected_rows["Date"].dt.year.astype(int)
            selected_rows["Month"] = selected_rows["Date"].dt.month.astype(int)
            selected_rows["Ridership"] = pd.to_numeric(
                selected_rows["Number of Exits"],
                errors="coerce",
            ).fillna(0)
            monthly_chunks.append(
                selected_rows.groupby(
                    ["Year", "Month", "Entry Station"],
                    as_index=False,
                )["Ridership"].sum()
            )

    result = _combine_query_chunks(monthly_chunks)
    elapsed_seconds = perf_counter() - start
    return _operation_result(rows_read, result, elapsed_seconds)


def benchmark_processed_csv(scenario, csv_file):
    """Read the processed CSV summary and filter the requested query."""
    start = perf_counter()
    summary = pd.read_csv(csv_file)
    result = _filter_processed_summary(summary, scenario)
    elapsed_seconds = perf_counter() - start
    return _operation_result(len(summary), result, elapsed_seconds)


def benchmark_processed_parquet(scenario, parquet_file):
    """Read the processed Parquet summary and filter the requested query."""
    start = perf_counter()
    summary = pd.read_parquet(parquet_file)
    result = _filter_processed_summary(summary, scenario)
    elapsed_seconds = perf_counter() - start
    return _operation_result(len(summary), result, elapsed_seconds)
# endregion


# region Filtering Helpers
def _filter_raw_chunk(chunk, scenario):
    """Apply year/month filters to one raw CSV chunk before aggregation."""
    date_text = chunk["Date"].astype(str)
    year_prefixes = tuple(f"{year}-" for year in scenario.years)
    selected_rows = chunk[date_text.str.startswith(year_prefixes)].copy()

    if scenario.months is not None and not selected_rows.empty:
        period_prefixes = tuple(
            f"{year}-{month:02d}"
            for year in scenario.years
            for month in scenario.months
        )
        selected_rows = selected_rows[
            selected_rows["Date"].astype(str).str.startswith(period_prefixes)
        ].copy()

    return selected_rows


def _filter_processed_summary(summary, scenario):
    """Filter a processed monthly station summary to the requested query."""
    result = summary[summary["Year"].astype(int).isin(scenario.years)].copy()

    if scenario.months is not None:
        result = result[result["Month"].astype(int).isin(scenario.months)].copy()
    if scenario.station_codes:
        result = result[result["Entry Station"].isin(scenario.station_codes)].copy()

    return (
        result.groupby(["Year", "Month", "Entry Station"], as_index=False)["Ridership"]
        .sum()
        .sort_values(["Year", "Month", "Entry Station"])
        .reset_index(drop=True)
    )


def _combine_query_chunks(chunks):
    """Combine raw chunk aggregates into one comparable result table."""
    if not chunks:
        return pd.DataFrame(columns=["Year", "Month", "Entry Station", "Ridership"])

    return (
        pd.concat(chunks, ignore_index=True)
        .groupby(["Year", "Month", "Entry Station"], as_index=False)["Ridership"]
        .sum()
        .sort_values(["Year", "Month", "Entry Station"])
        .reset_index(drop=True)
    )


def _operation_result(rows_read, result, elapsed_seconds):
    """Return standard benchmark operation metrics."""
    return {
        "rows_read": rows_read,
        "result_rows": len(result),
        "ridership_total": float(result["Ridership"].sum()) if not result.empty else 0.0,
        "elapsed_seconds": elapsed_seconds,
    }
# endregion


# region Result Assembly
def run_benchmarks(year, month, runs, suite, warmup=True):
    """Run benchmark workflows and return one summary row per query/workflow."""
    scenarios = (
        representative_scenarios()
        if suite == "representative"
        else (single_month_scenario(year, month),)
    )
    if warmup:
        _warm_processed_readers()

    workflows = (
        (
            "raw_hourly_od",
            lambda scenario: _raw_files_for_scenario(scenario),
            benchmark_raw_hourly_od,
        ),
        (
            "processed_csv_summary",
            lambda _scenario: (PROCESSED_STATION_RIDERSHIP_MONTHLY_SUMMARY_CSV,),
            lambda scenario: benchmark_processed_csv(
                scenario,
                PROCESSED_STATION_RIDERSHIP_MONTHLY_SUMMARY_CSV,
            ),
        ),
        (
            "processed_parquet_summary",
            lambda _scenario: (PROCESSED_STATION_RIDERSHIP_MONTHLY_SUMMARY_PARQUET,),
            lambda scenario: benchmark_processed_parquet(
                scenario,
                PROCESSED_STATION_RIDERSHIP_MONTHLY_SUMMARY_PARQUET,
            ),
        ),
    )

    rows = []
    for scenario in scenarios:
        for workflow, input_paths, operation in workflows:
            file_paths = tuple(input_paths(scenario))
            missing_paths = [path for path in file_paths if not path.exists()]
            if missing_paths:
                rows.append(_missing_file_row(workflow, scenario, file_paths, runs))
                continue

            run_results = [operation(scenario) for _ in range(runs)]
            elapsed_values = [result["elapsed_seconds"] for result in run_results]
            representative = run_results[elapsed_values.index(median(elapsed_values))]
            rows.append(_benchmark_row(workflow, scenario, file_paths, runs, representative, elapsed_values))

    benchmark_df = pd.DataFrame(rows)
    benchmark_df["speedup_vs_raw"] = pd.NA
    for query, query_rows in benchmark_df.groupby("query"):
        raw_elapsed = query_rows.loc[
            query_rows["workflow"] == "raw_hourly_od",
            "median_elapsed_seconds",
        ]
        if raw_elapsed.empty or pd.isna(raw_elapsed.iloc[0]) or raw_elapsed.iloc[0] <= 0:
            continue
        benchmark_df.loc[benchmark_df["query"] == query, "speedup_vs_raw"] = (
            raw_elapsed.iloc[0]
            / benchmark_df.loc[benchmark_df["query"] == query, "median_elapsed_seconds"]
        )

    return benchmark_df


def _benchmark_row(workflow, scenario, file_paths, runs, representative, elapsed_values):
    """Build one benchmark result row."""
    return {
        "query": scenario.query,
        "workflow": workflow,
        "period": scenario.period,
        "station_codes": ", ".join(scenario.station_codes or ()),
        "runs": int(runs),
        "file_path": "; ".join(path.relative_to(PROJECT_ROOT).as_posix() for path in file_paths),
        "file_size_mb": sum(path.stat().st_size for path in file_paths) / 1_000_000,
        "rows_read": representative["rows_read"],
        "result_rows": representative["result_rows"],
        "ridership_total": representative["ridership_total"],
        "median_elapsed_seconds": median(elapsed_values),
        "min_elapsed_seconds": min(elapsed_values),
        "max_elapsed_seconds": max(elapsed_values),
        "status": "ok",
    }


def _missing_file_row(workflow, scenario, file_paths, runs):
    """Return a benchmark row for workflows skipped because input is missing."""
    return {
        "query": scenario.query,
        "workflow": workflow,
        "period": scenario.period,
        "station_codes": ", ".join(scenario.station_codes or ()),
        "runs": int(runs),
        "file_path": "; ".join(path.relative_to(PROJECT_ROOT).as_posix() for path in file_paths),
        "file_size_mb": pd.NA,
        "rows_read": pd.NA,
        "result_rows": pd.NA,
        "ridership_total": pd.NA,
        "median_elapsed_seconds": pd.NA,
        "min_elapsed_seconds": pd.NA,
        "max_elapsed_seconds": pd.NA,
        "status": "missing_input",
    }


def _raw_files_for_scenario(scenario):
    """Return raw hourly OD files needed by a scenario."""
    return tuple(hourly_od_path_for_year(year) for year in scenario.years)


def _available_hourly_years():
    """Return sorted years available as raw hourly OD files."""
    years = []
    for path in discover_hourly_od_files():
        try:
            years.append(int(path.stem.rsplit("-", 1)[-1]))
        except ValueError:
            continue
    return tuple(sorted(years))


def _warm_processed_readers():
    """Warm pandas CSV/Parquet readers so query timings do not include import setup."""
    if PROCESSED_STATION_RIDERSHIP_MONTHLY_SUMMARY_CSV.exists():
        pd.read_csv(PROCESSED_STATION_RIDERSHIP_MONTHLY_SUMMARY_CSV, nrows=1)
    if PROCESSED_STATION_RIDERSHIP_MONTHLY_SUMMARY_PARQUET.exists():
        pd.read_parquet(PROCESSED_STATION_RIDERSHIP_MONTHLY_SUMMARY_PARQUET).head(1)
# endregion


# region CLI
def parse_args():
    """Parse command-line arguments."""
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--year", type=int, default=DEFAULT_YEAR)
    parser.add_argument("--month", type=int, default=DEFAULT_MONTH)
    parser.add_argument("--runs", type=int, default=3)
    parser.add_argument("--suite", choices=("single", "representative"), default="single")
    parser.add_argument("--no-warmup", action="store_true")
    parser.add_argument("--output", type=Path, default=BENCHMARK_OUTPUT_CSV)
    return parser.parse_args()


def main():
    """Run benchmarks and write the result CSV."""
    args = parse_args()
    if args.runs < 1:
        raise ValueError("--runs must be at least 1.")

    benchmark_df = run_benchmarks(
        args.year,
        args.month,
        args.runs,
        args.suite,
        warmup=not args.no_warmup,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    benchmark_df.to_csv(args.output, index=False)

    display_columns = [
        "query",
        "workflow",
        "period",
        "station_codes",
        "file_size_mb",
        "rows_read",
        "result_rows",
        "median_elapsed_seconds",
        "speedup_vs_raw",
        "status",
    ]
    print(benchmark_df[display_columns].to_string(index=False))
    print(f"\nBenchmark written to {args.output}")


if __name__ == "__main__":
    main()
# endregion
