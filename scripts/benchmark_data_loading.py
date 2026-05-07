"""Benchmark raw CSV aggregation against processed CSV and Parquet summaries.

This script produces local, machine-specific timing evidence for the data
engineering payoff of pre-aggregated Parquet summaries. It does not run during
app startup.
"""

# region Imports And Path Setup
from argparse import ArgumentParser
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
from src.data_loader import HOURLY_OD_COLUMNS, hourly_od_path_for_year  # noqa: E402
from src.station_mapping import STATION_MAPPING, normalize_hourly_od_station_code  # noqa: E402
# endregion


# region Constants
BENCHMARK_OUTPUT_CSV = PROCESSED / "benchmark_data_loading.csv"
DEFAULT_YEAR = 2025
DEFAULT_MONTH = 5
# endregion


# region Benchmark Operations
def benchmark_raw_hourly_od(year, month, raw_file):
    """Read a raw hourly OD year CSV, filter one month, and aggregate stations."""
    period_prefix = f"{int(year)}-{int(month):02d}"
    rows_read = 0
    monthly_chunks = []
    start = perf_counter()

    for chunk in pd.read_csv(
        raw_file,
        header=None,
        names=HOURLY_OD_COLUMNS,
        usecols=["Date", "ORIGIN", "Number of Exits"],
        chunksize=750_000,
    ):
        rows_read += len(chunk)
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
        result = (
            pd.concat(monthly_chunks, ignore_index=True)
            .groupby("Entry Station", as_index=False)["Ridership"]
            .sum()
        )
    else:
        result = pd.DataFrame(columns=["Entry Station", "Ridership"])

    elapsed_seconds = perf_counter() - start
    return {
        "rows_read": rows_read,
        "result_rows": len(result),
        "ridership_total": float(result["Ridership"].sum()) if not result.empty else 0.0,
        "elapsed_seconds": elapsed_seconds,
    }


def benchmark_processed_csv(year, month, csv_file):
    """Read the processed CSV summary and filter one period."""
    start = perf_counter()
    summary = pd.read_csv(csv_file)
    result = summary[
        (summary["Year"].astype(int) == int(year))
        & (summary["Month"].astype(int) == int(month))
    ].copy()
    elapsed_seconds = perf_counter() - start
    return {
        "rows_read": len(summary),
        "result_rows": len(result),
        "ridership_total": float(result["Ridership"].sum()) if not result.empty else 0.0,
        "elapsed_seconds": elapsed_seconds,
    }


def benchmark_processed_parquet(year, month, parquet_file):
    """Read the processed Parquet summary and filter one period."""
    start = perf_counter()
    summary = pd.read_parquet(parquet_file)
    result = summary[
        (summary["Year"].astype(int) == int(year))
        & (summary["Month"].astype(int) == int(month))
    ].copy()
    elapsed_seconds = perf_counter() - start
    return {
        "rows_read": len(summary),
        "result_rows": len(result),
        "ridership_total": float(result["Ridership"].sum()) if not result.empty else 0.0,
        "elapsed_seconds": elapsed_seconds,
    }
# endregion


# region Result Assembly
def run_benchmarks(year, month, runs):
    """Run all benchmark workflows and return one summary row per workflow."""
    raw_file = hourly_od_path_for_year(year)
    workflows = [
        (
            "raw_hourly_od",
            raw_file,
            lambda: benchmark_raw_hourly_od(year, month, raw_file),
        ),
        (
            "processed_csv_summary",
            PROCESSED_STATION_RIDERSHIP_MONTHLY_SUMMARY_CSV,
            lambda: benchmark_processed_csv(
                year,
                month,
                PROCESSED_STATION_RIDERSHIP_MONTHLY_SUMMARY_CSV,
            ),
        ),
        (
            "processed_parquet_summary",
            PROCESSED_STATION_RIDERSHIP_MONTHLY_SUMMARY_PARQUET,
            lambda: benchmark_processed_parquet(
                year,
                month,
                PROCESSED_STATION_RIDERSHIP_MONTHLY_SUMMARY_PARQUET,
            ),
        ),
    ]

    rows = []
    for workflow, file_path, operation in workflows:
        if not file_path.exists():
            rows.append(_missing_file_row(workflow, year, month, file_path, runs))
            continue

        run_results = [operation() for _ in range(runs)]
        elapsed_values = [result["elapsed_seconds"] for result in run_results]
        representative = run_results[elapsed_values.index(median(elapsed_values))]
        rows.append(
            {
                "workflow": workflow,
                "period": f"{int(year)}-{int(month):02d}",
                "runs": int(runs),
                "file_path": file_path.relative_to(PROJECT_ROOT).as_posix(),
                "file_size_mb": file_path.stat().st_size / 1_000_000,
                "rows_read": representative["rows_read"],
                "result_rows": representative["result_rows"],
                "ridership_total": representative["ridership_total"],
                "median_elapsed_seconds": median(elapsed_values),
                "min_elapsed_seconds": min(elapsed_values),
                "max_elapsed_seconds": max(elapsed_values),
                "status": "ok",
            }
        )

    benchmark_df = pd.DataFrame(rows)
    raw_elapsed = benchmark_df.loc[
        benchmark_df["workflow"] == "raw_hourly_od",
        "median_elapsed_seconds",
    ]
    if not raw_elapsed.empty and pd.notna(raw_elapsed.iloc[0]) and raw_elapsed.iloc[0] > 0:
        benchmark_df["speedup_vs_raw"] = raw_elapsed.iloc[0] / benchmark_df["median_elapsed_seconds"]
    else:
        benchmark_df["speedup_vs_raw"] = pd.NA

    return benchmark_df


def _missing_file_row(workflow, year, month, file_path, runs):
    """Return a benchmark row for workflows skipped because input is missing."""
    return {
        "workflow": workflow,
        "period": f"{int(year)}-{int(month):02d}",
        "runs": int(runs),
        "file_path": file_path.relative_to(PROJECT_ROOT).as_posix(),
        "file_size_mb": pd.NA,
        "rows_read": pd.NA,
        "result_rows": pd.NA,
        "ridership_total": pd.NA,
        "median_elapsed_seconds": pd.NA,
        "min_elapsed_seconds": pd.NA,
        "max_elapsed_seconds": pd.NA,
        "status": "missing_input",
    }
# endregion


# region CLI
def parse_args():
    """Parse command-line arguments."""
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--year", type=int, default=DEFAULT_YEAR)
    parser.add_argument("--month", type=int, default=DEFAULT_MONTH)
    parser.add_argument("--runs", type=int, default=3)
    parser.add_argument("--output", type=Path, default=BENCHMARK_OUTPUT_CSV)
    return parser.parse_args()


def main():
    """Run benchmarks and write the result CSV."""
    args = parse_args()
    if args.runs < 1:
        raise ValueError("--runs must be at least 1.")

    benchmark_df = run_benchmarks(args.year, args.month, args.runs)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    benchmark_df.to_csv(args.output, index=False)

    display_columns = [
        "workflow",
        "period",
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
