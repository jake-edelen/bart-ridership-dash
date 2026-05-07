import sys
import unittest
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.config import (
    HOURLY_OD_2018_CSV,
    PROCESSED_HOURLY_COMPLETENESS_CSV,
    PROCESSED_HOURLY_COMPLETENESS_PARQUET,
    PROCESSED_HOURLY_STATION_MONTHLY_SUMMARY,
    PROCESSED_HOURLY_STATION_MONTHLY_SUMMARY_PARQUET,
    PROCESSED_HOURLY_VALIDATION_CSV,
    PROCESSED_HOURLY_VALIDATION_PARQUET,
    PROCESSED_STATION_RIDERSHIP_MONTHLY_SUMMARY_PARQUET,
    RIDERSHIP_2018_DIR,
)
from src.data_loader import (
    build_monthly_station_ridership_summaries,
    discover_ridership_workbooks,
    extract_ridership_station_codes,
    find_ridership_workbook,
    get_available_ridership_periods,
    get_ridership_station_codes_by_year,
    load_hourly_station_ridership_for_month,
    load_monthly_ridership_long,
    load_monthly_station_ridership_summary,
    load_station_ridership_for_period,
    parse_ridership_period,
    summarize_station_ridership,
)
from src.station_mapping import normalize_workbook_station_code


class RidershipWorkbookTests(unittest.TestCase):
    def test_ridership_2022_is_discovered(self):
        workbooks = discover_ridership_workbooks()
        files_2022 = [path for path in workbooks if "ridership_2022" in str(path)]
        self.assertEqual(12, len(files_2022))

    def test_old_workbook_format_headers(self):
        codes = set(extract_ridership_station_codes("data/raw/ridership_2022/Ridership_202201.xlsx"))
        self.assertTrue({"RM", "EN", "EP", "PITT", "WP", "WD", "PC", "AN", "ML", "BE"} <= codes)
        self.assertNotIn("WS", codes)

    def test_new_workbook_format_headers(self):
        codes = set(extract_ridership_station_codes("data/raw/ridership_2024/Ridership_202401.xlsx"))
        self.assertTrue({"RM", "EN", "EP", "PITT", "WP", "WD", "PC", "AN", "ML", "BE", "12th", "19th"} <= codes)
        self.assertNotIn("WS", codes)

    def test_total_sheet_order_does_not_matter(self):
        codes = set(extract_ridership_station_codes("data/raw/ridership_OD_2025/Ridership_202505.xlsx"))
        self.assertTrue({"RM", "EN", "EP", "PITT", "WP", "WD", "PC", "AN", "ML", "BE"} <= codes)
        self.assertNotIn("WS", codes)

    def test_workbook_station_aliases_use_canonical_app_codes(self):
        self.assertEqual("PITT", normalize_workbook_station_code("WP"))
        self.assertEqual("WP", normalize_workbook_station_code("WD"))
        self.assertEqual("WD", normalize_workbook_station_code("WS"))
        self.assertEqual("OA", normalize_workbook_station_code("OA"))

    def test_codes_by_year_includes_2022(self):
        codes_by_year = get_ridership_station_codes_by_year()
        self.assertIn("2022", codes_by_year)
        self.assertTrue({"PITT", "WP", "WD", "PC", "AN", "ML", "BE"} <= set(codes_by_year["2022"]))
        self.assertNotIn("WS", set(codes_by_year["2022"]))

    def test_available_ridership_periods_include_project_range(self):
        periods = set(get_available_ridership_periods())

        self.assertIn((2018, 1), periods)
        self.assertIn((2025, 5), periods)

    def test_processed_ridership_summary_drives_period_loading(self):
        monthly_summary = load_monthly_station_ridership_summary()
        periods = get_available_ridership_periods()
        station_summary = load_station_ridership_for_period(2018, 1)

        self.assertIn((2018, 1), periods)
        self.assertIn("Source Type", monthly_summary.columns)
        self.assertEqual(9785965, int(station_summary["Ridership"].sum()))

    def test_processed_hourly_monthly_summary_artifact_exists(self):
        hourly_summary = pd.read_csv(PROCESSED_HOURLY_STATION_MONTHLY_SUMMARY)

        self.assertTrue(
            {
                "Year",
                "Month",
                "Entry Station",
                "station_id",
                "Ridership",
                "Full Station Name",
                "display_name",
                "Period",
                "Source Type",
                "Source File",
            }
            <= set(hourly_summary.columns)
        )
        self.assertIn(
            (2018, 1),
            {
                (int(row["Year"]), int(row["Month"]))
                for _, row in hourly_summary[["Year", "Month"]].drop_duplicates().iterrows()
            },
        )
        january_2018 = hourly_summary[
            (hourly_summary["Year"].astype(int) == 2018)
            & (hourly_summary["Month"].astype(int) == 1)
        ]
        self.assertEqual(9785965, int(january_2018["Ridership"].sum()))

    def test_processed_parquet_artifacts_exist_and_match_csv_rows(self):
        artifact_pairs = (
            (
                PROCESSED_HOURLY_STATION_MONTHLY_SUMMARY,
                PROCESSED_HOURLY_STATION_MONTHLY_SUMMARY_PARQUET,
            ),
            (PROCESSED_HOURLY_VALIDATION_CSV, PROCESSED_HOURLY_VALIDATION_PARQUET),
            (PROCESSED_HOURLY_COMPLETENESS_CSV, PROCESSED_HOURLY_COMPLETENESS_PARQUET),
        )

        for csv_path, parquet_path in artifact_pairs:
            csv_data = pd.read_csv(csv_path)
            parquet_data = pd.read_parquet(parquet_path)

            self.assertTrue(parquet_path.exists())
            self.assertEqual(len(csv_data), len(parquet_data))

        self.assertTrue(PROCESSED_STATION_RIDERSHIP_MONTHLY_SUMMARY_PARQUET.exists())

    def test_hourly_workbook_validation_artifact_flags_station_differences(self):
        validation = pd.read_csv(PROCESSED_HOURLY_VALIDATION_CSV)

        self.assertTrue(
            {
                "Year",
                "Month",
                "Period",
                "Entry Station",
                "station_id",
                "display_name",
                "Workbook Ridership",
                "Hourly OD Ridership",
                "Difference",
                "Has Difference",
            }
            <= set(validation.columns)
        )
        self.assertTrue(validation["Has Difference"].isin([True, False]).all())
        self.assertTrue(
            ((validation["Year"].astype(int) == 2018) & (validation["Month"].astype(int) == 2)).any()
        )
        self.assertTrue(validation["Has Difference"].any())

    def test_hourly_completeness_artifact_documents_quality_flags(self):
        quality = pd.read_csv(PROCESSED_HOURLY_COMPLETENESS_CSV)

        self.assertTrue(
            {
                "Period",
                "Missing Full Days",
                "Missing Full Day List",
                "Raw Hourly Ridership",
                "Imputed Hourly Ridership",
                "Imputed Ridership Added",
                "Quality Flag",
            }
            <= set(quality.columns)
        )
        self.assertTrue(
            (quality["Quality Flag"] == "missing_full_days_imputed").any()
        )
        self.assertGreater(
            float(
                quality.loc[
                    quality["Quality Flag"] == "missing_full_days_imputed",
                    "Imputed Ridership Added",
                ].sum()
            ),
            0,
        )

    def test_data_loading_benchmark_artifact_has_expected_workflows(self):
        benchmark_path = PROJECT_ROOT / "data" / "processed" / "benchmark_data_loading.csv"
        benchmark = pd.read_csv(benchmark_path)

        self.assertTrue(
            {
                "workflow",
                "period",
                "file_size_mb",
                "rows_read",
                "median_elapsed_seconds",
                "speedup_vs_raw",
                "status",
            }
            <= set(benchmark.columns)
        )
        self.assertEqual(
            {
                "raw_hourly_od",
                "processed_csv_summary",
                "processed_parquet_summary",
            },
            set(benchmark["workflow"]),
        )

    def test_ridership_period_parser_and_finder(self):
        workbook = find_ridership_workbook(2025, 5)

        self.assertEqual((2025, 5), parse_ridership_period(workbook))

    def test_load_station_ridership_for_selected_period(self):
        station_summary = load_station_ridership_for_period(2025, 5)

        self.assertFalse(station_summary.empty)
        self.assertTrue(
            {
                "Entry Station",
                "station_id",
                "Ridership",
                "Full Station Name",
                "display_name",
            }
            <= set(station_summary.columns)
        )
        self.assertTrue({"ML", "BE"} <= set(station_summary["Entry Station"]))
        self.assertTrue({"mlpt", "bery"} <= set(station_summary["station_id"]))

    def test_january_2018_hourly_od_fallback_matches_monthly_total(self):
        station_summary = load_hourly_station_ridership_for_month(2018, 1, HOURLY_OD_2018_CSV)

        self.assertAlmostEqual(9785965, float(station_summary["Ridership"].sum()))
        self.assertIn("PITT", set(station_summary["Entry Station"]))

    def test_hourly_od_extension_station_aliases_are_normalized(self):
        station_summary = load_hourly_station_ridership_for_month(2025, 5)

        self.assertTrue({"ML", "BE"} <= set(station_summary["Entry Station"]))
        self.assertGreater(
            float(station_summary.loc[station_summary["Entry Station"] == "ML", "Ridership"].iloc[0]),
            0,
        )
        self.assertGreater(
            float(station_summary.loc[station_summary["Entry Station"] == "BE", "Ridership"].iloc[0]),
            0,
        )

    def test_ridership_parser_excludes_total_rows_and_columns(self):
        station_summary = summarize_station_ridership(
            load_monthly_ridership_long("data/raw/ridership_OD_2025/Ridership_202505.xlsx")
        )
        station_codes = set(station_summary["Entry Station"].astype(str))

        self.assertNotIn("nan", station_codes)
        self.assertNotIn("Grand Total", station_codes)
        self.assertAlmostEqual(4670575, float(station_summary["Ridership"].sum()))

    def test_monthly_summary_builder_adds_period_columns(self):
        hourly_summary = pd.concat(
            [
                load_hourly_station_ridership_for_month(2018, 1, HOURLY_OD_2018_CSV).assign(
                    Year=2018,
                    Month=1,
                    Period="2018-01",
                    **{"Source Type": "hourly_od", "Source File": str(HOURLY_OD_2018_CSV)},
                ),
                load_hourly_station_ridership_for_month(2018, 2, HOURLY_OD_2018_CSV).assign(
                    Year=2018,
                    Month=2,
                    Period="2018-02",
                    **{"Source Type": "hourly_od", "Source File": str(HOURLY_OD_2018_CSV)},
                ),
            ],
            ignore_index=True,
        )
        station_summary = build_monthly_station_ridership_summaries(
            raw_root=RIDERSHIP_2018_DIR,
            hourly_summary=hourly_summary,
        )

        self.assertFalse(station_summary.empty)
        self.assertEqual({2018}, set(station_summary["Year"]))
        self.assertTrue({1, 12} <= set(station_summary["Month"]))
        self.assertIn("Source Type", station_summary.columns)
        january_summary = station_summary[station_summary["Month"] == 1]
        self.assertEqual({"hourly_od"}, set(january_summary["Source Type"]))
        self.assertAlmostEqual(9785965, float(january_summary["Ridership"].sum()))
        february_summary = station_summary[station_summary["Month"] == 2]
        self.assertEqual({"hourly_od"}, set(february_summary["Source Type"]))
        self.assertEqual(
            148391,
            int(february_summary.loc[february_summary["Entry Station"] == "PITT", "Ridership"].iloc[0]),
        )
        self.assertEqual(
            76895,
            int(february_summary.loc[february_summary["Entry Station"] == "WP", "Ridership"].iloc[0]),
        )
        self.assertEqual(
            70757,
            int(february_summary.loc[february_summary["Entry Station"] == "WD", "Ridership"].iloc[0]),
        )
        self.assertNotIn("WS", set(february_summary["Entry Station"]))


if __name__ == "__main__":
    unittest.main()
