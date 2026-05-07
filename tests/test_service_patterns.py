import sys
import unittest
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.service_patterns import (
    get_active_station_codes_for_range,
    get_service_patterns_for_date,
    get_service_patterns_for_range,
    validate_service_pattern_station_codes,
)


class ServicePatternTests(unittest.TestCase):
    def _active_codes(self, service_date):
        codes = set()
        for pattern in get_service_patterns_for_date(service_date):
            codes.update(pattern.station_codes)
        return codes

    def test_service_date_station_activation(self):
        self.assertFalse({"PC", "AN", "ML", "BE"} & self._active_codes("2018-04-01"))
        self.assertTrue({"PC", "AN"} <= self._active_codes("2018-06-01"))
        self.assertFalse({"ML", "BE"} & self._active_codes("2018-06-01"))
        self.assertTrue({"PC", "AN", "ML", "BE"} <= self._active_codes("2020-07-01"))
        self.assertTrue({"CL", "OA"} <= self._active_codes("2018-04-01"))

    def test_post_2023_patterns_are_returned(self):
        service_ids = {pattern.service_id for pattern in get_service_patterns_for_date("2024-01-01")}
        self.assertIn("red_2023_plus", service_ids)
        self.assertNotIn("red_pre_2023", service_ids)

    def test_red_line_sfo_millbrae_order_changes_in_2023(self):
        pre_2023_red = next(
            pattern
            for pattern in get_service_patterns_for_date("2023-02-15")
            if pattern.service_id == "red_pre_2023"
        )
        post_2023_red = next(
            pattern
            for pattern in get_service_patterns_for_date("2024-01-01")
            if pattern.service_id == "red_2023_plus"
        )

        self.assertLess(
            pre_2023_red.station_codes.index("MB"),
            pre_2023_red.station_codes.index("SO"),
        )
        self.assertLess(
            post_2023_red.station_codes.index("SO"),
            post_2023_red.station_codes.index("MB"),
        )

    def test_range_modes(self):
        end_codes = {
            code
            for pattern in get_service_patterns_for_range("2018-01-01", "2025-12-31", mode="end")
            for code in pattern.station_codes
        }
        start_codes = {
            code
            for pattern in get_service_patterns_for_range("2018-01-01", "2025-12-31", mode="start")
            for code in pattern.station_codes
        }
        all_changes = get_service_patterns_for_range(
            "2018-01-01",
            "2025-12-31",
            mode="all_changes",
        )

        self.assertTrue({"PC", "AN", "ML", "BE"} <= end_codes)
        self.assertFalse({"PC", "AN", "ML", "BE"} & start_codes)
        self.assertGreater(len(all_changes), len(get_service_patterns_for_range("2018-01-01", "2025-12-31")))

    def test_active_station_codes_for_range(self):
        codes = set(get_active_station_codes_for_range("2018-01-01", "2025-12-31"))
        self.assertTrue({"PC", "AN", "ML", "BE"} <= codes)

    def test_service_pattern_codes_are_mapped(self):
        self.assertEqual((), validate_service_pattern_station_codes())


if __name__ == "__main__":
    unittest.main()
