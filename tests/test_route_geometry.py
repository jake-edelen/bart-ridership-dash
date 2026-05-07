import sys
import unittest
from pathlib import Path

from shapely.geometry import Point

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.data_loader import load_raw_routes, load_stations
from src.route_builder import (
    build_current_routes,
    validate_extension_station_geometry,
    validate_station_geometry_coverage,
)
from src.service_patterns import get_active_station_codes_for_range
from src.station_mapping import STATION_MAPPING


class GeometryCoverageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.stations_gdf = load_stations()
        cls.raw_routes_gdf = load_raw_routes()

    def test_service_pattern_station_geometry_coverage(self):
        codes = get_active_station_codes_for_range("2018-01-01", "2025-12-31")
        coverage = validate_station_geometry_coverage(
            self.stations_gdf,
            STATION_MAPPING,
            codes,
        )
        self.assertEqual((), coverage["missing_mapping"])
        self.assertEqual((), coverage["missing_geometry"])
        self.assertEqual((), coverage["null_geometry"])

    def test_extension_stations_touch_raw_route_geometry(self):
        problems = validate_extension_station_geometry(
            self.raw_routes_gdf,
            self.stations_gdf,
            STATION_MAPPING,
            ("PC", "AN", "ML", "BE", "WD", "OA"),
        )
        self.assertEqual((), problems)

    def test_processed_station_and_route_artifacts_drive_geometry_loading(self):
        self.assertFalse(self.stations_gdf.empty)
        self.assertFalse(self.raw_routes_gdf.empty)
        self.assertTrue({"Name", "lat", "lon", "geometry"} <= set(self.stations_gdf.columns))
        self.assertTrue(self.stations_gdf.geometry.notna().all())
        self.assertTrue(self.raw_routes_gdf.geometry.notna().all())
        self.assertIn("Oakland International Airport", set(self.stations_gdf["Name"]))
        self.assertNotIn("Coliseum/Airport Connector", set(self.stations_gdf["Name"]))

    def test_route_builder_outputs_current_routes(self):
        routes = build_current_routes(self.raw_routes_gdf, self.stations_gdf, service_date="2025-01-01")
        self.assertGreaterEqual(len(routes), 6)
        self.assertTrue({"route", "Color", "service_id", "geometry"} <= set(routes.columns))
        self.assertFalse(routes.geometry.isna().any())

    def test_route_builder_outputs_airport_connector_route(self):
        routes = build_current_routes(self.raw_routes_gdf, self.stations_gdf, service_date="2025-01-01")
        airport_route = routes[routes["service_id"] == "oak_airport_connector_all"]

        self.assertEqual(1, len(airport_route))
        self.assertEqual(("CL", "OA"), airport_route.iloc[0]["station_codes"])
        self.assertEqual(1, len(airport_route.iloc[0].geometry.geoms))

    def test_route_builder_keeps_berryessa_routes_continuous_through_oakland(self):
        routes = build_current_routes(self.raw_routes_gdf, self.stations_gdf, service_date="2025-01-01")
        target_routes = routes[
            routes["route"].isin(
                [
                    "Berryessa/North San José to Richmond",
                    "Berryessa/North San José to Daly City",
                ]
            )
        ]

        self.assertEqual(2, len(target_routes))
        for _, route in target_routes.iterrows():
            route_parts = list(route.geometry.geoms)
            gaps = [
                Point(first.coords[-1]).distance(Point(second.coords[0]))
                for first, second in zip(route_parts, route_parts[1:])
            ]
            self.assertLess(max(gaps), 1e-3)

    def test_route_builder_representative_eras_align_to_station_sequences(self):
        for service_date in ("2018-04-01", "2018-06-01", "2020-07-01", "2025-01-01"):
            routes = build_current_routes(self.raw_routes_gdf, self.stations_gdf, service_date=service_date)
            self.assertFalse(routes.empty)
            self.assertFalse(routes.geometry.isna().any())

    def test_route_builder_routes_pass_service_pattern_stations(self):
        routes = build_current_routes(self.raw_routes_gdf, self.stations_gdf, service_date="2025-01-01")
        station_points = {
            row["Name"]: row.geometry
            for _, row in self.stations_gdf.iterrows()
        }

        for _, route in routes.iterrows():
            for station_code in route["station_codes"]:
                station_name = STATION_MAPPING[station_code]
                self.assertLess(
                    route.geometry.distance(station_points[station_name]),
                    1e-3,
                    msg=f"{route['route']} does not pass near {station_name}.",
                )


if __name__ == "__main__":
    unittest.main()
