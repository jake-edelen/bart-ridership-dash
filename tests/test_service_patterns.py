import unittest
import sys
from pathlib import Path

from shapely.geometry import Point

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.data_loader import (
    build_monthly_station_ridership_summaries,
    discover_ridership_workbooks,
    extract_ridership_station_codes,
    find_ridership_workbook,
    get_available_ridership_periods,
    get_ridership_station_codes_by_year,
    load_station_ridership_for_period,
    load_monthly_ridership_long,
    load_raw_routes,
    load_stations,
    parse_ridership_period,
    summarize_station_ridership,
)
from src.config import RIDERSHIP_2018_DIR
from src.route_builder import (
    build_current_routes,
    validate_extension_station_geometry,
    validate_station_geometry_coverage,
)
from src.service_patterns import (
    get_active_station_codes_for_range,
    get_service_patterns_for_date,
    get_service_patterns_for_range,
    validate_service_pattern_station_codes,
)
from src.station_mapping import STATION_MAPPING


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


class RidershipWorkbookTests(unittest.TestCase):
    def test_ridership_2022_is_discovered(self):
        workbooks = discover_ridership_workbooks()
        files_2022 = [path for path in workbooks if "ridership_2022" in str(path)]
        self.assertEqual(12, len(files_2022))

    def test_old_workbook_format_headers(self):
        codes = set(extract_ridership_station_codes("data/raw/ridership_2022/Ridership_202201.xlsx"))
        self.assertTrue({"RM", "EN", "EP", "PC", "AN", "ML", "BE"} <= codes)

    def test_new_workbook_format_headers(self):
        codes = set(extract_ridership_station_codes("data/raw/ridership_2024/Ridership_202401.xlsx"))
        self.assertTrue({"RM", "EN", "EP", "PC", "AN", "ML", "BE", "12th", "19th"} <= codes)

    def test_total_sheet_order_does_not_matter(self):
        codes = set(extract_ridership_station_codes("data/raw/ridership_OD_2025/Ridership_202505.xlsx"))
        self.assertTrue({"RM", "EN", "EP", "PC", "AN", "ML", "BE"} <= codes)

    def test_codes_by_year_includes_2022(self):
        codes_by_year = get_ridership_station_codes_by_year()
        self.assertIn("2022", codes_by_year)
        self.assertTrue({"PC", "AN", "ML", "BE"} <= set(codes_by_year["2022"]))

    def test_available_ridership_periods_include_project_range(self):
        periods = set(get_available_ridership_periods())

        self.assertIn((2018, 1), periods)
        self.assertIn((2025, 5), periods)

    def test_ridership_period_parser_and_finder(self):
        workbook = find_ridership_workbook(2025, 5)

        self.assertEqual((2025, 5), parse_ridership_period(workbook))

    def test_load_station_ridership_for_selected_period(self):
        station_summary = load_station_ridership_for_period(2025, 5)

        self.assertFalse(station_summary.empty)
        self.assertTrue({"Entry Station", "Ridership", "Full Station Name"} <= set(station_summary.columns))
        self.assertTrue({"ML", "BE"} <= set(station_summary["Entry Station"]))

    def test_ridership_parser_excludes_total_rows_and_columns(self):
        station_summary = summarize_station_ridership(
            load_monthly_ridership_long("data/raw/ridership_OD_2025/Ridership_202505.xlsx")
        )
        station_codes = set(station_summary["Entry Station"].astype(str))

        self.assertNotIn("nan", station_codes)
        self.assertNotIn("Grand Total", station_codes)
        self.assertAlmostEqual(4670575, float(station_summary["Ridership"].sum()))

    def test_monthly_summary_builder_adds_period_columns(self):
        station_summary = build_monthly_station_ridership_summaries(
            raw_root=RIDERSHIP_2018_DIR.parent / "ridership_2018"
        )

        self.assertFalse(station_summary.empty)
        self.assertEqual({2018}, set(station_summary["Year"]))
        self.assertTrue({1, 12} <= set(station_summary["Month"]))


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


class DashRouteYearTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import app as dash_app

        cls.dash_app = dash_app

    def _route_names_for_year(self, year):
        return set(self.dash_app._routes_for_year(year)["route"].tolist())

    def test_route_years_expose_different_service_patterns(self):
        routes_2018 = self._route_names_for_year(2018)
        routes_2025 = self._route_names_for_year(2025)

        self.assertIn("Fremont to Richmond", routes_2018)
        self.assertNotIn("Berryessa/North San José to Richmond", routes_2018)
        self.assertIn("Berryessa/North San José to Richmond", routes_2025)
        self.assertIn("Antioch to SFO International Airport", routes_2025)
        self.assertIn("Coliseum to Oakland Airport", routes_2025)
        self.assertIn("Richmond to Millbrae + SFO", routes_2025)

    def test_route_dropdown_resets_invalid_route_for_selected_year(self):
        options, value = self.dash_app.update_route_dropdown(
            2018,
            "Berryessa/North San José to Richmond",
        )

        self.assertEqual("all", value)
        self.assertIn({"label": "All Routes", "value": "all"}, options)
        self.assertNotIn(
            "Berryessa/North San José to Richmond",
            {option["value"] for option in options},
        )

    def test_route_dropdown_keeps_valid_route_for_selected_year(self):
        _, value = self.dash_app.update_route_dropdown(
            2025,
            "Berryessa/North San José to Richmond",
        )

        self.assertEqual("Berryessa/North San José to Richmond", value)

    def test_route_dropdown_includes_airport_connector(self):
        for year in (2018, 2025):
            options, value = self.dash_app.update_route_dropdown(year, "Coliseum to Oakland Airport")
            option_values = {option["value"] for option in options}

            self.assertIn("Coliseum to Oakland Airport", option_values)
            self.assertEqual("Coliseum to Oakland Airport", value)

    def test_ridership_month_dropdown_matches_selected_year(self):
        options, value = self.dash_app.update_ridership_month_dropdown(2025, 13)
        option_values = {option["value"] for option in options}

        self.assertIn(5, option_values)
        self.assertNotIn(13, option_values)
        self.assertEqual(1, value)

    def test_map_callback_returns_figures_for_selected_year(self):
        figures = self.dash_app.update_maps(2025, "all", 2018, 1, None, None)

        self.assertEqual(3, len(figures))
        self.assertTrue(all(figure.__class__.__name__ == "Figure" for figure in figures))

    def test_map_callback_uses_selected_ridership_period(self):
        _, ridership_map, bar_chart = self.dash_app.update_maps(2025, "all", 2025, 5, None, None)

        self.assertEqual("Station Ridership Map (May 2025)", ridership_map.layout.title.text)
        self.assertEqual("Top 10 Stations by May 2025 Ridership", bar_chart.layout.title.text)

    def test_ridership_map_marker_sizes_are_capped(self):
        routes_gdf = self.dash_app._routes_for_year(2025)
        stations_gdf = self.dash_app._stations_for_ridership_period(2025, 5)
        figure = self.dash_app._build_ridership_route_map(
            routes_gdf,
            "all",
            stations_gdf,
            "May 2025",
        )
        station_trace = figure.data[-1]

        self.assertLessEqual(max(station_trace.marker.size), 34)
        self.assertGreaterEqual(min(station_trace.marker.size), 5)

    def test_route_legend_has_one_clickable_item_per_route(self):
        routes_gdf = self.dash_app._routes_for_year(2025)
        figure = self.dash_app._build_colored_route_map(routes_gdf)
        visible_legend_names = [
            trace.name
            for trace in figure.data
            if trace.showlegend is True
        ]

        self.assertEqual(len(set(routes_gdf["route"])) + 1, len(visible_legend_names))
        self.assertEqual(len(set(visible_legend_names)), len(visible_legend_names))
        self.assertIn("BART Stations", visible_legend_names)
        self.assertEqual("toggle", figure.layout.legend.itemclick)
        self.assertEqual("toggleothers", figure.layout.legend.itemdoubleclick)

    def test_colored_route_map_offsets_route_geometries_only_for_display(self):
        routes_gdf = self.dash_app._routes_for_year(2025)
        visible_routes_gdf = self.dash_app._offset_colored_route_geometries(routes_gdf)

        for _, route in routes_gdf.iterrows():
            visible_route = visible_routes_gdf[
                visible_routes_gdf["service_id"] == route["service_id"]
            ].iloc[0]
            offset_meters = self.dash_app.COLOR_ROUTE_OFFSETS_METERS[route["Color"]]

            if offset_meters == 0:
                self.assertTrue(route.geometry.equals_exact(visible_route.geometry, 0))
            else:
                self.assertFalse(route.geometry.equals_exact(visible_route.geometry, 0))

        original_routes_gdf = self.dash_app._routes_for_year(2025)
        self.assertTrue(routes_gdf.geometry.geom_equals_exact(original_routes_gdf.geometry, 0).all())

    def test_colored_route_map_uses_softened_display_colors(self):
        routes_gdf = self.dash_app._routes_for_year(2025)
        figure = self.dash_app._build_colored_route_map(routes_gdf)
        route_traces = [
            trace
            for trace in figure.data
            if trace.name != "BART Stations"
        ]

        for trace in route_traces:
            route_color = routes_gdf[routes_gdf["route"] == trace.name].iloc[0]["Color"]
            self.assertEqual(
                self.dash_app.DISPLAY_ROUTE_COLORS[route_color],
                trace.line.color,
            )

    def test_black_route_map_uses_single_ungrouped_trace_per_route(self):
        routes_gdf = self.dash_app._routes_for_year(2025)
        figure = self.dash_app._build_ridership_route_map(
            routes_gdf,
            "all",
            self.dash_app.stations_gdf,
        )
        route_traces = [
            trace
            for trace in figure.data
            if trace.name != "BART Stations"
        ]

        self.assertEqual(len(set(routes_gdf["route"])), len(route_traces))
        self.assertTrue(all(trace.showlegend is True for trace in route_traces))
        self.assertTrue(all(trace.legendgroup is None for trace in route_traces))
        self.assertTrue(all(None in trace.lat for trace in route_traces))
        self.assertEqual(1.01, figure.layout.legend.x)
        self.assertEqual(0.96, figure.layout.legend.y)
        self.assertEqual("left", figure.layout.legend.xanchor)

    def test_ridership_bar_chart_defaults_to_top_10_full_station_names(self):
        figure = self.dash_app._build_ridership_bar_chart(self.dash_app.station_ridership)

        self.assertEqual("Top 10 Stations by January 2018 Ridership", figure.layout.title.text)
        self.assertEqual(10, len(figure.data[0].y))
        self.assertTrue(all(name not in set(self.dash_app.station_mapping) for name in figure.data[0].y))

    def test_ridership_bar_chart_compares_selected_full_station_names(self):
        filtered_ridership = self.dash_app.station_ridership[
            self.dash_app.station_ridership["Full Station Name"].isin(
                ["Embarcadero", "Montgomery St"]
            )
        ]
        figure = self.dash_app._build_ridership_bar_chart(filtered_ridership)

        self.assertEqual("Selected Station Ridership (January 2018)", figure.layout.title.text)
        self.assertEqual({"Embarcadero", "Montgomery St"}, set(figure.data[0].y))

    def test_ridership_source_remains_january_2018(self):
        expected = summarize_station_ridership(
            load_monthly_ridership_long(RIDERSHIP_2018_DIR / "Ridership_201801.xlsx")
        )
        actual = self.dash_app.app_data.station_ridership

        self.assertFalse(actual.empty)
        self.assertAlmostEqual(float(expected["Ridership"].sum()), float(actual["Ridership"].sum()))


if __name__ == "__main__":
    unittest.main()
