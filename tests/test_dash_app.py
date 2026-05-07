import sys
import unittest
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.data_loader import load_station_ridership_for_period


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

    def test_ridership_year_drives_route_dropdown_options(self):
        options, _ = self.dash_app.update_route_dropdown(2018, "all")
        option_values = {option["value"] for option in options}

        self.assertIn("Fremont to Richmond", option_values)
        self.assertNotIn("Berryessa/North San José to Richmond", option_values)

    def test_ridership_month_dropdown_matches_selected_year(self):
        options, value = self.dash_app.update_ridership_month_dropdown(2025, 13)
        option_values = {option["value"] for option in options}

        self.assertIn(5, option_values)
        self.assertNotIn(13, option_values)
        self.assertEqual(1, value)

    def test_station_dropdown_options_match_selected_ridership_period(self):
        january_options = {
            option["value"]
            for option in self.dash_app._station_options_for_ridership_period(2018, 1)
        }
        may_2025_options = {
            option["value"]
            for option in self.dash_app._station_options_for_ridership_period(2025, 5)
        }

        self.assertNotIn("Milpitas", january_options)
        self.assertIn("Milpitas", may_2025_options)

    def test_station_dropdown_clears_stations_missing_from_period(self):
        _, station1, _, station2 = self.dash_app.update_station_dropdowns(
            2018,
            1,
            "Milpitas",
            "Embarcadero",
        )

        self.assertIsNone(station1)
        self.assertEqual("Embarcadero", station2)

    def test_map_callback_returns_figures_for_selected_year(self):
        figures = self.dash_app.update_maps("all", 2018, 1, None, None)

        self.assertEqual(3, len(figures))
        self.assertTrue(all(figure.__class__.__name__ == "Figure" for figure in figures))

    def test_map_callback_uses_selected_ridership_period(self):
        _, ridership_map, bar_chart = self.dash_app.update_maps("all", 2025, 5, None, None)

        self.assertEqual("Station Ridership Map (May 2025)", ridership_map.layout.title.text)
        self.assertEqual("Top 10 Stations by May 2025 Ridership", bar_chart.layout.title.text)

    def test_map_callback_uses_ridership_year_for_route_geometry(self):
        colored_map, _, _ = self.dash_app.update_maps("all", 2018, 1, None, None)
        route_names = {trace.name for trace in colored_map.data if trace.name != "BART Stations"}

        self.assertIn("Fremont to Richmond", route_names)
        self.assertNotIn("Berryessa/North San José to Richmond", route_names)

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

    def test_selected_station_color_scale_uses_period_range(self):
        _, ridership_map, bar_chart = self.dash_app.update_maps(
            "all",
            2018,
            1,
            "19th St/Oakland",
            "16th St/Mission",
        )
        station_trace = ridership_map.data[-1]
        period_max = float(
            self.dash_app._station_ridership_for_period(2018, 1)["Ridership"].max()
        )

        self.assertEqual(0, station_trace.marker.cmin)
        self.assertEqual(period_max, station_trace.marker.cmax)
        self.assertEqual(0, bar_chart.layout.coloraxis.cmin)
        self.assertEqual(period_max, bar_chart.layout.coloraxis.cmax)

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

        self.assertEqual("Selected Station Comparison (January 2018)", figure.layout.title.text)
        self.assertEqual({"Embarcadero", "Montgomery St"}, set(figure.data[0].y))
        self.assertEqual(1, len(figure.layout.annotations))
        self.assertIn("riders higher than", figure.layout.annotations[0].text)

    def test_colored_route_map_title_has_no_explanatory_subheader(self):
        routes_gdf = self.dash_app._routes_for_year(2025)
        figure = self.dash_app._build_colored_route_map(routes_gdf, 2025)

        self.assertEqual("Route Reference Map (2025 service)", figure.layout.title.text)

    def test_ridership_source_remains_january_2018(self):
        expected = load_station_ridership_for_period(2018, 1)
        actual = self.dash_app.app_data.station_ridership

        self.assertFalse(actual.empty)
        self.assertAlmostEqual(float(expected["Ridership"].sum()), float(actual["Ridership"].sum()))

    def test_app_entrypoint_does_not_reference_raw_inputs(self):
        app_source = (PROJECT_ROOT / "app.py").read_text(encoding="utf-8")
        forbidden_runtime_tokens = (
            "data/raw",
            "RAW",
            "BART_STATIONS_DIR",
            "BART_ROUTES_GDB",
            "HOURLY_OD_DIR",
            "HOURLY_OD_2018_CSV",
            "RIDERSHIP_2018_DIR",
            "read_excel",
        )

        for token in forbidden_runtime_tokens:
            self.assertNotIn(token, app_source)


if __name__ == "__main__":
    unittest.main()
