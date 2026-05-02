"""Dash entrypoint for the BART ridership dashboard.

This module should stay focused on presentation: layout, callback wiring, and
Plotly figure construction. Data loading, station mapping, and route geometry
corrections live under `src/` so future data/modeling work can scale without
turning the Dash file back into a preprocessing script.
"""

# region Imports
from functools import lru_cache

import dash
import geopandas as gpd
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash import dcc, html
from dash.dependencies import Input, Output, State
from shapely.geometry import LineString, MultiLineString

from src.data_loader import load_app_data
from src.route_builder import build_current_routes
# endregion


# region Runtime Data
AVAILABLE_ROUTE_YEARS = tuple(range(2018, 2026))
DEFAULT_ROUTE_YEAR = 2025

app_data = load_app_data()
stations_gdf = app_data.stations_gdf
raw_routes_gdf = app_data.raw_routes_gdf
station_ridership = app_data.station_ridership
station_mapping = app_data.station_mapping

year_options = [
    {"label": str(year), "value": year}
    for year in AVAILABLE_ROUTE_YEARS
]


@lru_cache(maxsize=None)
def _routes_for_year(year):
    """Return route geometries active at the end of the selected year."""
    service_date = f"{int(year)}-12-31"
    return build_current_routes(raw_routes_gdf, stations_gdf, service_date=service_date)


def _route_options_for_routes(routes_gdf):
    """Build route dropdown options from the route table currently in scope."""
    options = [
        {"label": route, "value": route}
        for route in routes_gdf["route"].unique()
    ]
    options.insert(0, {"label": "All Routes", "value": "all"})
    return options


route_to_gdf = _routes_for_year(DEFAULT_ROUTE_YEAR)
route_options = _route_options_for_routes(route_to_gdf)
# endregion


# region Figure Styling
ROUTE_LEGEND = dict(
    title=dict(text="Routes (click to toggle)", font=dict(size=12)),
    bgcolor="rgba(255,255,255,0.96)",
    bordercolor="#9AA6B2",
    borderwidth=1,
    font=dict(size=12),
    itemclick="toggle",
    itemdoubleclick="toggleothers",
    groupclick="togglegroup",
    itemsizing="constant",
)

COLOR_ROUTE_OFFSETS_METERS = {
    "Yellow": -200,
    "Red": -100,
    "Blue": 0,
    "Green": 100,
    "Orange": 200,
    "Gray": 0,
}

DISPLAY_ROUTE_COLORS = {
    "Yellow": "Yellow",
    "Red": "Red",
    "Blue": "#007EF5",
    "Green": "Green",
    "Orange": "Orange",
    "Gray": "#9AA6B2",
}

# endregion


# region Layout
def create_layout():
    """Build the Dash layout from already-loaded app data.

    The layout references globally loaded station and route options so Dash can
    render dropdowns immediately. The current layout is intentionally preserved
    from the original app; visual responsiveness is a separate cleanup step.
    """
    return html.Div(
        [
            html.H1(
                "BART System Map & Ridership Data",
                style={"textAlign": "center", "color": "black"},
            ),
            html.Div(
                [
                    html.Div(
                        [
                            html.H4("Dash App Description"),
                            html.P("This Dash App displays BART routes with interactive elements."),
                            html.P(
                                "The dropdown menus control displayed routes and compare "
                                "January 2018 ridership by station."
                            ),
                            html.P(
                                "Ridership intensity is mapped using marker size and color, "
                                "adjusted by selected stations."
                            ),
                        ],
                        style={
                            "padding": "20px",
                            "backgroundColor": "lightgrey",
                            "border": "1px solid grey",
                            "borderRadius": "8px",
                            "alignSelf": "start",
                        },
                    ),
                    html.Div(
                        [
                            dcc.Graph(
                                id="bart-map-black",
                                style={
                                    "height": "520px",
                                    "border": "1px solid lightgrey",
                                },
                            ),
                        ],
                        style={
                            "minWidth": 0,
                        },
                    ),
                    html.Div(
                        [
                            html.Label("Select Route Year:"),
                            dcc.Dropdown(
                                id="route-year-dropdown",
                                options=year_options,
                                value=DEFAULT_ROUTE_YEAR,
                                clearable=False,
                            ),
                            html.Br(),
                            html.Label("Select Route:"),
                            dcc.Dropdown(
                                id="route-dropdown",
                                options=route_options,
                                value="all",
                                clearable=False,
                            ),
                            html.Br(),
                            html.Label("Select Two Stations:"),
                            dcc.Dropdown(
                                id="station-1-dropdown",
                                options=[
                                    {"label": name, "value": name}
                                    for name in stations_gdf["Name"]
                                ],
                                clearable=True,
                            ),
                            html.Br(),
                            dcc.Dropdown(
                                id="station-2-dropdown",
                                options=[
                                    {"label": name, "value": name}
                                    for name in stations_gdf["Name"]
                                ],
                                clearable=True,
                            ),
                        ],
                        style={
                            "padding": "12px",
                            "alignSelf": "start",
                        },
                    ),
                ],
                style={
                    "display": "grid",
                    "gridTemplateColumns": "320px minmax(760px, 1fr) 300px",
                    "gap": "16px",
                    "alignItems": "start",
                    "padding": "0 16px",
                },
            ),
            html.Br(),
            html.H3("Route Map (Colored) and Ridership Comparison"),
            html.Div(
                [
                    dcc.Graph(
                        id="bart-map-colored",
                        style={"width": "50%", "display": "inline-block"},
                    ),
                    dcc.Graph(
                        id="ridership-bar-chart",
                        style={"width": "50%", "display": "inline-block"},
                    ),
                ],
                style={"display": "flex", "justifyContent": "center"},
            ),
        ]
    )


app = dash.Dash(__name__)
app.layout = create_layout
# endregion


# region Callbacks
@app.callback(
    [
        Output("route-dropdown", "options"),
        Output("route-dropdown", "value"),
    ],
    [Input("route-year-dropdown", "value")],
    [State("route-dropdown", "value")],
)
def update_route_dropdown(selected_year, selected_route):
    """Refresh route choices so they match the selected service year."""
    routes_gdf = _routes_for_year(selected_year)
    options = _route_options_for_routes(routes_gdf)
    valid_routes = {option["value"] for option in options}
    value = selected_route if selected_route in valid_routes else "all"
    return options, value


@app.callback(
    [
        Output("bart-map-colored", "figure"),
        Output("bart-map-black", "figure"),
        Output("ridership-bar-chart", "figure"),
    ],
    [
        Input("route-year-dropdown", "value"),
        Input("route-dropdown", "value"),
        Input("station-1-dropdown", "value"),
        Input("station-2-dropdown", "value"),
    ],
)
def update_maps(selected_year, selected_route, station1, station2):
    """Update route maps and ridership chart from selected UI filters.

    Args:
        selected_year: Service-pattern year used for route geometry display.
        selected_route: Route name from the route dropdown, or `"all"`.
        station1: Optional full station name from the first station dropdown.
        station2: Optional full station name from the second station dropdown.

    Returns:
        A tuple of Plotly figures for the colored route map, ridership map, and
        ridership bar chart.
    """
    reverse_mapping = {v: k for k, v in station_mapping.items()}

    if station1 and station2:
        station1_abbr = reverse_mapping.get(station1)
        station2_abbr = reverse_mapping.get(station2)
        filtered_stations = stations_gdf[stations_gdf["Name"].isin([station1, station2])]
        filtered_ridership = station_ridership[
            station_ridership["Entry Station"].isin([station1_abbr, station2_abbr])
        ]
    else:
        filtered_stations = stations_gdf
        filtered_ridership = station_ridership

    selected_routes_gdf = _routes_for_year(selected_year)
    fig_colored = _build_colored_route_map(selected_routes_gdf)
    fig_black = _build_ridership_route_map(
        selected_routes_gdf,
        selected_route,
        filtered_stations,
    )
    bar_fig = _build_ridership_bar_chart(filtered_ridership)

    return fig_colored, fig_black, bar_fig
# endregion


# region Figure Builders
def _add_route_line_traces(fig, routes_gdf, line_width, line_color=None):
    """Add one route trace per route so legend toggles remain reliable."""
    for _, row in routes_gdf.iterrows():
        route_name = row["route"]
        lat = []
        lon = []

        for line in row["geometry"].geoms:
            lat.extend(point[1] for point in line.coords)
            lon.extend(point[0] for point in line.coords)
            lat.append(None)
            lon.append(None)

        fig.add_trace(
            go.Scattermap(
                lat=lat,
                lon=lon,
                mode="lines",
                line=dict(
                    width=line_width,
                    color=line_color or DISPLAY_ROUTE_COLORS.get(row["Color"], row["Color"]),
                ),
                hovertemplate=f"{route_name}<extra></extra>",
                name=route_name,
                showlegend=True,
            )
        )                                                                             

def _build_colored_route_map(routes_gdf):
    """Create the bottom map with route-offset lines in route colors."""
    fig_colored = go.Figure()
    visible_routes_gdf = _offset_colored_route_geometries(routes_gdf)
    _add_route_line_traces(fig_colored, visible_routes_gdf, line_width=4)

    fig_colored.add_trace(
        go.Scattermap(
            lat=stations_gdf["lat"],
            lon=stations_gdf["lon"],
            mode="markers+text",
            marker=dict(size=6, color="black"),
            text=stations_gdf["Name"],
            textposition="top right",
            hoverinfo="text",
            name="BART Stations",
            showlegend=True, 
        )
    )

    fig_colored.update_layout(
        map_style="open-street-map",
        map_zoom=9,
        map_center={"lat": 37.7749, "lon": -122.4194},
        margin=dict(l=0, r=0, t=0, b=0),
        legend=dict(
            **ROUTE_LEGEND,
            x=1.01,
            y=1,
            xanchor="left",
            yanchor="top",
        ),
    )
    return fig_colored


def _offset_colored_route_geometries(routes_gdf):
    """Return a copy of route geometries offset for colored-map readability."""
    visible_routes_gdf = routes_gdf.copy()
    source_crs = visible_routes_gdf.crs or "EPSG:4326"

    visible_routes_gdf["geometry"] = [
        _offset_route_geometry(row["geometry"], COLOR_ROUTE_OFFSETS_METERS.get(row["Color"], 0), source_crs)
        for _, row in visible_routes_gdf.iterrows()
    ]
    return visible_routes_gdf


def _offset_route_geometry(geometry, offset_meters, source_crs):
    """Offset route geometry in meters, then convert back to map coordinates."""
    if offset_meters == 0:
        return geometry

    projected_geometry = gpd.GeoSeries([geometry], crs=source_crs).to_crs(epsg=3857).iloc[0]
    offset_lines = []

    for line in _iter_plot_lines(projected_geometry):
        offset_line = line.offset_curve(offset_meters, join_style=2)
        offset_lines.extend(_iter_plot_lines(offset_line))

    offset_geometry = MultiLineString(offset_lines)
    return gpd.GeoSeries([offset_geometry], crs="EPSG:3857").to_crs(source_crs).iloc[0]


def _iter_plot_lines(geometry):
    """Yield LineString parts from plotted line geometry."""
    if isinstance(geometry, LineString):
        yield geometry
    elif isinstance(geometry, MultiLineString):
        yield from geometry.geoms
    else:
        raise TypeError(f"Unsupported plot geometry type: {type(geometry)!r}")


def _build_ridership_route_map(routes_gdf, selected_route, filtered_stations):
    """Create the top map with black routes and ridership-scaled station markers.

    Args:
        routes_gdf: Route GeoDataFrame for the selected service year.
        selected_route: Route name to display, or `"all"` for every route.
        filtered_stations: Station GeoDataFrame after station dropdown filtering.
    """
    fig_black = go.Figure()
    filtered_routes = (
        routes_gdf
        if selected_route == "all"
        else routes_gdf[routes_gdf["route"] == selected_route]
    )
    _add_route_line_traces(fig_black, filtered_routes, line_width=2, line_color="black")

    fig_black.add_trace(
        go.Scattermap(
            lat=filtered_stations["lat"],
            lon=filtered_stations["lon"],
            mode="markers+text",
            marker=dict(
                size=filtered_stations["Ridership"] / 1000 + 3,
                color=filtered_stations["Ridership"],
                colorscale="YlOrRd",
                showscale=True,
                colorbar=dict(
                    title="Ridership",
                    x=0,
                    xanchor="right",
                    thickness=20,
                ),
            ),
            text=filtered_stations["Name"],
            textposition="top right",
            hoverinfo="text",
            name="BART Stations",
            showlegend=False,
        )
    )

    fig_black.update_layout(
        map_style="open-street-map",
        map_zoom=9,
        map_center={"lat": 37.7749, "lon": -122.4194},
        margin=dict(l=0, r=300, t=0, b=0),
        legend=dict(
            **ROUTE_LEGEND,
            x=1.01,
            y=.96,
            xanchor="left",
            yanchor="top",
        ),
    )
    return fig_black


def _build_ridership_bar_chart(filtered_ridership):
    """Create the station ridership summary bar chart.

    Args:
        filtered_ridership: Aggregated ridership rows for either all stations or
            the two selected stations.
    """
    if filtered_ridership.empty:
        bar_fig = go.Figure()
        bar_fig.update_layout(
            title="No Ridership Data Available",
            xaxis_title="Ridership",
            yaxis_title="Station Name",
        )
        return bar_fig

    chart_data = filtered_ridership.copy()
    chart_data["Ridership"] = pd.to_numeric(chart_data["Ridership"], errors="coerce").fillna(0)
    chart_data["Station"] = chart_data["Full Station Name"].fillna(chart_data["Entry Station"])

    if len(chart_data) > 2:
        chart_data = chart_data.nlargest(10, "Ridership").sort_values("Ridership")
        title = "Top 10 Stations by January 2018 Ridership"
    else:
        chart_data = chart_data.sort_values("Ridership")
        title = "Selected Station Ridership (January 2018)"

    return px.bar(
        chart_data,
        x="Ridership",
        y="Station",
        title=title,
        labels={"Ridership": "Total Ridership", "Station": "Station Name"},
        color="Ridership",
        color_continuous_scale="blues",
        template="plotly_white",
        orientation="h",
    )
# endregion


# region Entrypoint
if __name__ == "__main__":
    app.run(debug=True)
# endregion
