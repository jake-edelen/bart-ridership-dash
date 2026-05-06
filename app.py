"""Dash entrypoint for the BART ridership dashboard.

This module should stay focused on presentation: layout, callback wiring, and
Plotly figure construction. Data loading, station mapping, and route geometry
corrections live under `src/` so future data/modeling work can scale without
turning the Dash file back into a preprocessing script.
"""

# region Imports
from calendar import month_name
from functools import lru_cache

import dash
import geopandas as gpd
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash import dcc, html
from dash.dependencies import Input, Output, State
from shapely.geometry import LineString, MultiLineString

from src.data_loader import (
    attach_station_ridership,
    get_available_ridership_periods,
    load_station_ridership_for_period,
    load_app_data,
)
from src.route_builder import build_current_routes
# endregion


# region Runtime Data
DEFAULT_RIDERSHIP_YEAR = 2018
DEFAULT_RIDERSHIP_MONTH = 1

app_data = load_app_data()
stations_gdf = app_data.stations_gdf
raw_routes_gdf = app_data.raw_routes_gdf
station_ridership = app_data.station_ridership
station_mapping = app_data.station_mapping
available_ridership_periods = get_available_ridership_periods()

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


def _ridership_month_options_for_year(year):
    """Build month dropdown options for a selected ridership year."""
    months = [
        month
        for period_year, month in available_ridership_periods
        if period_year == int(year)
    ]
    return [
        {"label": month_name[month], "value": month}
        for month in sorted(months)
    ]


@lru_cache(maxsize=None)
def _station_ridership_for_period(year, month):
    """Return station-level ridership for the selected monthly period."""
    return load_station_ridership_for_period(int(year), int(month))


def _stations_for_ridership_period(year, month):
    """Return station geometry with ridership values for one monthly period."""
    base_stations = stations_gdf.drop(
        columns=["Entry Station", "Ridership", "Full Station Name"],
        errors="ignore",
    )
    return attach_station_ridership(
        base_stations,
        _station_ridership_for_period(year, month),
    )


def _ridership_period_label(year, month):
    """Return a readable month/year label for figure titles and hovers."""
    return f"{month_name[int(month)]} {int(year)}"


def _station_options_for_ridership_period(year, month):
    """Build station dropdown options from stations present in a ridership period."""
    station_names = (
        _station_ridership_for_period(year, month)["Full Station Name"]
        .dropna()
        .unique()
    )
    available_station_names = set(stations_gdf["Name"])
    visible_station_names = sorted(
        name
        for name in station_names
        if name in available_station_names
    )
    return [
        {"label": name, "value": name}
        for name in visible_station_names
    ]


route_to_gdf = _routes_for_year(DEFAULT_RIDERSHIP_YEAR)
route_options = _route_options_for_routes(route_to_gdf)
station_options = _station_options_for_ridership_period(
    DEFAULT_RIDERSHIP_YEAR,
    DEFAULT_RIDERSHIP_MONTH,
)

ridership_year_options = [
    {"label": str(year), "value": year}
    for year in sorted({year for year, _ in available_ridership_periods})
]
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

RIDERSHIP_COLOR_SCALE = "YlOrRd"

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
                                "monthly station ridership."
                            ),
                            html.P(
                                "Ridership year also controls the route service year shown on the maps."
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
                            html.Label("Select Ridership Year:"),
                            dcc.Dropdown(
                                id="ridership-year-dropdown",
                                options=ridership_year_options,
                                value=DEFAULT_RIDERSHIP_YEAR,
                                clearable=False,
                            ),
                            html.Br(),
                            html.Label("Select Ridership Month:"),
                            dcc.Dropdown(
                                id="ridership-month-dropdown",
                                options=_ridership_month_options_for_year(DEFAULT_RIDERSHIP_YEAR),
                                value=DEFAULT_RIDERSHIP_MONTH,
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
                                options=station_options,
                                clearable=True,
                            ),
                            html.Br(),
                            dcc.Dropdown(
                                id="station-2-dropdown",
                                options=station_options,
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
app.validation_layout = create_layout()
# endregion


# region Callbacks
@app.callback(
    [
        Output("route-dropdown", "options"),
        Output("route-dropdown", "value"),
    ],
    [Input("ridership-year-dropdown", "value")],
    [State("route-dropdown", "value")],
)
def update_route_dropdown(ridership_year, selected_route):
    """Refresh route choices so they match the selected ridership year."""
    routes_gdf = _routes_for_year(ridership_year)
    options = _route_options_for_routes(routes_gdf)
    valid_routes = {option["value"] for option in options}
    value = selected_route if selected_route in valid_routes else "all"
    return options, value


@app.callback(
    [
        Output("ridership-month-dropdown", "options"),
        Output("ridership-month-dropdown", "value"),
    ],
    [Input("ridership-year-dropdown", "value")],
    [State("ridership-month-dropdown", "value")],
)
def update_ridership_month_dropdown(selected_year, selected_month):
    """Refresh month choices so they match the selected ridership year."""
    options = _ridership_month_options_for_year(selected_year)
    valid_months = {option["value"] for option in options}
    value = selected_month if selected_month in valid_months else options[0]["value"]
    return options, value


@app.callback(
    [
        Output("station-1-dropdown", "options"),
        Output("station-1-dropdown", "value"),
        Output("station-2-dropdown", "options"),
        Output("station-2-dropdown", "value"),
    ],
    [
        Input("ridership-year-dropdown", "value"),
        Input("ridership-month-dropdown", "value"),
    ],
    [
        State("station-1-dropdown", "value"),
        State("station-2-dropdown", "value"),
    ],
)
def update_station_dropdowns(ridership_year, ridership_month, station1, station2):
    """Refresh station choices so they match the selected ridership period."""
    ridership_year = ridership_year or DEFAULT_RIDERSHIP_YEAR
    ridership_month = ridership_month or DEFAULT_RIDERSHIP_MONTH
    options = _station_options_for_ridership_period(ridership_year, ridership_month)
    valid_stations = {option["value"] for option in options}
    station1 = station1 if station1 in valid_stations else None
    station2 = station2 if station2 in valid_stations else None
    return options, station1, options, station2


@app.callback(
    [
        Output("bart-map-colored", "figure"),
        Output("bart-map-black", "figure"),
        Output("ridership-bar-chart", "figure"),
    ],
    [
        Input("route-dropdown", "value"),
        Input("ridership-year-dropdown", "value"),
        Input("ridership-month-dropdown", "value"),
        Input("station-1-dropdown", "value"),
        Input("station-2-dropdown", "value"),
    ],
)
def update_maps(
    selected_route,
    ridership_year=None,
    ridership_month=None,
    station1=None,
    station2=None,
):
    """Update route maps and ridership chart from selected UI filters.

    Args:
        selected_route: Route name from the route dropdown, or `"all"`.
        ridership_year: Year used for station ridership and route geometry.
        ridership_month: Month used for station ridership display.
        station1: Optional full station name from the first station dropdown.
        station2: Optional full station name from the second station dropdown.

    Returns:
        A tuple of Plotly figures for the colored route map, ridership map, and
        ridership bar chart.
    """
    ridership_year = ridership_year or DEFAULT_RIDERSHIP_YEAR
    ridership_month = ridership_month or DEFAULT_RIDERSHIP_MONTH
    period_label = _ridership_period_label(ridership_year, ridership_month)
    selected_station_ridership = _station_ridership_for_period(ridership_year, ridership_month)
    stations_for_period = _stations_for_ridership_period(ridership_year, ridership_month)
    reverse_mapping = {v: k for k, v in station_mapping.items()}

    if station1 and station2:
        station1_abbr = reverse_mapping.get(station1)
        station2_abbr = reverse_mapping.get(station2)
        filtered_stations = stations_for_period[stations_for_period["Name"].isin([station1, station2])]
        filtered_ridership = selected_station_ridership[
            selected_station_ridership["Entry Station"].isin([station1_abbr, station2_abbr])
        ]
    else:
        filtered_stations = stations_for_period
        filtered_ridership = selected_station_ridership

    selected_routes_gdf = _routes_for_year(ridership_year)
    valid_routes = set(selected_routes_gdf["route"])
    if selected_route != "all" and selected_route not in valid_routes:
        selected_route = "all"

    fig_colored = _build_colored_route_map(selected_routes_gdf, ridership_year)
    fig_black = _build_ridership_route_map(
        selected_routes_gdf,
        selected_route,
        filtered_stations,
        period_label,
        _ridership_color_range(stations_for_period["Ridership"]),
    )
    bar_fig = _build_ridership_bar_chart(
        filtered_ridership,
        period_label,
        _ridership_color_range(selected_station_ridership["Ridership"]),
    )

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

def _build_colored_route_map(routes_gdf, selected_year=None):
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
            hovertemplate="%{text}<extra></extra>",
            name="BART Stations",
            showlegend=True, 
        )
    )

    title_year = f"{selected_year} " if selected_year else ""
    fig_colored.update_layout(
        map_style="open-street-map",
        map_zoom=9,
        map_center={"lat": 37.7749, "lon": -122.4194},
        title=dict(
            text=f"Route Reference Map ({title_year}service)",
            x=0.5,
            xanchor="center",
            font=dict(size=16),
        ),
        margin=dict(l=0, r=0, t=58, b=0),
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


def _scale_ridership_marker_sizes(ridership_values, min_size=5, max_size=34):
    """Scale station marker sizes by relative ridership without runaway circles."""
    ridership = pd.to_numeric(ridership_values, errors="coerce").fillna(0)
    max_ridership = ridership.max()

    if max_ridership <= 0:
        return pd.Series(min_size, index=ridership.index)

    return min_size + (ridership / max_ridership).pow(0.5) * (max_size - min_size)


def _ridership_color_range(station_ridership_values):
    """Return a stable color range for one ridership period."""
    ridership = pd.to_numeric(station_ridership_values, errors="coerce").fillna(0)
    max_ridership = ridership.max()

    if max_ridership <= 0:
        return 0, 1

    return 0, max_ridership


def _build_ridership_route_map(
    routes_gdf,
    selected_route,
    filtered_stations,
    ridership_period_label="January 2018",
    color_range=None,
):
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
    color_min, color_max = color_range or _ridership_color_range(filtered_stations["Ridership"])
    _add_route_line_traces(
        fig_black,
        filtered_routes,
        line_width=2,
        line_color="rgba(35, 35, 35, 0.45)",
    )

    fig_black.add_trace(
        go.Scattermap(
            lat=filtered_stations["lat"],
            lon=filtered_stations["lon"],
            mode="markers+text",
            marker=dict(
                size=_scale_ridership_marker_sizes(filtered_stations["Ridership"]),
                color=filtered_stations["Ridership"],
                colorscale=RIDERSHIP_COLOR_SCALE,
                cmin=color_min,
                cmax=color_max,
                showscale=True,
                colorbar=dict(
                    title=f"{ridership_period_label}<br>Ridership",
                    x=0,
                    xanchor="right",
                    thickness=20,
                ),
            ),
            text=filtered_stations["Name"],
            customdata=filtered_stations["Ridership"],
            textposition="top right",
            hovertemplate=f"%{{text}}<br>{ridership_period_label} ridership: %{{customdata:,.0f}}<extra></extra>",
            name="BART Stations",
            showlegend=False,
        )
    )

    fig_black.update_layout(
        map_style="open-street-map",
        map_zoom=9,
        map_center={"lat": 37.7749, "lon": -122.4194},
        title=dict(
            text=f"Station Ridership Map ({ridership_period_label})",
            x=0.5,
            xanchor="center",
            font=dict(size=16),
        ),
        margin=dict(l=0, r=300, t=42, b=0),
        legend=dict(
            **ROUTE_LEGEND,
            x=1.01,
            y=.96,
            xanchor="left",
            yanchor="top",
        ),
    )
    return fig_black


def _build_ridership_bar_chart(
    filtered_ridership,
    ridership_period_label="January 2018",
    color_range=None,
):
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
    color_min, color_max = color_range or _ridership_color_range(chart_data["Ridership"])

    if len(chart_data) > 2:
        chart_data = chart_data.nlargest(10, "Ridership").sort_values("Ridership")
        title = f"Top 10 Stations by {ridership_period_label} Ridership"
        comparison_note = None
    else:
        chart_data = chart_data.sort_values("Ridership")
        title = f"Selected Station Comparison ({ridership_period_label})"
        comparison_note = _selected_station_difference_note(chart_data)

    bar_fig = px.bar(
        chart_data,
        x="Ridership",
        y="Station",
        title=title,
        labels={"Ridership": "Total Ridership", "Station": "Station Name"},
        color="Ridership",
        color_continuous_scale=RIDERSHIP_COLOR_SCALE,
        range_color=(color_min, color_max),
        template="plotly_white",
        orientation="h",
    )
    bar_fig.update_layout(
        xaxis_title=f"{ridership_period_label} Ridership",
        yaxis_title="Station",
        coloraxis_colorbar=dict(title="Ridership"),
        margin=dict(l=120, r=20, t=58, b=45),
    )
    if comparison_note:
        bar_fig.add_annotation(
            text=comparison_note,
            xref="paper",
            yref="paper",
            x=1,
            y=1.12,
            xanchor="right",
            yanchor="bottom",
            showarrow=False,
            font=dict(size=12),
        )
    return bar_fig


def _selected_station_difference_note(chart_data):
    """Return a concise ridership difference note for selected station comparisons."""
    if len(chart_data) != 2:
        return None

    sorted_data = chart_data.sort_values("Ridership", ascending=False)
    high_station, low_station = sorted_data["Station"].tolist()
    high_ridership, low_ridership = sorted_data["Ridership"].tolist()
    difference = high_ridership - low_ridership

    return f"{high_station} is {difference:,.0f} riders higher than {low_station}"
# endregion


# region Entrypoint
if __name__ == "__main__":
    app.run(debug=True)
# endregion
