import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import dash
from dash import Dash, html, dash_table, dcc
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc
import geopandas as gpd
from shapely.geometry import LineString, MultiLineString, Point
from shapely.ops import nearest_points
from math import isclose
from pathlib import Path

DEBUG = False
ROOT = Path(__file__).parent.resolve()

DATA = ROOT / "data"
RAW = DATA / "raw"
PROCESSED = DATA / "processed"
PROCESSED.mkdir(parents=True, exist_ok=True)




bart = pd.read_csv(RAW / 'date-hour-soo-dest-2018.csv')
bart.head()


bart.loc[-1] = bart.columns
bart = bart.sort_index().reset_index(drop=True)
new_cols = ['Date', 'Hour', 'ORIGIN', 'DESTINATION', 'Number of Exits']
bart.columns = new_cols
bart.head()


bart['DESTINATION'].value_counts()

stations_gdf = gpd.read_file(RAW / 'BART_Stations_2025' / 'doc.kml', driver='KML')

stations_gdf.to_file(PROCESSED / 'stations.geojson', driver='GeoJSON')
stations_gdf.to_csv(PROCESSED / 'stations.csv', index=False)
stations_gdf.head()


stations_gdf = stations_gdf.drop(columns=['description', 'timestamp', 'begin', 'end', 'altitudeMode', 'tessellate', 'extrude', 'visibility', 'drawOrder',
                           'icon', 'snippet'], errors='ignore')
stations_gdf.head()


gdb_path = RAW / 'BART_Routes' / 'p20' / 'shortexercise1.gdb'
route_to_gdf = gpd.read_file(gdb_path, layer='BARTLine')

route_to_gdf.head()
route_to_gdf.to_file(PROCESSED / 'bart_routes.geojson', driver='GeoJSON')

routes_gdf = gpd.read_file(RAW / 'BART_Stations_2025' / 'bart_geojson.json')
routes_gdf.to_file(PROCESSED / 'bart_geojson.geojson', driver='GeoJSON')
routes_gdf.sort_values(by=['shape_leng'])

first_geom_0 = routes_gdf['geometry'].iloc[0]

new_subroute = first_geom_0.geoms[0]


first_geom_1 = route_to_gdf.iloc[0]['geometry']
last_two = [first_geom_1.geoms[-1].coords[-2], first_geom_1.geoms[-1].coords[-1]]


def find_nearest_coord(target_coord, line):
    """Finds the nearest point in a LineString to a given station coordinate."""
    point = Point(target_coord)  # Convert station coord to a point
    nearest_point = nearest_points(line, point)[1]  # Find closest point on line
    return nearest_point.coords[0]  # Return the coordinate

def drop_z(coord):
    return (coord[0], coord[1])  

def create_subroute(start_station, end_station):

    start_coords = drop_z(stations_gdf.loc[stations_gdf['Name'] == start_station, 'geometry'].values[0].coords[0])
    end_coords = drop_z(stations_gdf.loc[stations_gdf['Name'] == end_station, 'geometry'].values[0].coords[0])

    if DEBUG: print(f'Start Station: {start_station} -> {start_coords}')
    if DEBUG: print(f'End Station: {end_station} -> {end_coords}')

    matching_segments = []

    for _, row in route_to_gdf.iterrows():
        for line in row['geometry'].geoms:
            closest_start = find_nearest_coord(start_coords, line)
            closest_end = find_nearest_coord(end_coords, line)
            coords = list(line.coords)
            if closest_start in coords and closest_end in coords:
                if DEBUG: print(f"Found matching segment in route: {row['route']}")
                matching_segments.append(line)
                
    if not matching_segments:
        if DEBUG: print('No matching routes found')
        return None

    subroute_geom = MultiLineString(matching_segments)

    new_route_df = pd.DataFrame({
        'route': ['New Subroute'],
        'Color': ['Purple'],
        'geometry': [subroute_geom]
    })

    new_route_gdf = gpd.GeoDataFrame(new_route_df, geometry="geometry", crs=route_to_gdf.crs)

    return new_route_gdf


new_subroute = create_subroute("Oakland International Airport", "Fruitvale")  # Replace with desired stations

# Append to route dataset if successful
if new_subroute is not None:
    route_to_gdf = pd.concat([route_to_gdf, new_subroute], ignore_index=True)
    if DEBUG: print("Subroute successfully added!")

route_to_gdf = route_to_gdf.sort_values(by=['Shape_Length'])


color = ['Orange', 'Blue', 'Green', 'Red', 'Yellow']
route = ['Berryessa/North San Jose to Richmond', 'Dublin/Pleasanton to Daly City', 'Berryessa/North San Jose to Daly City', 'Richmond to Millbrae', 'Antioch to SFO International Airport']
route_to_gdf['Color'] = color
route_to_gdf['route'] = route


new_line = LineString(last_two)
new_multi = MultiLineString([new_line])

new_route_df = pd.DataFrame({
    "route": ["Coliseum to Oakland Int'l Airport"],  
    "Color": ["Grey"],  # Assign a unique color (e.g., pink)
    "geometry": [new_multi]
})

new_route_gdf = gpd.GeoDataFrame(new_route_df, geometry="geometry", crs=route_to_gdf.crs)

# Append the new route to route_to_gdf
route_to_gdf = pd.concat([route_to_gdf, new_route_gdf], ignore_index=True)


def coord_close(coord1, coord2, tol=1e-6):
    return abs(coord1[0] - coord2[0]) < tol and abs(coord1[1] - coord2[1]) < tol

station_coords = {row["Name"]: (row.geometry.x, row.geometry.y) for _, row in stations_gdf.iterrows()}

stations_to_add = {
    "Dublin/Pleasanton to Daly City": ["Lake Merritt", "12th St/Oakland City Center"]
}

updated_routes = []

for _, row in route_to_gdf.iterrows():
    route_name = row["route"].strip()
    geometry = row["geometry"]

    if route_name in stations_to_add:
        if DEBUG: print(f"Modifying route: {route_name}")
        
        line_strings = list(geometry.geoms)
        last_line = line_strings[-1]
        updated_coords = list(last_line.coords)
        
        if DEBUG: print(f"Before inserting stations: {updated_coords}")

        # Choose an insertion index
        insert_index = 0  # insert at beginning
        for i, station in enumerate(stations_to_add[route_name]):  
            if station in station_coords:
                coord = station_coords[station]
                if not any(coord_close(coord, existing) for existing in updated_coords):
                    updated_coords.insert(insert_index, coord)
                    if DEBUG: print(f"Inserted station {station} at index {insert_index}: {coord}")
                else:
                    if DEBUG: print(f"Skipped duplicate station {station}")

        new_line = LineString(updated_coords)
        updated_route = MultiLineString(line_strings[:-1] + [new_line])
        updated_routes.append(updated_route)

        if DEBUG: print(f"Updated geometry: {updated_route}")
    else:
        updated_routes.append(geometry)

assert len(updated_routes) == len(route_to_gdf), "Geometry mismatch after update!"
route_to_gdf["geometry"] = updated_routes

if DEBUG: print("Updated routes with inserted station coordinates successfully!")


stations_remove = {
    'Richmond to Millbrae': ['Fruitvale', 'Lake Merritt', 'Oakland International Airport']
}

routes_remove = []

for _, row in route_to_gdf.iterrows():
    route_name = row['route'].strip()
    geometry = row['geometry']

    if route_name in stations_remove:
        if DEBUG: print(f"Modifying route: {route_name}")
        line_strings = list(geometry.geoms)
        update_segs = []

        for line in line_strings:
            coords = list(line.coords)
            if DEBUG: print(f"Original segment coords: {coords[:5]} ...")
            station_coords_2D = {st: (coord[0], coord[1]) for st, coord in station_coords.items()}

            # Filter out station coordinates that need to be removed
            remove_coords = [station_coords_2D[st] for st in stations_remove[route_name] if st in station_coords_2D]
            if DEBUG: print(f"Coordinates to remove: {remove_coords}")  # Debugging output

            # Filter out station coordinates that need to be removed
            filtered_coords = [coord for coord in coords if not any(coord_close(coord, rem) for rem in remove_coords)]

            if DEBUG: print(f"Filtered segment coords: {filtered_coords[:5]} ...")

            if len(filtered_coords) > 1:
                update_segs.append(LineString(filtered_coords))
            else:
                if DEBUG: print(f"Skipping segment due to insufficient points: {filtered_coords}")
        route_remove = MultiLineString(update_segs) if update_segs else None
        routes_remove.append(route_remove)
        if DEBUG: print(f"Updated geometry: {route_remove}")
    else:
        routes_remove.append(geometry)

route_to_gdf['geometry'] = routes_remove


file_path = RAW / 'ridership_2018' / 'Ridership_201801.xlsx'
# Load dataset WITHOUT setting a header (treat all as raw data)
df_ridership = pd.read_excel(file_path, sheet_name=0, header=None)


entry_stations = df_ridership.iloc[1, 1:].tolist()  # Second row (entry stations), excluding first column
exit_stations = df_ridership.iloc[2:, 0].tolist()  # First column (exit stations), excluding first two rows

# Convert numeric station names (e.g., 16 → "16th")
entry_stations = [
    str(int(name)) + "th" if isinstance(name, (int, float)) and not pd.isna(name) else str(name)
    for name in entry_stations
]
exit_stations = [
    str(int(name)) + "th" if isinstance(name, (int, float)) and not pd.isna(name) else str(name)
    for name in exit_stations
]

# **Ensure station names are unique**
seen = set()
for i, name in enumerate(entry_stations):
    while name in seen:
        name += "_dup"  # Append '_dup' to handle duplicates
    seen.add(name)
    entry_stations[i] = name

# **Step 2: Remove header rows & columns, keep only data**
df_ridership_clean = df_ridership.iloc[2:, 1:].copy()  # Remove first two rows (labels) and first column (labels)
df_ridership_clean.columns = entry_stations  # Assign correct entry station names
df_ridership_clean.insert(0, "Exit Station", exit_stations)  # Insert exit station column

# Convert ridership data to numeric
df_ridership_clean.iloc[:, 1:] = df_ridership_clean.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")

if DEBUG: print("Cleaned Data (First 5 Rows):\n", df_ridership_clean.head())

# **Step 3: Convert to long format for visualization**
df_ridership_long = df_ridership_clean.melt(
    id_vars=["Exit Station"],
    var_name="Entry Station",
    value_name="Ridership"
)

# Convert to categorical for proper ordering
df_ridership_long["Entry Station"] = pd.Categorical(df_ridership_long["Entry Station"], categories=entry_stations, ordered=True)
df_ridership_long["Exit Station"] = pd.Categorical(df_ridership_long["Exit Station"], categories=exit_stations, ordered=True)

# **Final Long Format Data**
if DEBUG: print("Final Long Format Data (First 5 Rows):\n", df_ridership_long.head())


station_mapping = {
    "RM": "Richmond",
    "EN": "El Cerrito del Norte",
    "EP": "El Cerrito Plaza",
    "NB": "North Berkeley",
    "BK": "Berkeley",
    "AS": "Ashby",
    "MA": "MacArthur",
    "19th": "19th St/Oakland",
    "12th": "12th St/Oakland City Center",
    "LM": "Lake Merritt",
    "FV": "Fruitvale",
    "CL": "Coliseum",
    "SL": "San Leandro",
    "BF": "Bay Fair",
    "HY": "Hayward",
    "SH": "South Hayward",
    "UC": "Union City",
    "FM": "Fremont",
    "CN": "Concord",
    "PH": "Pleasant Hill/Contra Costa Centre",
    "WC": "Walnut Creek",
    "LF": "Lafayette",
    "OR": "Orinda",
    "RR": "Rockridge",
    "OW": "West Oakland",
    "EM": "Embarcadero",
    "MT": "Montgomery St",
    "PL": "Powell St",
    "CC": "Civic Center/UN Plaza",
    "16th": "16th St/Mission",
    "24th": "24th St/Mission",
    "GP": "Glen Park",
    "BP": "Balboa Park",
    "DC": "Daly City",
    "CM": "Colma",
    "CV": "Castro Valley",
    "ED": "Dublin/Pleasanton",
    "NC": "North Concord/Martinez",
    "WP": "West Dublin/Pleasanton",
    "SS": "South San Francisco",
    "SB": "San Bruno",
    "SO": "San Francisco International Airport",
    "MB": "Millbrae",
    "WD": "Warm Springs/South Fremont",
    "OA": "Oakland International Airport",
    "WS": "Coliseum/Airport Connector",
    #"Exits": "Exits"  # Keeping "Exits" as-is
}



full_station_names = stations_gdf["Name"].unique().tolist()
station_ridership = df_ridership_long.groupby("Entry Station")["Ridership"].sum().reset_index()
station_ridership = station_ridership[station_ridership["Entry Station"] != "Exits"]
station_ridership['Full Station Name'] = station_ridership['Entry Station'].map(station_mapping)
# Get station abbreviations from station_ridership
station_abbreviations = station_ridership["Entry Station"].unique().tolist()

stations_gdf["lat"] = stations_gdf["geometry"].y  # Extract latitude
stations_gdf["lon"] = stations_gdf["geometry"].x  # Extract longitude

stations_gdf = stations_gdf.merge(station_ridership, left_on="Name", right_on="Full Station Name", how="left")
stations_gdf["Ridership"] = stations_gdf["Ridership"].fillna(0) 
# Route selection dropdown options
route_options = [{"label": route, "value": route} for route in route_to_gdf["route"].unique()]
route_options.insert(0, {"label": "All Routes", "value": "all"})  # Add "All Routes" option


# Initialize Dash app
app = dash.Dash(__name__)

app.layout = html.Div([
    html.H1("BART System Map & Ridership Data", style={"textAlign": "center", "color": "black"}),

    # TOP: Black route map + dropdowns
    html.Div([
        html.Div([
            dcc.Graph(id="bart-map-black", style={"height": "500px", "border": "1px solid lightgrey", "position": "relative"}),
        ], style={"width": "75%", "display": "inline-block"}),#, "padding": "5px"}),

            html.Div([
                html.H4("Dash App Description"),
                html.P("This Dash App displays BART Routes with interactive elements."),
                html.P("The dropdown menus allows control for the routes displayed on the graph and a comparison between aggregate ridership over the 2018 calendar year."),
                html.P("The ridership intensity is mapped according to circle size and color, being commensurately adjusted as controlled by the selected stations."),
                html.P("Ridership is clearly most intense at the first few stops in the city displaying the overall commuting trends for BART ridership."),
            ], style={
                "width": "400px", 
                "padding": "20px", 
                "backgroundColor": "lightgrey", 
                "border": "1px solid grey",
                "borderRadius": "10px",
                "marginLeft": "15px",
            }),

        
        html.Div([
            html.Div(style={'height': '200px'}),  # Empty space for legend alignment

            html.Label("Select Route:"),
            dcc.Dropdown(
                id="route-dropdown",
                options=route_options,
                value="all",
                clearable=False
            ),
            html.Br(),

            html.Label("Select Two Stations:"),
            dcc.Dropdown(
                id="station-1-dropdown",
                options=[{"label": name, "value": name} for name in stations_gdf["Name"]],
                clearable=True
            ),
            dcc.Dropdown(
                id="station-2-dropdown",
                options=[{"label": name, "value": name} for name in stations_gdf["Name"]],
                clearable=True
            ),
        ], style={
            "width": "250px",
            "display": "inline-block",
            "verticalAlign": "top",
            "padding": "10px",
            "position": "relative",
            "left": "-650px"
        }),
    ], style={"position": "relative", "height": "500px", "display": "flex", "justifyContent": "space-between", "alignItems": "flex-start"}),

    html.Br(),
    html.H3("Route Map (Colored) and Ridership Comparison"),

    # BOTTOM: Colored route map + bar chart
    html.Div([
        dcc.Graph(id="bart-map-colored", style={"width": "50%", "display": "inline-block"}),
        dcc.Graph(id="ridership-bar-chart", style={"width": "50%", "display": "inline-block"}),
    ], style={"display": "flex", "justifyContent": "center"}),
])



@app.callback(
    [
        Output("bart-map-colored", "figure"),
        Output("bart-map-black", "figure"),
        Output("ridership-bar-chart", "figure"),
    ],
    [
        Input("route-dropdown", "value"),
        Input("station-1-dropdown", "value"),
        Input("station-2-dropdown", "value")
    ]
)
def update_maps(selected_route, station1, station2):
    # **Initialize Figures**
    fig_colored = go.Figure()

    # **Apply Station Filtering (Only show rides between selected stations)**
    # Reverse the station_mapping dictionary to get full name -> abbreviation
    reverse_mapping = {v: k for k, v in station_mapping.items()}

    # In your callback, right before filtering station_ridership:
    if station1 and station2:
        station1_abbr = reverse_mapping.get(station1)
        station2_abbr = reverse_mapping.get(station2)

        filtered_stations = stations_gdf[stations_gdf["Name"].isin([station1, station2])]
        filtered_ridership = station_ridership[station_ridership["Entry Station"].isin([station1_abbr, station2_abbr])]
    else:
        filtered_stations = stations_gdf
        filtered_ridership = station_ridership


    # filtered_routes = route_to_gdf if selected_route == "all" else route_to_gdf[route_to_gdf["route"] == selected_route]

    # **Add Route Lines (Colored & Black)**
    for _, row in route_to_gdf.iterrows():
        for line in row["geometry"].geoms:
            lat_values = [p[1] for p in line.coords]
            lon_values = [p[0] for p in line.coords]

            # Colored Routes
            fig_colored.add_trace(go.Scattermapbox(
                lat=lat_values,
                lon=lon_values,
                mode="lines",
                line=dict(width=3, color=row["Color"]),
                text=row["route"],
                hoverinfo="text",
                name=row["route"]
            ))
    
    fig_colored.add_trace(go.Scattermapbox(
        lat=stations_gdf['lat'],
        lon=stations_gdf['lon'],
        mode='markers+text',
        marker=dict(
            size=6,
            color='black'
        ),
        text=stations_gdf['Name'],
        textposition='top right',
        hoverinfo='text',
        name='BART Stations'
    ))

    fig_colored.update_layout(
        mapbox_style="open-street-map",
        mapbox_zoom=9,
        mapbox_center={"lat": 37.7749, "lon": -122.4194},
        margin=dict(l=0, r=0, t=0, b=0)
    )
    
    fig_black = go.Figure()
    
    filtered_routes = route_to_gdf if selected_route == "all" else route_to_gdf[route_to_gdf["route"] == selected_route]
    
    for _, row in filtered_routes.iterrows():
        for line in row["geometry"].geoms:
            lat_values = [p[1] for p in line.coords]
            lon_values = [p[0] for p in line.coords]

            fig_black.add_trace(go.Scattermapbox(
                lat=lat_values,
                lon=lon_values,
                mode="lines",
                line=dict(width=2, color="black"),
                text=row["route"],
                hoverinfo="text",
                name=row["route"]
            ))
    
    fig_black.add_trace(go.Scattermapbox(
        lat=filtered_stations['lat'],
        lon=filtered_stations['lon'],
        mode='markers+text',
        marker=dict(
            size=filtered_stations['Ridership'] / 1000 + 3,
            color=filtered_stations['Ridership'],
            colorscale='YlOrRd',
            showscale=True,
            colorbar=dict(
                title='Ridership',
                x=1.01,
                xanchor='left',
                thickness=20
            ),
        ),
        text=filtered_stations['Name'],
        textposition='top right',
        hoverinfo='text',
        name='BART Stations'
    ))

    fig_black.update_layout(
        mapbox_style="open-street-map",
        mapbox_zoom=9,
        mapbox_center={"lat": 37.7749, "lon": -122.4194},
        margin=dict(l=0, r=0, t=0, b=0),
        legend=dict(
            x=1.25,
            y=1,
            xanchor='left',
            yanchor='top',
            bgcolor='rgba(255,255,255,0.9)',
            bordercolor='lightgrey',
            borderwidth=1,
            font=dict(size=12)
        )
    )
    # **Update Bar Chart (Filtered for Two Stations)**
    if filtered_ridership.empty:
        bar_fig = go.Figure()
        bar_fig.update_layout(title="No Ridership Data Available", xaxis_title="Station", yaxis_title="Ridership")
    else:
        bar_fig = px.bar(
            filtered_ridership,
            x="Entry Station",
            y="Ridership",
            title="Ridership Comparison by Station",
            labels={"Ridership": "Total Ridership", "Entry Station": "Station"},
            color="Ridership",
            color_continuous_scale="blues",
            template='plotly_white'
        )

    return fig_colored, fig_black, bar_fig

if __name__ == '__main__':
    app.run(debug=True)
