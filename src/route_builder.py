"""Curated route geometry builder for the BART dashboard.

The raw GDB route source has useful newer geometry, but no route names or
colors. This module treats those rows as infrastructure-like segments, labels
them by their station endpoints, and composes visible service rows from
date-aware service patterns.
"""

# region Imports
import geopandas as gpd
from shapely.geometry import LineString, MultiLineString, Point
from shapely.ops import substring

from .service_patterns import get_service_patterns_for_date
from .station_mapping import STATION_MAPPING, normalize_station_code
# endregion


# region Segment Definitions
SEGMENT_ENDPOINTS = {
    frozenset(("12th St/Oakland City Center", "Richmond")): "richmond_oakland",
    frozenset(("Fruitvale", "Dublin/Pleasanton")): "dublin_fruitvale",
    frozenset(("Bay Fair", "Berryessa/North San José")): "berryessa_bay_fair",
    frozenset(("Millbrae", "Oakland International Airport")): "millbrae_airport_oakland",
    frozenset(("MacArthur", "Antioch")): "antioch_macarthur",
}

# endregion


# region Public Builder
def build_current_routes(raw_routes_gdf, stations_gdf, service_date="2025-01-01"):
    """Build visible route rows for service active on `service_date`.

    Args:
        raw_routes_gdf: Raw GDB route geometries with `Shape_Length` and geometry.
        stations_gdf: Station point GeoDataFrame used to classify endpoints.
        service_date: Date used to select active service patterns.

    Returns:
        GeoDataFrame with one row per visible service pattern.
    """
    route_line_index = build_route_line_index(raw_routes_gdf, stations_gdf)
    service_patterns = get_service_patterns_for_date(service_date)
    _raise_for_missing_pattern_station_geometry(stations_gdf, service_patterns)
    route_rows = []

    for pattern in service_patterns:
        route_rows.append(
            {
                "service_id": pattern.service_id,
                "route": pattern.display_name,
                "Color": pattern.line_color,
                "valid_from": pattern.valid_from,
                "valid_to": pattern.valid_to,
                "station_codes": pattern.station_codes,
                "geometry": build_pattern_geometry(pattern, route_line_index),
            }
        )

    return gpd.GeoDataFrame(route_rows, geometry="geometry", crs=raw_routes_gdf.crs)


def build_pattern_geometry(pattern, route_line_index):
    """Build a route geometry by walking consecutive stations in a pattern."""
    station_names = [_station_name_for_code(code) for code in pattern.station_codes]
    lines = []

    for start_name, end_name in zip(station_names, station_names[1:]):
        lines.append(
            _extract_line_between_stations(
                start_name,
                end_name,
                route_line_index,
                pattern.service_id,
            )
        )

    return MultiLineString(lines)


def build_route_line_index(raw_routes_gdf, stations_gdf, tolerance=1e-3):
    """Index raw infrastructure lines by station projections along each line.

    The raw GDB rows are coarse line strings that span many stations. This
    index records where each station falls along each raw line so route patterns
    can be assembled by extracting the actual sub-line from the previous
    station to the next station.
    """
    station_lookup = _station_lookup(stations_gdf)
    line_index = []

    for _, row in raw_routes_gdf.iterrows():
        for line in _iter_lines(row["geometry"]):
            station_projections = {}

            for station_name, station_point in station_lookup.items():
                if line.distance(station_point) <= tolerance:
                    station_projections[station_name] = line.project(station_point)

            if len(station_projections) >= 2:
                line_index.append(
                    {
                        "line": line,
                        "station_projections": station_projections,
                    }
                )

    line_index.extend(_synthetic_connector_lines(station_lookup))
    return tuple(line_index)


def build_infrastructure_segments(raw_routes_gdf, stations_gdf):
    """Classify raw GDB geometries as reusable infrastructure segments.

    Segment identity comes from nearest station endpoints, not row order or
    `Shape_Length`.
    """
    station_lookup = _station_lookup(stations_gdf)
    segments = {}

    for _, row in raw_routes_gdf.iterrows():
        endpoint_names = _nearest_endpoint_station_names(row["geometry"], station_lookup)
        segment_id = SEGMENT_ENDPOINTS.get(frozenset(endpoint_names))

        if segment_id is None:
            raise ValueError(
                "Unrecognized route geometry endpoints: "
                f"{endpoint_names[0]} <-> {endpoint_names[1]}"
            )

        segments[segment_id] = row["geometry"]

    missing_segments = set(SEGMENT_ENDPOINTS.values()) - set(segments)
    if missing_segments:
        raise ValueError(f"Missing expected infrastructure segments: {sorted(missing_segments)}")

    return segments
# endregion


# region Validation
def validate_station_geometry_coverage(stations_gdf, station_mapping, station_codes):
    """Return station-code mapping/geometry coverage problems.

    Args:
        stations_gdf: Station GeoDataFrame from `doc.kml`.
        station_mapping: Canonical code-to-full-name mapping.
        station_codes: Codes that should be backed by station coordinates.

    Returns:
        Dictionary with `missing_mapping`, `missing_geometry`, and
        `null_geometry` lists.
    """
    station_names = set(stations_gdf["Name"].dropna())
    missing_mapping = []
    missing_geometry = []
    null_geometry = []

    for code in station_codes:
        station_name = station_mapping.get(code)
        if station_name is None:
            missing_mapping.append(code)
            continue

        matches = stations_gdf[stations_gdf["Name"] == station_name]
        if matches.empty:
            missing_geometry.append((code, station_name))
            continue

        if matches["geometry"].isna().any() or any(geom is None for geom in matches["geometry"]):
            null_geometry.append((code, station_name))

    return {
        "missing_mapping": tuple(sorted(missing_mapping)),
        "missing_geometry": tuple(sorted(missing_geometry)),
        "null_geometry": tuple(sorted(null_geometry)),
        "available_station_count": len(station_names),
    }


def validate_extension_station_geometry(raw_routes_gdf, stations_gdf, station_mapping, station_codes, tolerance=1e-8):
    """Return extension stations that are not on or near any raw route geometry."""
    problems = []

    for code in station_codes:
        station_name = station_mapping[code]
        station_rows = stations_gdf[stations_gdf["Name"] == station_name]
        if station_rows.empty:
            problems.append((code, station_name, None))
            continue

        station_point = station_rows.iloc[0]["geometry"]
        distance = min(float(geometry.distance(station_point)) for geometry in raw_routes_gdf.geometry)
        if distance > tolerance:
            problems.append((code, station_name, distance))

    return tuple(problems)
# endregion


# region Helpers
def _raise_for_missing_pattern_station_geometry(stations_gdf, service_patterns):
    """Fail route building if active pattern station coordinates are incomplete."""
    station_codes = sorted(
        {
            normalize_station_code(code)
            for pattern in service_patterns
            for code in pattern.station_codes
        }
    )
    coverage = validate_station_geometry_coverage(stations_gdf, STATION_MAPPING, station_codes)
    if coverage["missing_mapping"] or coverage["missing_geometry"] or coverage["null_geometry"]:
        raise ValueError(f"Missing service-pattern station geometry coverage: {coverage}")


def _station_lookup(stations_gdf):
    """Return station name to point geometry mapping."""
    return {
        row["Name"]: row.geometry
        for _, row in stations_gdf.iterrows()
        if row.geometry is not None
    }


def _station_name_for_code(code):
    """Return the full station name for a canonical service-pattern code."""
    normalized_code = normalize_station_code(code)
    station_name = STATION_MAPPING.get(normalized_code)
    if station_name is None:
        raise ValueError(f"No station mapping found for service-pattern code {code!r}.")
    return station_name


def _extract_line_between_stations(start_name, end_name, route_line_index, service_id):
    """Extract the infrastructure sub-line connecting two station names."""
    candidates = []

    for record in route_line_index:
        station_projections = record["station_projections"]
        if start_name in station_projections and end_name in station_projections:
            start_distance = station_projections[start_name]
            end_distance = station_projections[end_name]
            candidates.append(
                (
                    abs(end_distance - start_distance),
                    record["line"],
                    start_distance,
                    end_distance,
                )
            )

    if not candidates:
        raise ValueError(
            f"No infrastructure line found for {service_id}: "
            f"{start_name} -> {end_name}."
        )

    _, line, start_distance, end_distance = min(candidates, key=lambda candidate: candidate[0])
    segment = substring(line, start_distance, end_distance)

    if segment.is_empty or not isinstance(segment, LineString):
        raise ValueError(
            f"Could not extract route geometry for {service_id}: "
            f"{start_name} -> {end_name}."
        )

    return segment


def _iter_lines(geometry):
    """Yield LineString parts from LineString or MultiLineString geometry."""
    if isinstance(geometry, LineString):
        yield geometry
    elif isinstance(geometry, MultiLineString):
        yield from geometry.geoms
    else:
        raise TypeError(f"Unsupported route geometry type: {type(geometry)!r}")


def _synthetic_connector_lines(station_lookup):
    """Return station-correct connector lines missing from the raw GDB source."""
    connector_pairs = (
        ("Coliseum", "Oakland International Airport"),
    )
    lines = []

    for start_name, end_name in connector_pairs:
        start_point = station_lookup[start_name]
        end_point = station_lookup[end_name]
        line = LineString([start_point, end_point])
        lines.append(
            {
                "line": line,
                "station_projections": {
                    start_name: 0.0,
                    end_name: line.length,
                },
            }
        )

    return lines


def _nearest_endpoint_station_names(geometry, station_lookup):
    """Identify the two farthest endpoint stations for a MultiLineString."""
    endpoints = []
    for line in geometry.geoms:
        endpoints.append(Point(line.coords[0]))
        endpoints.append(Point(line.coords[-1]))

    farthest_pair = max(
        ((a, b) for a in endpoints for b in endpoints),
        key=lambda pair: pair[0].distance(pair[1]),
    )
    return tuple(_nearest_station_name(point, station_lookup) for point in farthest_pair)


def _nearest_station_name(point, station_lookup):
    """Return the nearest station name to a point."""
    return min(
        station_lookup,
        key=lambda station_name: station_lookup[station_name].distance(point),
    )


def _combine_multiline_geometries(geometries):
    """Combine LineString/MultiLineString geometries into one MultiLineString."""
    lines = []
    for geometry in geometries:
        if isinstance(geometry, LineString):
            lines.append(geometry)
        elif isinstance(geometry, MultiLineString):
            lines.extend(list(geometry.geoms))
        else:
            raise TypeError(f"Unsupported route geometry type: {type(geometry)!r}")
    return MultiLineString(lines)
# endregion
