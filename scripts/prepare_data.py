"""Export raw geospatial inputs to `data/processed/`.

This script is intentionally separate from `app.py` so running or importing the
Dash app does not mutate files on disk. Use it when processed GeoJSON/CSV
artifacts need to be regenerated from raw data.
"""

# region Imports And Path Setup
from pathlib import Path
import sys

import geopandas as gpd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.config import BART_ROUTES_GDB, BART_STATIONS_DIR, PROCESSED  # noqa: E402
from src.data_loader import (  # noqa: E402
    write_hourly_station_monthly_summaries,
    write_hourly_workbook_validation,
    write_monthly_station_ridership_summaries,
)
# endregion


# region Script Entrypoint
def main():
    """Regenerate processed station and route exports from raw geospatial files."""
    PROCESSED.mkdir(parents=True, exist_ok=True)

    stations_gdf = gpd.read_file(BART_STATIONS_DIR / "doc.kml", driver="KML")
    stations_gdf.to_file(PROCESSED / "stations.geojson", driver="GeoJSON")
    stations_gdf.to_csv(PROCESSED / "stations.csv", index=False)

    route_to_gdf = gpd.read_file(BART_ROUTES_GDB, layer="BARTLine")
    route_to_gdf.to_file(PROCESSED / "bart_routes.geojson", driver="GeoJSON")

    routes_gdf = gpd.read_file(BART_STATIONS_DIR / "bart_geojson.json")
    routes_gdf.to_file(PROCESSED / "bart_geojson.geojson", driver="GeoJSON")

    hourly_summary_path = write_hourly_station_monthly_summaries()
    print(f"Hourly OD monthly station summaries written to {hourly_summary_path}")

    ridership_summary_path = write_monthly_station_ridership_summaries(
        hourly_summary_path=hourly_summary_path,
    )
    print(f"Monthly ridership summaries written to {ridership_summary_path}")

    validation_path = write_hourly_workbook_validation()
    print(f"Hourly OD validation written to {validation_path}")

    print(f"Processed data written to {PROCESSED}")


if __name__ == "__main__":
    main()
# endregion
