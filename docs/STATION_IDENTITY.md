# Station Identity

This project separates station identity from source-specific station labels.
The dashboard should not expose temporal aliases or raw source codes to users.

## Policy

Use one stable station ID per real-world ridership station:

```text
Internal identity: official BART station abbreviation
Input cleanup: source-specific alias mapping
UI display: readable station name
Timeline behavior: date-aware station availability
```

The app should display station names, not raw codes. If a historical or raw
source uses a different code for the same physical station, that code belongs
in an alias table and should normalize to the same canonical `station_id`.

## Reference Tables

`data/reference/stations_master.csv` is the station dimension table. Its
columns are:

```text
station_id
display_name
latitude
longitude
opened_date
closed_date
geom_source_name
```

`station_id` uses lowercase BART station abbreviations such as `embr`, `warm`,
`wdub`, `mlpt`, and `bery`. These are stable internal IDs, not display labels.

`opened_date` and `closed_date` control whether a station should be available
for a selected time period. Blank date values currently mean the station is
treated as open for the app's 2018-2025 analysis window.

`geom_source_name` preserves the station name from the geometry source so the
reference table can be joined back to station point geometry without relying on
display text alone.

## Alias Handling

Raw data sources should normalize into `station_id` before aggregation.
Examples:

```text
Hourly OD WARM -> warm
Workbook WS -> warm
Hourly OD WDUB -> wdub
Workbook WD -> wdub
Deprecated ASBY -> ashb
Current ASHB -> ashb
```

The dashboard should not show different station codes for different years.
Historical code differences are data-cleaning details, not UI concepts.

## Exclusions

`eBART Transfer` and `Coliseum/Airport Connector` are not included in
`stations_master.csv` because they are not standalone ridership stations in the
current station-ridership dashboard model. They may still appear in geometry
sources and can be handled separately if a future feature needs transfer or
facility-level infrastructure points.
