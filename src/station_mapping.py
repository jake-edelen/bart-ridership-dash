"""Station abbreviation mapping for BART ridership workbooks.

The current ridership Excel files use compact station labels while the station
geometry file uses full display names. This mapping is the bridge between those
sources and should be audited when adding new years or station datasets.
"""

# region Station Abbreviation Mapping
STATION_MAPPING = {
    "RM": "Richmond",
    "EN": "El Cerrito del Norte",
    "EP": "El Cerrito Plaza",
    "NB": "North Berkeley",
    "BK": "Downtown Berkeley",
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
    "ML": "Milpitas",
    "BE": "Berryessa/North San José",
    "PC": "Pittsburg Center",
    "AN": "Antioch",
}
# endregion


# region Station Code Normalization
CANONICAL_CODE_ALIASES = {
    "12": "12th",
    "12.0": "12th",
    "16": "16th",
    "16.0": "16th",
    "19": "19th",
    "19.0": "19th",
    "24": "24th",
    "24.0": "24th",
}

IGNORED_STATION_NAMES = {
    "eBART Transfer",
}


def normalize_station_code(code):
    """Return the canonical ridership station code used by service patterns."""
    text_code = str(code).strip()
    return CANONICAL_CODE_ALIASES.get(text_code, text_code)
# endregion
