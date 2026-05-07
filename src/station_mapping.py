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
    "PITT": "Pittsburg/Bay Point",
}
# endregion


# region Station Code Normalization
HOURLY_OD_STATION_ALIASES = {
    "12TH": "12th",
    "16TH": "16th",
    "19TH": "19th",
    "24TH": "24th",
    "ANTC": "AN",
    "ASHB": "AS",
    "BALB": "BP",
    "BAYF": "BF",
    "BERY": "BE",
    "CAST": "CV",
    "CIVC": "CC",
    "COLM": "CM",
    "COLS": "CL",
    "CONC": "CN",
    "DALY": "DC",
    "DBRK": "BK",
    "DELN": "EN",
    "DUBL": "ED",
    "EMBR": "EM",
    "FRMT": "FM",
    "FTVL": "FV",
    "GLEN": "GP",
    "HAYW": "HY",
    "LAFY": "LF",
    "LAKE": "LM",
    "MCAR": "MA",
    "MLPT": "ML",
    "MLBR": "MB",
    "MONT": "MT",
    "NBRK": "NB",
    "NCON": "NC",
    "OAKL": "OA",
    "ORIN": "OR",
    "PCTR": "PC",
    "PHIL": "PH",
    "PITT": "PITT",
    "PLZA": "EP",
    "POWL": "PL",
    "RICH": "RM",
    "ROCK": "RR",
    "SANL": "SL",
    "SBRN": "SB",
    "SFIA": "SO",
    "SHAY": "SH",
    "SSAN": "SS",
    "UCTY": "UC",
    "WARM": "WD",
    "WCRK": "WC",
    "WDUB": "WP",
    "WOAK": "OW",
}

WORKBOOK_STATION_ALIASES = {
    "WP": "PITT",
    "WD": "WP",
    "WS": "WD",
}

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
    "Coliseum/Airport Connector",
}


def normalize_station_code(code):
    """Return the canonical ridership station code used by service patterns."""
    text_code = str(code).strip()
    return CANONICAL_CODE_ALIASES.get(text_code, text_code)


def normalize_hourly_od_station_code(code):
    """Return the canonical station code for four-character hourly OD labels."""
    text_code = str(code).strip().upper()
    return HOURLY_OD_STATION_ALIASES.get(text_code, normalize_station_code(text_code))


def normalize_workbook_station_code(code):
    """Return the canonical station code for BART ridership workbook labels.

    BART workbook labels use legacy abbreviations for some stations. For
    example, workbook `WP` means Pittsburg/Bay Point, while the app's canonical
    `WP` means West Dublin/Pleasanton.
    """
    text_code = normalize_station_code(code)
    return WORKBOOK_STATION_ALIASES.get(text_code, text_code)
# endregion
