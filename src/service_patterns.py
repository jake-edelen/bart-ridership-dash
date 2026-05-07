"""Date-aware BART service pattern definitions.

Routes are modeled as service patterns over station codes because BART route
identity changes over time while the physical track geometry is reused. These
definitions are intentionally independent from Plotly and GeoPandas so they can
support future map filtering, validation, and modeling workflows.
"""

# region Imports
from dataclasses import dataclass
from datetime import date, datetime
from typing import Iterable

from .station_mapping import STATION_MAPPING, normalize_station_code
# endregion


# region Constants
ANTIOCH_START = date(2018, 5, 26)
BERRYESSA_START = date(2020, 6, 13)
REIMAGINED_SCHEDULE_START = date(2023, 9, 11)

MIN_SERVICE_DATE = date(2018, 1, 1)
OPEN_ENDED = None

BASE_NORTH_SOUTH = [
    "RM",
    "EN",
    "EP",
    "NB",
    "BK",
    "AS",
    "MA",
    "19th",
    "12th",
    "OW",
    "EM",
    "MT",
    "PL",
    "CC",
    "16th",
    "24th",
    "GP",
    "BP",
    "DC",
]

RICHMOND_BRANCH = ["RM", "EN", "EP", "NB", "BK", "AS", "MA"]
SF_PENINSULA = ["OW", "EM", "MT", "PL", "CC", "16th", "24th", "GP", "BP", "DC"]
MILLBRAE_BRANCH = ["CM", "SS", "SB", "SO", "MB"]
MILLBRAE_THEN_SFO = ["CM", "SS", "SB", "MB", "SO"]
DUBLIN_BRANCH = ["ED", "WP", "CV", "BF", "SL", "CL", "FV", "LM"]
FREMONT_BRANCH = ["FM", "UC", "SH", "HY", "BF", "SL", "CL", "FV", "LM"]
BERRYESSA_BRANCH = ["BE", "ML", "WD", "FM", "UC", "SH", "HY", "BF", "SL", "CL", "FV", "LM"]
ANTIOCH_BRANCH = ["AN", "PC", "NC", "CN", "PH", "WC", "LF", "OR", "RR", "MA"]
OAKLAND_NORTH_TRUNK = ["12th", "19th", "MA"]
OAKLAND_SF_TRUNK = ["12th", "OW"]
OAKLAND_AIRPORT_CONNECTOR = ["CL", "OA"]
# endregion


# region Data Model
@dataclass(frozen=True)
class ServicePattern:
    """A date-valid BART service pattern over ordered station codes."""

    service_id: str
    line_color: str
    display_name: str
    valid_from: date
    valid_to: date | None
    station_codes: tuple[str, ...]
    notes: str = ""

    def is_active_on(self, service_date):
        """Return whether this service pattern is active on `service_date`."""
        service_date = parse_service_date(service_date)
        if service_date < self.valid_from:
            return False
        return self.valid_to is None or service_date <= self.valid_to

    def intersects_range(self, start_date, end_date):
        """Return whether this service pattern overlaps a date range."""
        start_date = parse_service_date(start_date)
        end_date = parse_service_date(end_date)
        pattern_end = self.valid_to or date.max
        return self.valid_from <= end_date and pattern_end >= start_date
# endregion


# region Service Pattern Definitions
SERVICE_PATTERNS = tuple(
    ServicePattern(
        service_id="orange_2018_pre_antioch",
        line_color="Orange",
        display_name="Fremont to Richmond",
        valid_from=MIN_SERVICE_DATE,
        valid_to=ANTIOCH_START.replace(day=25),
        station_codes=tuple(FREMONT_BRANCH + OAKLAND_NORTH_TRUNK + list(reversed(RICHMOND_BRANCH[:-1]))),
        notes="Pre-eBART service era before Pittsburg Center/Antioch ridership appears.",
    )
    for _ in [None]
) + (
    ServicePattern(
        service_id="green_2018_pre_antioch",
        line_color="Green",
        display_name="Fremont to Daly City",
        valid_from=MIN_SERVICE_DATE,
        valid_to=ANTIOCH_START.replace(day=25),
        station_codes=tuple(FREMONT_BRANCH + OAKLAND_SF_TRUNK + SF_PENINSULA[1:]),
        notes="Pre-Berryessa Fremont-origin green service.",
    ),
    ServicePattern(
        service_id="blue_all",
        line_color="Blue",
        display_name="Dublin/Pleasanton to Daly City",
        valid_from=MIN_SERVICE_DATE,
        valid_to=OPEN_ENDED,
        station_codes=tuple(DUBLIN_BRANCH + OAKLAND_SF_TRUNK + SF_PENINSULA[1:]),
        notes="Blue line pattern is stable for this project's high-level route model.",
    ),
    ServicePattern(
        service_id="oak_airport_connector_all",
        line_color="Gray",
        display_name="Coliseum to Oakland Airport",
        valid_from=MIN_SERVICE_DATE,
        valid_to=OPEN_ENDED,
        station_codes=tuple(OAKLAND_AIRPORT_CONNECTOR),
        notes="BART to OAK connector service opened before this project's 2018 start date.",
    ),
    ServicePattern(
        service_id="red_pre_2023",
        line_color="Red",
        display_name="Richmond to SFO/Millbrae",
        valid_from=MIN_SERVICE_DATE,
        valid_to=REIMAGINED_SCHEDULE_START.replace(day=10),
        station_codes=tuple(RICHMOND_BRANCH + list(reversed(OAKLAND_NORTH_TRUNK[:-1])) + OAKLAND_SF_TRUNK[1:] + SF_PENINSULA[1:] + MILLBRAE_THEN_SFO),
        notes="Pre-September-2023 red line pattern.",
    ),
    ServicePattern(
        service_id="yellow_2018_antioch_to_sfo",
        line_color="Yellow",
        display_name="Antioch to SFO International Airport",
        valid_from=ANTIOCH_START,
        valid_to=OPEN_ENDED,
        station_codes=tuple(ANTIOCH_BRANCH + ["19th"] + OAKLAND_SF_TRUNK + SF_PENINSULA[1:] + MILLBRAE_BRANCH[:-1]),
        notes="Antioch/Pittsburg Center active from May 26, 2018 onward.",
    ),
    ServicePattern(
        service_id="orange_2018_antioch_to_2020",
        line_color="Orange",
        display_name="Fremont to Richmond",
        valid_from=ANTIOCH_START,
        valid_to=BERRYESSA_START.replace(day=12),
        station_codes=tuple(FREMONT_BRANCH + OAKLAND_NORTH_TRUNK + list(reversed(RICHMOND_BRANCH[:-1]))),
        notes="Post-eBART, pre-Berryessa orange service era.",
    ),
    ServicePattern(
        service_id="green_2018_antioch_to_2020",
        line_color="Green",
        display_name="Fremont to Daly City",
        valid_from=ANTIOCH_START,
        valid_to=BERRYESSA_START.replace(day=12),
        station_codes=tuple(FREMONT_BRANCH + OAKLAND_SF_TRUNK + SF_PENINSULA[1:]),
        notes="Post-eBART, pre-Berryessa green service era.",
    ),
    ServicePattern(
        service_id="orange_2020_plus",
        line_color="Orange",
        display_name="Berryessa/North San José to Richmond",
        valid_from=BERRYESSA_START,
        valid_to=OPEN_ENDED,
        station_codes=tuple(BERRYESSA_BRANCH + OAKLAND_NORTH_TRUNK + list(reversed(RICHMOND_BRANCH[:-1]))),
        notes="Milpitas/Berryessa active from June 13, 2020 onward.",
    ),
    ServicePattern(
        service_id="green_2020_plus",
        line_color="Green",
        display_name="Berryessa/North San José to Daly City",
        valid_from=BERRYESSA_START,
        valid_to=OPEN_ENDED,
        station_codes=tuple(BERRYESSA_BRANCH + OAKLAND_SF_TRUNK + SF_PENINSULA[1:]),
        notes="Milpitas/Berryessa active from June 13, 2020 onward.",
    ),
    ServicePattern(
        service_id="red_2023_plus",
        line_color="Red",
        display_name="Richmond to Millbrae + SFO",
        valid_from=REIMAGINED_SCHEDULE_START,
        valid_to=OPEN_ENDED,
        station_codes=tuple(RICHMOND_BRANCH + list(reversed(OAKLAND_NORTH_TRUNK[:-1])) + OAKLAND_SF_TRUNK[1:] + SF_PENINSULA[1:] + MILLBRAE_BRANCH),
        notes="Post-reimagined-schedule red line service model.",
    ),
)
# endregion


# region Query API
def parse_service_date(value):
    """Parse a date-like value into a `datetime.date`."""
    if isinstance(value, date) and not isinstance(value, datetime):
        return value
    if isinstance(value, datetime):
        return value.date()
    return date.fromisoformat(str(value))


def get_service_patterns_for_date(service_date):
    """Return service patterns active on a single date."""
    service_date = parse_service_date(service_date)
    return tuple(pattern for pattern in SERVICE_PATTERNS if pattern.is_active_on(service_date))


def get_service_patterns_for_range(start_date, end_date, mode="end"):
    """Return service patterns for a selected date range.

    Args:
        start_date: Inclusive range start.
        end_date: Inclusive range end.
        mode: `end`, `start`, or `all_changes`.
    """
    start_date = parse_service_date(start_date)
    end_date = parse_service_date(end_date)
    if start_date > end_date:
        raise ValueError("start_date must be on or before end_date.")

    if mode == "end":
        return get_service_patterns_for_date(end_date)
    if mode == "start":
        return get_service_patterns_for_date(start_date)
    if mode == "all_changes":
        return tuple(pattern for pattern in SERVICE_PATTERNS if pattern.intersects_range(start_date, end_date))
    raise ValueError("mode must be one of: end, start, all_changes.")


def get_active_station_codes_for_range(start_date, end_date):
    """Return canonical station codes used by any pattern active in a range."""
    codes = set()
    for pattern in get_service_patterns_for_range(start_date, end_date, mode="all_changes"):
        codes.update(normalize_station_code(code) for code in pattern.station_codes)
    return tuple(sorted(codes))


def iter_service_station_codes(patterns: Iterable[ServicePattern] = SERVICE_PATTERNS):
    """Yield canonical station codes referenced by service patterns."""
    for pattern in patterns:
        for code in pattern.station_codes:
            yield normalize_station_code(code)
# endregion


# region Validation
def validate_service_pattern_station_codes(patterns: Iterable[ServicePattern] = SERVICE_PATTERNS):
    """Return station codes referenced by patterns but missing from mapping."""
    pattern_codes = set(iter_service_station_codes(patterns))
    return tuple(sorted(code for code in pattern_codes if code not in STATION_MAPPING))
# endregion
