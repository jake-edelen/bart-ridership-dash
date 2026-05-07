"""Reference-table helpers for canonical BART station identity.

The dashboard still has some legacy app-code compatibility paths, but processed
data should expose official `station_id` and `display_name` fields from the
reference tables in `data/reference/`.
"""

# region Imports
from functools import lru_cache
from pathlib import Path

import pandas as pd

from .config import STATION_ALIASES_CSV, STATIONS_MASTER_CSV
# endregion


# region Loaders
@lru_cache(maxsize=None)
def load_station_master(path=None):
    """Load the canonical station dimension table."""
    csv_path = Path(path) if path else STATIONS_MASTER_CSV
    station_master = pd.read_csv(csv_path, dtype={"station_id": "string"})
    station_master["station_id"] = station_master["station_id"].astype(str)
    return station_master


@lru_cache(maxsize=None)
def load_station_aliases(path=None):
    """Load source-specific station-code aliases."""
    csv_path = Path(path) if path else STATION_ALIASES_CSV
    aliases = pd.read_csv(
        csv_path,
        dtype={
            "source_system": "string",
            "raw_code": "string",
            "station_id": "string",
        },
    )
    aliases["source_system"] = aliases["source_system"].astype(str)
    aliases["raw_code"] = aliases["raw_code"].astype(str)
    aliases["station_id"] = aliases["station_id"].astype(str)
    return aliases
# endregion


# region Normalization
def normalize_station_id(raw_code, source_system):
    """Normalize one source-specific station code into canonical `station_id`."""
    aliases = load_station_aliases()
    source_aliases = aliases[aliases["source_system"] == str(source_system)].copy()
    lookup = {
        _normalize_raw_code(row["raw_code"]): row["station_id"]
        for _, row in source_aliases.iterrows()
    }
    key = _normalize_raw_code(raw_code)
    if key not in lookup:
        raise KeyError(f"No station alias for {source_system}:{raw_code}")
    return lookup[key]


def add_station_identity(dataframe, code_column="Entry Station", source_system="legacy_app_code"):
    """Add canonical `station_id` and `display_name` columns to a DataFrame.

    Args:
        dataframe: DataFrame containing source-specific station codes.
        code_column: Column with raw/source station code values.
        source_system: Alias source system to use from `station_aliases.csv`.

    Returns:
        Copy of `dataframe` with `station_id` and `display_name` columns.

    Raises:
        ValueError: If any source code does not map to `station_id` or any
            `station_id` does not exist in `stations_master.csv`.
    """
    result = dataframe.copy()
    aliases = load_station_aliases()
    station_master = load_station_master()

    source_aliases = aliases[aliases["source_system"] == str(source_system)].copy()
    source_aliases["_raw_code_key"] = source_aliases["raw_code"].map(_normalize_raw_code)

    result["_raw_code_key"] = result[code_column].map(_normalize_raw_code)
    result = result.merge(
        source_aliases[["_raw_code_key", "station_id"]],
        on="_raw_code_key",
        how="left",
        validate="many_to_one",
    )

    missing_aliases = sorted(
        result.loc[result["station_id"].isna(), code_column].dropna().astype(str).unique()
    )
    if missing_aliases:
        raise ValueError(
            f"Missing station aliases for source_system={source_system}: {missing_aliases}"
        )

    result = result.merge(
        station_master[["station_id", "display_name"]],
        on="station_id",
        how="left",
        validate="many_to_one",
    )
    missing_station_ids = sorted(
        result.loc[result["display_name"].isna(), "station_id"].dropna().astype(str).unique()
    )
    if missing_station_ids:
        raise ValueError(f"Station IDs missing from stations_master.csv: {missing_station_ids}")

    return result.drop(columns=["_raw_code_key"])


def _normalize_raw_code(raw_code):
    """Return a stable comparison key for source station labels."""
    return str(raw_code).strip().upper()
# endregion
