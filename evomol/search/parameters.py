"""
Module for global search parameters.
"""

from dataclasses import dataclass


@dataclass
class SearchParameters:
    """Global search parameters."""

    fitness_criteria: str = "QED"


search_parameters = SearchParameters()
