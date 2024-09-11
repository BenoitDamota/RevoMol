"""
Tests for the molecular graph.
"""

import os

import pytest

from evomol import default_parameters as dp
from evomol.representation import Molecule, MolecularGraph


@pytest.fixture(autouse=True, scope="module")
def setup_parameters() -> None:
    """set the parameters for the tests."""
    dp.setup_default_parameters(
        accepted_atoms=["C", "O", "N", "F"],
        max_heavy_atoms=38,
    )


########################################
#           explicit_valences
########################################
@pytest.mark.parametrize(
    ("smiles", "explicit_valences"),
    (
        pytest.param(
            "",
            [],
            id="explicit valence of empty smiles",
        ),
        pytest.param(
            "C1CCCCC1",
            [2, 2, 2, 2, 2, 2],
            id="explicit valence of C1CCCCC1",
        ),
        pytest.param(
            "c1ccccccc1",
            [3, 3, 3, 3, 3, 3, 3, 3],
            id="explicit valence of c1ccccccc1",
        ),
        pytest.param(
            "c1ccccccccc1",
            [3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
            id="explicit valence of c1ccccccccc1",
        ),
    ),
)
def test_molecular_graph_explicit_valence(
    smiles: str, explicit_valences: list[int]
) -> None:
    """Test the explicit valences of the molecular graph."""
    assert MolecularGraph(smiles).explicit_valences == explicit_valences
