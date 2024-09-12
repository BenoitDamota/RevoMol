"""
Tests for the molecular graph.
"""

import pytest

from evomol import default_parameters as dp
from evomol.representation import MolecularGraph


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
            id="explicit valences of empty smiles",
        ),
        pytest.param(
            "C1CCCCC1",
            [2, 2, 2, 2, 2, 2],
            id="explicit valences of C1CCCCC1",
        ),
        pytest.param(
            "c1ccccccc1",
            [3, 3, 3, 3, 3, 3, 3, 3],
            id="explicit valences of c1ccccccc1",
        ),
        pytest.param(
            "c1ccccccccc1",
            [3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
            id="explicit valences of c1ccccccccc1",
        ),
    ),
)
def test_molecular_graph_explicit_valence(
    smiles: str, explicit_valences: list[int]
) -> None:
    """Test the explicit valences of the molecular graph."""
    assert MolecularGraph(smiles).explicit_valences == explicit_valences


########################################
#             formal_charges
########################################
@pytest.mark.parametrize(
    ("smiles", "formal_charges"),
    (
        pytest.param(
            "",
            [],
            id="formal charges of empty smiles",
        ),
        pytest.param(
            "C1CCCCC1",
            [0, 0, 0, 0, 0, 0],
            id="formal charges of C1CCCCC1",
        ),
        pytest.param(
            "c1ccccccc1",
            [0, 0, 0, 0, 0, 0, 0, 0],
            id="formal charges of c1ccccccc1",
        ),
        pytest.param(
            "c1ccccccccc1",
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            id="formal charges of c1ccccccccc1",
        ),
        pytest.param(
            "C[N+](=O)[O-]",
            [0, 1, 0, -1],
            id="formal charges of C[N+](=O)[O-]",
        ),
    ),
)
def test_molecular_graph_formal_charges(smiles: str, formal_charges: list[int]) -> None:
    """Test the formal charges of the molecular graph."""
    assert MolecularGraph(smiles).formal_charges == formal_charges


########################################
#          implicit_valences
########################################
@pytest.mark.parametrize(
    ("smiles", "implicit_valences"),
    (
        pytest.param(
            "",
            [],
            id="implicit valences of empty smiles",
        ),
        pytest.param(
            "C1CCCCC1",
            [2, 2, 2, 2, 2, 2],
            id="implicit valences of C1CCCCC1",
        ),
        pytest.param(
            "c1ccccccc1",
            [1, 1, 1, 1, 1, 1, 1, 1],
            id="implicit valences of c1ccccccc1",
        ),
        pytest.param(
            "c1ccccccccc1",
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            id="implicit valences of c1ccccccccc1",
        ),
        pytest.param(
            "C[N+](=O)[O-]",
            [3, 0, 0, 0],
            id="implicit valences of C[N+](=O)[O-]",
        ),
    ),
)
def test_molecular_graph_implicit_valences(
    smiles: str, implicit_valences: list[int]
) -> None:
    """Test the implicit valences of the molecular graph."""
    assert MolecularGraph(smiles).implicit_valences == implicit_valences
