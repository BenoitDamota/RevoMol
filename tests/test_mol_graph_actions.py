"""
Tests for the molecular graph actions.
"""

from evomol.representation.molecular_graph import MolecularGraph
from evomol.representation.molecular_graph_actions import (
    AddAtomMolGraph,
)

# CutAtomMolGraph,
from evomol.representation.molecule import Molecule
from evomol.representation.smiles import SMILES


def set_parameters() -> None:
    """set the parameters for the tests."""
    Molecule.id_representation_class = SMILES
    Molecule.representations_class = [MolecularGraph]
    Molecule.max_heavy_atoms = 38
    Molecule.accepted_atoms = ["C", "O", "N", "F"]


def test_add_atom():
    """test the add atom action."""
    set_parameters()
    mol = Molecule("")
    action_space = Molecule.get_action_space(MolecularGraph, AddAtomMolGraph)
    actions = action_space.list_actions(mol)
    print(actions)
