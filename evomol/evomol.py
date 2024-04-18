"""Main module."""

from evomol.neighborhood_strategy import (
    NeighborhoodStrategy,
    RandomNeighborhoodStrategy,
)
from evomol.representation.molecular_graph import MolecularGraph
from evomol.representation.molecular_graph_actions import (
    ActionSpaceAddAtomMolGraph,
    ActionSpaceChangeBondMolGraph,
    ActionSpaceCutAtomMolGraph,
    ActionSpaceInsertCarbonMolGraph,
    ActionSpaceMoveFunctionalGroupMolGraph,
    ActionSpaceRemoveAtomMolGraph,
    ActionSpaceRemoveGroupMolGraph,
    ActionSpaceSubstituteAtomMolGraph,
)
from evomol.representation.molecule import Molecule, pprint_action_space
from evomol.representation.smiles import SMILES


def main() -> None:
    """Main function."""
    Molecule.id_representation_class = SMILES
    Molecule.representations_class = [MolecularGraph]
    Molecule.max_heavy_atoms = 38

    accepted_atoms: list[str] = ["C", "O", "N"]
    accepted_substitutions: dict[str, list[str]] = {}

    for accepted_atom in accepted_atoms:
        if accepted_atom != "H":
            curr_atom_subst: list[str] = []
            for other_accepted_atom in accepted_atoms:
                if other_accepted_atom != accepted_atom:
                    curr_atom_subst.append(other_accepted_atom)
            accepted_substitutions[accepted_atom] = curr_atom_subst

    accepted_substitutions = {"C": ["O", "N"], "O": ["C"], "N": ["C"]}

    MolecularGraph.action_space = [
        ActionSpaceAddAtomMolGraph(accepted_atoms=accepted_atoms, allow_bonding=True),
        ActionSpaceChangeBondMolGraph(prevent_removing_creating_bonds=False),
        ActionSpaceCutAtomMolGraph(),
        ActionSpaceInsertCarbonMolGraph(),
        ActionSpaceMoveFunctionalGroupMolGraph(),
        ActionSpaceRemoveAtomMolGraph(),
        ActionSpaceRemoveGroupMolGraph(only_remove_smallest_group=False),
        ActionSpaceSubstituteAtomMolGraph(
            accepted_substitutions=accepted_substitutions
        ),
    ]

    mol = Molecule("CC(=O)OC1=CC=CC=C1C(=O)O")
    # mol = Molecule("CCO")
    # mol = Molecule("")
    actions = mol.list_all_possible_actions()
    pprint_action_space(actions)

    NeighborhoodStrategy.depth = 10
    RandomNeighborhoodStrategy().mutate(mol)

    # garder voisins d'une molécule dans la mémoire.
