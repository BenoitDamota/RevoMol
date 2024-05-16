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
from evomol.representation.molecule import Molecule, pprint_action_space, max_valence
from evomol.representation.smiles import SMILES


def main() -> None:
    """Main function."""
    Molecule.id_representation_class = SMILES
    Molecule.representations_class = [MolecularGraph]
    Molecule.max_heavy_atoms = 38

    accepted_atoms: list[str] = ["C", "O", "N", "F"]
    accepted_substitutions: dict[str, list[str]] = {}

    for accepted_atom in accepted_atoms:
        if accepted_atom != "H":
            curr_atom_subst: list[str] = []
            for other_accepted_atom in accepted_atoms:
                if other_accepted_atom != accepted_atom:
                    curr_atom_subst.append(other_accepted_atom)
            accepted_substitutions[accepted_atom] = curr_atom_subst

    # accepted_substitutions = {"C": ["O", "N"], "O": ["C", "N"], "N": ["C", "O"]}

    MolecularGraph.action_space = [
        ActionSpaceAddAtomMolGraph(accepted_atoms=accepted_atoms),
        ActionSpaceChangeBondMolGraph(avoid_break_bond=False),
        ActionSpaceCutAtomMolGraph(),
        ActionSpaceInsertCarbonMolGraph(),
        ActionSpaceMoveFunctionalGroupMolGraph(),
        ActionSpaceRemoveAtomMolGraph(),
        ActionSpaceRemoveGroupMolGraph(only_remove_smallest_group=False),
        ActionSpaceSubstituteAtomMolGraph(
            accepted_substitutions=accepted_substitutions
        ),
    ]

    # for smiles in ("C", "C(C)(C)(C)(C)", "O=C=O", "N"):
    #     mol = Molecule(smiles)
    #     mol_graph = mol.get_representation(MolecularGraph)
    #     for atom in range(mol_graph.nb_atoms):
    #         print(
    #             mol,
    #             atom,
    #             max_valence(mol_graph.atom_type(atom)),
    #             mol_graph.atom_degree(atom, as_multigraph=True),
    #             mol_graph.atom_degree(atom, as_multigraph=False),
    #         )

    # mol = Molecule("CC(=O)OC1=CC=CC=C1C(=O)O")
    mol = Molecule("O=CC1=NC=CC=C1")
    # mol = Molecule("CCO")
    # mol = Molecule("")
    actions = mol.list_all_possible_actions()
    pprint_action_space(actions)
    for representation, action_space in actions.items():
        for _, actions_list in action_space.items():
            for action in actions_list:
                print(action, "->", action.apply(mol))

    NeighborhoodStrategy.depth = 1
    RandomNeighborhoodStrategy().mutate(mol)

    # # garder voisins d'une molécule dans la mémoire.
