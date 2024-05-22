"""Main module."""

from evomol.neighborhood_strategy import (
    NeighborhoodStrategy,
    RandomNeighborhoodStrategy,
)
from evomol.representation.molecular_graph import MolecularGraph
from evomol.representation.molecular_graph_actions import (
    AddAtomMolGraph,
    ChangeBondMolGraph,
    CutAtomMolGraph,
    InsertCarbonMolGraph,
    MoveFunctionalGroupMolGraph,
    RemoveAtomMolGraph,
    RemoveGroupMolGraph,
    SubstituteAtomMolGraph,
)
from evomol.representation.molecule import Molecule, pprint_action_space
from evomol.representation.smiles import SMILES


def main() -> None:
    """Main function."""
    Molecule.id_representation_class = SMILES
    Molecule.representations_class = [MolecularGraph]
    Molecule.max_heavy_atoms = 38

    Molecule.accepted_atoms = ["C", "O", "N", "F"]
    accepted_substitutions: dict[str, list[str]] = {}

    for accepted_atom in Molecule.accepted_atoms:
        if accepted_atom != "H":
            curr_atom_subst: list[str] = []
            for other_accepted_atom in Molecule.accepted_atoms:
                if other_accepted_atom != accepted_atom:
                    curr_atom_subst.append(other_accepted_atom)
            accepted_substitutions[accepted_atom] = curr_atom_subst

    # accepted_substitutions = {"C": ["O", "N"], "O": ["C", "N"], "N": ["C", "O"]}

    ChangeBondMolGraph.avoid_bond_breaking = False
    RemoveGroupMolGraph.only_remove_smallest_group = False
    SubstituteAtomMolGraph.accepted_substitutions = accepted_substitutions
    MolecularGraph.action_space = [
        # AddAtomMolGraph,
        # ChangeBondMolGraph,
        CutAtomMolGraph,
        # InsertCarbonMolGraph,
        # MoveFunctionalGroupMolGraph,
        # RemoveAtomMolGraph,
        # RemoveGroupMolGraph,
        # SubstituteAtomMolGraph,
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
    # mol = Molecule("O=CC1=NC=CC=C1")
    # mol = Molecule("CCO")
    # mol = Molecule("")

    for smiles in [
        # "C",
        # "C([H])",
        # "C([H])([H])",
        # "C([H])([H])([H])",
        # "C([H])([H])([H])([H])",
        # "COO",
        # "N",
        # "NO",
        # "N=O",
        # "C[N+](=O)[O-]",
        # "C[N+]([O-])O",
        # "C[NH+]([O-])O",
        # "[OH-]",
        # "[N+]",
        # "N(=O)O",
        # "c1ccc(cc1)[N+](=O)[O-]",
        # "S",
        # "SOO",
        # "O=S=O",
        # "O1S(=O)(=O)O1",
        # "OS(=O)(=O)O",
        # "OS(=O)(O)O",
        # "SO",
        # "OSO",
        # "O=SO",
        # "O=S[O+]",
        # "O=[S+]O",
        # "O=O",
        # "C#N",
        # "C:1:C:C:C:C:C1",
        # "C1=CC=C[CH]=C1",
        # "C1=CC=CN=C1",
        # "[CH3]",
        # "[CH2]O",
        # "O[CH2]O",
        # "O[CH2]OO",
        # tests cut atom
        "CCC",  # cut middle C
        "C=C=C",  # cut middle C then simple bond
        "[C+]#CC",  # cut middle C then triple bond
        "O[CH2]O",  # no cut (C is not mutable)
        "O[CH2]OO",  # one cut on 2nd O
        "O[CH]N=O",  # one cut on N
        "O[C+]=C[C+]=O",  # not cut (C+ are not mutable and not the same bond)
        "O[C+]=C=[C+]=O",  # one cut middle C with a double bond after
        "O[C+]=CC",  # one cut on middle C with a double bond after
        "O[C+]C=C",  # one cut on middle C with a simple bond after
        "OC=C[C+]",  # one cut on middle C with a simple bond after and one cut 1st C
    ]:
        mol = Molecule(smiles)
        print(mol)
        # mol_graph = mol.get_representation(MolecularGraph)
        # mol_graph.update_representation()
        # print(mol_graph.canonical_smiles)
        # print("atoms", mol_graph.atoms())
        # for atom in mol_graph.mol.GetAtoms():
        #     print(
        #         atom.GetSymbol(),
        #         atom.GetFormalCharge(),
        #         atom.GetImplicitValence(),
        #         atom.GetNoImplicit(),
        #         atom.GetNumRadicalElectrons(),
        #         atom.GetBoolProp("mutability") and not atom.GetNoImplicit(),
        #     )
        # for bond in mol_graph.mol.GetBonds():
        #     print(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondType())
        # print("formal_charge", mol_graph.formal_charge_vector())
        # print(
        #     "explicit valence",
        #     list(mol_graph._explicit_valence(i) for i in range(mol_graph.nb_atoms)),
        # )
        # print("implicit valence", mol_graph.implicit_valence_vector())
        # print()

        actions = mol.list_all_possible_actions()
        pprint_action_space(actions)
        for _, action_space in actions.items():
            for _, actions_list in action_space.items():
                for action in actions_list:
                    print(action, "->", end=" ")
                    new_mol = action.apply()
                    print(new_mol)
                    new_mol_graph = new_mol.get_representation(MolecularGraph)

        print()

    # NeighborhoodStrategy.depth = 1
    # RandomNeighborhoodStrategy().mutate(mol)

    # # garder voisins d'une molécule dans la mémoire.
