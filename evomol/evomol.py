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
    MoveGroupMolGraph,
    RemoveAtomMolGraph,
    RemoveGroupMolGraph,
    SubstituteAtomMolGraph,
)
from evomol.representation.molecule import Molecule, pprint_action_space, max_valence
from evomol.representation.smiles import SMILES


def main() -> None:
    neighborhood_tester()
    return

    """Main function."""
    Molecule.id_representation_class = SMILES
    Molecule.representations_class = [MolecularGraph]
    Molecule.max_heavy_atoms = 38

    Molecule.accepted_atoms = ["C", "O", "N", "F", "S"]

    ChangeBondMolGraph.avoid_bond_breaking = False
    RemoveGroupMolGraph.only_remove_smallest_group = False
    SubstituteAtomMolGraph.init_accepted_substitutions_from_accepted_atoms()
    MolecularGraph.action_space = [
        # AddAtomMolGraph,
        # ChangeBondMolGraph,
        # CutAtomMolGraph,
        InsertCarbonMolGraph,
        # MoveGroupMolGraph,
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

    for smiles in [
        # "C",
        # "C([H])([H])([H])",
        # "[CH3]",
        # "COO",
        # "N=O",
        # "C[N+](=O)[O-]",
        # "C[N+]([O-])O",
        # "C[NH+]([O-])O",
        # "[OH-]",
        # "[N+]",
        # "O=NO",
        # "c1ccc(cc1)[N+](=O)[O-]",
        # "CC(O)=CC",
        # "CC(=O)CC",
        # "C1(N)=NC=CC=C1",
        # "C1(=N)NC=CC=C1",
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
        # # test move functional group
        # "C1=CC=C[CH]=C1C1=CC=C[CH]=C1",
        # "c1ncncc1CCO",
        # "CCCO",
        # "c1ccc[N+]1c1ccoc1",  # 3 links with N for each C
        # "c1ccc[N+]1c1[CH]ccoc1",  # 3 links with N for each C and not [CH]
        # "c1ccc[S+](=O)1c1[CH]ccoc1",  # 3 links with N for each C and not [CH]
        # "c1ccc[S+](=O)1[N+]1[CH]ccocc1",  # no action
        # "C=C",
        # # "[19F][13C@H]([16OH])[35Cl]",
        # "N[C@@H](C)C(=O)O",
        # "N[C@H](C)C(=O)O",
        # "NC(N)C(=O)O",
        # tableau latex
        "COO",
        "[OH-]",
        "[N+]",
        "O=S=O",
        "[CH3]",
        "C[H]",
        "C[C@H](N)O",
        "[N+]C",
    ]:
        from rdkit import Chem

        mol = Molecule(smiles)
        print(mol)
        mol_graph = mol.get_representation(MolecularGraph)
        # print("original", mol)
        # mol_graph = mol.get_representation(MolecularGraph)
        # print("convertd", mol_graph.smiles)
        # # mol_graph.draw()
        # Chem.Kekulize(mol_graph.mol)
        # print("kekulize", mol_graph.smiles)
        # # mol_graph.draw()

        # mol = Molecule(smiles)
        # mol_graph = mol.get_representation(MolecularGraph)
        # Chem.Kekulize(mol_graph.mol, clearAromaticFlags=True)
        # print("kekuli 2", mol_graph.smiles)
        # Chem.Kekulize(mol_graph.mol)
        # print("kekuli 3", mol_graph.smiles)
        # print("aromatic", mol_graph.canonical_smiles)
        # # mol_graph.draw()

        # mol = Molecule(smiles)
        # mol_graph = mol.get_representation(MolecularGraph)
        # mol_graph.update_representation()
        # print("updatee ", mol_graph.smiles)
        # # mol_graph.draw()

        # mol = Molecule(smiles)
        # mol_graph = mol.get_representation(MolecularGraph)
        # mol_graph.update_representation(update_property_cache=False)
        # print("updnopc ", mol_graph.smiles)

        # # mol = Molecule(smiles)
        # mol_graph = MolecularGraph(smiles, sanitize=True)
        # mol_graph.update_representation()
        # print("sanitize", mol_graph.smiles)
        # # mol_graph.draw()

        print("atoms", mol_graph.atoms)
        for i, atom in enumerate(mol_graph.mol.GetAtoms()):
            print(
                atom.GetSymbol(),
                max_valence(mol_graph.atom_type(i)),
                atom.GetFormalCharge(),
                mol_graph._explicit_valence(i),
                # atom.GetImplicitValence(),
                mol_graph._implicit_valence(i),
                # atom.GetNoImplicit(),
                atom.GetNumRadicalElectrons(),
                atom.GetBoolProp("mutability") and not atom.GetNoImplicit(),
            )

        print()


def neighborhood_tester():
    Molecule.id_representation_class = SMILES
    Molecule.representations_class = [MolecularGraph]
    Molecule.max_heavy_atoms = 38

    Molecule.accepted_atoms = ["C", "O", "N", "F", "S"]

    ChangeBondMolGraph.avoid_bond_breaking = False
    RemoveGroupMolGraph.only_remove_smallest_group = False
    SubstituteAtomMolGraph.init_accepted_substitutions_from_accepted_atoms()
    MolecularGraph.action_space = [
        AddAtomMolGraph,
        ChangeBondMolGraph,
        CutAtomMolGraph,
        InsertCarbonMolGraph,
        MoveGroupMolGraph,
        RemoveAtomMolGraph,
        RemoveGroupMolGraph,
        SubstituteAtomMolGraph,
    ]

    smiles = "C[C@H](N)O"
    smiles = "C[C@@H](N)O"
    smiles = "[CH3]O"
    smiles = "[CH2]O"
    smiles = "O=S=O"
    smiles = "C[SH](=O)=O"
    smiles = "CS(O)O"
    smiles = "[CH3]O"
    smiles = "CS(=O)=O"
    # smiles = "[OH2]"
    smiles = "C[OH]"
    smiles = "C[NH2]"
    smiles = "C[SH2]"
    smiles = "C#S=O"
    # smiles = "N1C=CC=C1"
    # smiles = "B=O"
    smiles = "BrC"
    smiles = "FS(F)(F)(F)(F)F"
    smiles = "FS(F)(F)(F)F"
    smiles = "P(Cl)(Cl)(Cl)(Cl)Cl"
    smiles = "P(Cl)(Cl)(Cl)Cl"

    mol = Molecule(smiles)
    print(f"{mol=}")
    mol_graph = mol.get_representation(MolecularGraph)
    print("atoms", mol_graph.atoms)
    for atom in mol_graph.mol.GetAtoms():
        print(
            atom.GetSymbol(),
            atom.GetFormalCharge(),
            atom.GetImplicitValence(),
            atom.GetNoImplicit(),
            atom.GetNumRadicalElectrons(),
            atom.GetBoolProp("mutability") and not atom.GetNoImplicit(),
        )
    print(f"A {mol_graph.smiles=}")
    for atom in mol_graph.mol.GetAtoms():
        print(
            atom.GetSymbol(),
            atom.GetFormalCharge(),
            atom.GetImplicitValence(),
            atom.GetNoImplicit(),
            atom.GetNumRadicalElectrons(),
            atom.GetBoolProp("mutability") and not atom.GetNoImplicit(),
        )
    print(f"B {mol_graph.canonical_smiles=}")
    for atom in mol_graph.mol.GetAtoms():
        print(
            atom.GetSymbol(),
            atom.GetFormalCharge(),
            atom.GetImplicitValence(),
            atom.GetNoImplicit(),
            atom.GetNumRadicalElectrons(),
            atom.GetBoolProp("mutability") and not atom.GetNoImplicit(),
        )
    mol_graph.update_representation(True)
    print("UPDATE REPRESENTATION")
    for atom in mol_graph.mol.GetAtoms():
        print(
            atom.GetSymbol(),
            atom.GetFormalCharge(),
            atom.GetImplicitValence(),
            atom.GetNoImplicit(),
            atom.GetNumRadicalElectrons(),
            atom.GetBoolProp("mutability") and not atom.GetNoImplicit(),
        )
    print(f"C {mol_graph.smiles=}")
    print(f"D {mol_graph.canonical_smiles=}")

    return
    actions = mol.list_all_possible_actions()
    pprint_action_space(actions)
    for _, action_space in actions.items():
        for _, actions_list in action_space.items():
            for action in actions_list:
                print(action, "->", end=" ")
                new_mol = action.apply()
                print(new_mol)
                new_mol_graph = new_mol.get_representation(MolecularGraph)

                print("atoms", new_mol_graph.atoms)
                for atom in new_mol_graph.mol.GetAtoms():
                    print(
                        atom.GetSymbol(),
                        atom.GetFormalCharge(),
                        atom.GetImplicitValence(),
                        atom.GetNoImplicit(),
                        atom.GetNumRadicalElectrons(),
                        atom.GetBoolProp("mutability") and not atom.GetNoImplicit(),
                    )
                new_mol_graph = new_mol.get_representation(MolecularGraph)
                # new_mol_graph.draw()
                print()
