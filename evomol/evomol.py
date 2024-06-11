"""Main module."""

from evomol.representation.molecular_graph import MolecularGraph
from evomol.representation.molecular_graph_actions import (
    AddAtomMolGraph,
    AddGroupMolGraph,
    ChangeBondMolGraph,
    CutAtomMolGraph,
    InsertCarbonMolGraph,
    MoveGroupMolGraph,
    RemoveAtomMolGraph,
    RemoveGroupMolGraph,
    SubstituteAtomMolGraph,
)
from evomol.representation.molecule import Molecule, pprint_action_space
from evomol.representation.smiles import SMILES


def main() -> None:
    """Main function."""
    # neighborhood_tester()
    evaluation_tester()


def neighborhood_tester() -> None:
    """Test the neighborhood of atoms in a molecule."""
    Molecule.id_representation_class = SMILES
    Molecule.representations_class = [MolecularGraph]
    Molecule.max_heavy_atoms = 38

    Molecule.accepted_atoms = ["C", "O", "N", "F", "S"]

    ChangeBondMolGraph.avoid_bond_breaking = False
    RemoveGroupMolGraph.remove_only_smallest = False
    MolecularGraph.action_space = [
        # AddAtomMolGraph,
        AddGroupMolGraph,
        # ChangeBondMolGraph,
        # CutAtomMolGraph,
        # InsertCarbonMolGraph,
        # MoveGroupMolGraph,
        # RemoveAtomMolGraph,
        # RemoveGroupMolGraph,
        # SubstituteAtomMolGraph,
    ]

    smiles = "C[C@H](N)O"
    smiles = "C[C@@H](N)O"
    smiles = "[CH3]O"
    smiles = "[CH2]O"
    smiles = "O=S=O"
    smiles = "CCCC"
    # smiles = "O=C"
    # smiles = "C1=CSC(=C1)C=O"
    # smiles = "C[SH](=O)=O"
    # smiles = "CS(O)O"
    # smiles = "[CH3]O"
    # smiles = "CS(=O)=O"
    # smiles = "C1=CC=C[CH]=C1C1=CC=C[CH]=C1"
    # smiles = "C[CH2]C"
    # # smiles = "[OH2]"
    # smiles = "C[OH]"
    # smiles = "C[NH2]"
    # smiles = "C[SH2]"
    # smiles = "C#S=O"
    # # smiles = "N1C=CC=C1"
    # # smiles = "B=O"
    # smiles = "BrC"
    # smiles = "FS(F)(F)(F)(F)F"
    # smiles = "FS(F)(F)(F)F"
    # smiles = "P(Cl)(Cl)(Cl)(Cl)Cl"
    # smiles = "P(Cl)(Cl)(Cl)Cl"
    # smiles = "C1=CC=C[CH]=C1CC"
    smiles = "C"
    # smiles = "C([H])([H])([H])"
    # smiles = "[CH3]"
    # smiles = "COO"
    # smiles = "N=O"
    # smiles = "C[N+](=O)[O-]"
    # smiles = "C[N+]([O-])O"
    # smiles = "C[NH+]([O-])O"
    # smiles = "[OH-]"
    # smiles = "[N+]"
    # smiles = "O=NO"
    # smiles = "c1ccc(cc1)[N+](=O)[O-]"
    # smiles = "CC(O)=CC"
    # smiles = "CC(=O)CC"
    # smiles = "C1(N)=NC=CC=C1"
    # smiles = "C1(=N)NC=CC=C1"
    # smiles = "S"
    # smiles = "SOO"
    # smiles = "O=S=O"
    # smiles = "O1S(=O)(=O)O1"
    # smiles = "OS(=O)(=O)O"
    # smiles = "OS(=O)(O)O"
    # smiles = "SO"
    # smiles = "OSO"
    # smiles = "O=SO"
    # smiles = "O=S[O+]"
    # smiles = "O=[S+]O"
    # smiles = "O=O"
    # smiles = "C#N"
    # smiles = "C:1:C:C:C:C:C1"
    # smiles = "C1=CC=C[CH]=C1"
    # smiles = "C1=CC=CN=C1"
    # smiles = "[CH3]"
    # smiles = "[CH2]O"
    # smiles = "O[CH2]O"
    # smiles = "O[CH2]OO"
    # smiles = # test move functional grou
    # smiles = "C1=CC=C[CH]=C1C1=CC=C[CH]=C1"
    # smiles = "c1ncncc1CCO"
    # smiles = "CCCO"
    # smiles = "c1ccc[N+]1c1ccoc1",  # 3 links with N for each
    # smiles = "c1ccc[N+]1c1[CH]ccoc1",  # 3 links with N for each C and not [CH
    # smiles = "c1ccc[S+](=O)1c1[CH]ccoc1",  # 3 links with N for each C and not [CH
    # smiles = "c1ccc[S+](=O)1[N+]1[CH]ccocc1",  # no actio
    # smiles = "C=C"
    # smiles = # "[19F][13C@H]([16OH])[35Cl]"
    # smiles = "N[C@@H](C)C(=O)O"
    # smiles = "N[C@H](C)C(=O)O"
    # smiles = "NC(N)C(=O)O"
    # smiles = "COO"
    # smiles = "[OH-]"
    # smiles = "[N+]"
    # smiles = "O=S=O"
    # smiles = "[CH3]"
    # smiles = "C[H]"
    # smiles = "C[C@H](N)O"
    # smiles = "[N+]C"

    mol = Molecule(smiles)
    print(f"{mol=}")

    actions = mol.list_all_possible_actions()
    pprint_action_space(actions)
    for _, action_space in actions.items():
        for _, actions_list in action_space.items():
            for action in actions_list:
                print(action, "->", end=" ")
                new_mol = action.apply()
                print(new_mol)
                # new_mol_graph = new_mol.get_representation(MolecularGraph)

                # # print("atoms", new_mol_graph.atoms)
                # # for atom in new_mol_graph.mol.GetAtoms():
                # #     print(
                # #         atom.GetSymbol(),
                # #         atom.GetFormalCharge(),
                # #         atom.GetImplicitValence(),
                # #         atom.GetNoImplicit(),
                # #         atom.GetNumRadicalElectrons(),
                # #         atom.GetBoolProp("mutability") and not atom.GetNoImplicit(),
                # #     )
                # new_mol_graph = new_mol.get_representation(MolecularGraph)
                # new_mol_graph.draw()
                # print()


def evaluation_tester():
    import evomol.evaluation.evaluation as eval

    evaluations = [
        eval.LogP(),
        eval.SAScore(),
        eval.QED(),
        eval.CycleScore(),
        eval.CLScore(),
        eval.IsomerGuacamol(
            "C9H8O4"
        ),  # "CC(=O)OC1=CC=CC=C1C(=O)O", "CC(=O)NC1=CC=C(C=C1)O" "C"
        eval.Rediscovery("C1=CC=CC=C1"),
        eval.QED(),
        eval.UnkownGenericCyclicScaffolds("path"),
        eval.SillyWalks("path", radius=2),
        eval.NPerturbations(),
        eval.RDFilters(),
    ]

    smiles = "C"

    mol = Molecule(smiles)

    for evaluation in evaluations:
        print(evaluation, evaluation(mol))
