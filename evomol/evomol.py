"""Main module."""

import evomol.evaluation as evaluator
from evomol.action import molecular_graph as mg
from evomol.action import pprint_action_space
from evomol.representation import SMILES, MolecularGraph, Molecule


def main() -> None:
    """Main function."""
    Molecule.id_representation_class = SMILES
    Molecule.representations_class = [MolecularGraph]
    Molecule.max_heavy_atoms = 38

    Molecule.accepted_atoms = ["C", "O", "N", "F", "S"]

    neighborhood_tester()

    evaluation_tester()


def neighborhood_tester() -> None:
    """Test the neighborhood of atoms in a molecule."""

    mg.ChangeBondMolGraph.avoid_bond_breaking = False
    mg.RemoveGroupMolGraph.remove_only_smallest = False
    MolecularGraph.action_space = [
        mg.AddAtomMolGraph,
        mg.AddGroupMolGraph,
        mg.ChangeBondMolGraph,
        mg.CutAtomMolGraph,
        mg.InsertCarbonMolGraph,
        mg.MoveGroupMolGraph,
        mg.RemoveAtomMolGraph,
        mg.RemoveGroupMolGraph,
        mg.SubstituteAtomMolGraph,
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


def evaluation_tester() -> None:
    """Test the evaluation of a molecule."""

    evaluations = [
        evaluator.CLScore(),
        evaluator.CycleScore(),
        evaluator.GenericCyclicFeatures(),
        evaluator.GenericCyclicScaffolds(),
        evaluator.Isomer("C9H8O4"),
        evaluator.LogP(),
        evaluator.NPerturbations(),
        evaluator.QED(),
        evaluator.RDFilters(),
        evaluator.Rediscovery("C1=CC=CC=C1"),
        evaluator.SAScore(),
        evaluator.SillyWalks(radius=2),
    ]

    # il n'y a pas de fichier json sur https://github.com/PatWalters/silly_walks
    # qu'est ce que c'est "data/ECFP/complete_ChEMBL_ecfp4_dict.json" ?
    smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
    smiles = "CC(=O)NC1=CC=C(C=C1)O"
    smiles = "c1ccccc1"
    smiles = "C=1=CC1"
    # smiles = "C"

    mol = Molecule(smiles)

    for evaluation in evaluations:
        print(evaluation.name, evaluation.evaluate(mol))
