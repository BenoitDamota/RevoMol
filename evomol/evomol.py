"""Main module."""

import time

from evomol import evaluation as evaluator
from evomol.action import molecular_graph as mg
from evomol.action import pprint_action_space
from evomol.representation import SMILES, MolecularGraph, Molecule
from evomol.search.evolutionary_algorithm import (
    EvolutionaryAlgorithm,
    ParametersEvolutionaryAlgorithm,
)
from evomol.search.neighborhood_strategy import RandomNeighborhoodStrategy
from evomol.search.parameters import search_parameters


def main() -> None:
    """Main function."""
    Molecule.id_representation_class = SMILES
    Molecule.representations_class = [MolecularGraph]
    Molecule.max_heavy_atoms = 38

    Molecule.accepted_atoms = ["C", "O", "N", "F", "S", "P", "Cl", "Br"]
    Molecule.accepted_atoms = ["C", "O", "N", "F", "S"]

    # neighborhood_tester()

    # evaluation_tester()
    # test_speed()
    # evomol_tester()

    # diversity_tester()

    # neighbors_validity_tester()
    # from_c_change_heavy_atom()


def evomol_tester() -> None:
    """Test the evolutionary algorithm."""

    MolecularGraph.action_space = [
        mg.AddAtomMG,
        mg.AddGroupMG,
        mg.ChangeBondMG,
        mg.CutAtomMG,
        mg.InsertCarbonMG,
        mg.MoveGroupMG,
        mg.RemoveAtomMG,
        mg.RemoveGroupMG,
        mg.SubstituteAtomMG,
    ]

    search_parameters.fitness_criteria = "QED"

    evaluations = [
        evaluator.UnknownGCF(),
        evaluator.FilterUnknownGCF(threshold=0),
        evaluator.UnknownECFP(),
        evaluator.FilterUnknownECFP(threshold=0),
        evaluator.QED,
    ]

    molecules = [Molecule("")]

    population = [
        Molecule(molecule.get_representation(MolecularGraph).aromatic_canonical_smiles)
        for molecule in molecules
    ]

    for molecule in population:
        for evaluation in evaluations:
            evaluation.evaluate(molecule)

    evo_mol = EvolutionaryAlgorithm(
        population=population,
        parameters=ParametersEvolutionaryAlgorithm(),
        evaluations=evaluations,
        neighborhood_strategy=RandomNeighborhoodStrategy(),
    )
    evo_mol.run()


def diversity_tester() -> None:
    """Test the diversity of a molecule."""

    # TODO pour la taille voir publi Scalable estimator of the diversity ou
    # faire sans limite mais mémoire possiblement pas contigue
    descriptors = [
        # evaluator.Descriptor(  # TODO fonctionne pas
        #     "scaffolds",
        #     evaluator.scaffolds,
        #     1000,
        # ),
        evaluator.Descriptor(
            "gen_scaffolds",
            evaluator.gen_scaffolds,
            1000,
        ),
        evaluator.Descriptor(
            "compute_ifg",
            evaluator.compute_ifg,
            1000,
        ),
        evaluator.Descriptor(
            "atoms_list",
            evaluator.atoms_list,
            20,
        ),
        evaluator.Descriptor(
            "shg_1",
            evaluator.shg_1,
            1000,
        ),
        # evaluator.Descriptor(
        #     "checkmol",
        #     evaluator.checkmol,
        #     1000,
        # ),
        evaluator.Descriptor(
            "ECFP4",
            evaluator.ecfp4,
            1024,
        ),
    ]

    molecules = [
        Molecule("O=S(=O)(O)c1ccccc1-c1ccsn1"),
        Molecule("COS(=O)(=O)Nc1ccccc1"),
        Molecule("Oc1cnccc1-c1cccs1"),
        Molecule("Cc1cccc(CS(=O)(=O)O)c1"),
        Molecule("Oc1cccc(-c2cccs2)c1"),
        Molecule("NS(=O)(=O)Oc1ccccc1"),
        Molecule("Cc1ccc(F)cc1S(=O)(=O)O"),
        Molecule("O=S(=O)(O)Cc1ccccc1"),
        Molecule("Fc1ccc(-c2ccc[nH]2)s1"),
        Molecule("FCc1ccc(-c2ccccc2)cc1"),
        Molecule("Sc1cccc(-c2ccco2)c1"),
        Molecule("Sc1ccc(-c2cccs2)cc1"),
        Molecule("Sc1cccc(-c2cccs2)c1"),
        Molecule("O=S(=O)(O)Nc1ccccc1"),
        Molecule("Nc1cccc(-c2ccccc2)c1"),
        Molecule("Sc1cccc(-c2ccccc2)c1"),
        Molecule("c1c[nH]c(-c2cccs2)c1"),
        Molecule("c1csc(-c2cncs2)c1"),
        Molecule("c1ccc(-c2ccsn2)cc1"),
        Molecule("c1ccc(-c2cnsc2)cc1"),
        Molecule("c1cncc(-c2cccs2)c1"),
        Molecule("c1ccc(-c2ccns2)cc1"),
        Molecule("Fc1ccsc1-c1ccccc1"),
        Molecule("CS(=O)(=O)Oc1ccccc1"),
        Molecule("c1csc(-c2ncco2)c1"),
        Molecule("O=S(=O)(O)c1cc2ccc1-2"),
        Molecule("O=S(=O)(O)c1ccccc1"),
        Molecule("Fc1ccccc1-c1ccccc1"),
        Molecule("c1ccc(-c2ccncc2)cc1"),
        Molecule("c1ccc(Nc2ccc3cc2-3)cc1"),
        Molecule("CS(=O)(=O)c1ccccc1"),
        Molecule("c1csc(-c2ccsc2)c1"),
        Molecule("c1csc(-c2cccs2)c1"),
        Molecule("c1csc(-c2ccsc2)c1"),
        Molecule("c1csc(-c2cccs2)c1"),
        Molecule("c1csc(-c2ccsc2)c1"),
        Molecule("c1csc(-c2cccs2)c1"),
        Molecule("c1csc(-c2ccsc2)c1"),
        Molecule("c1csc(-c2cccs2)c1"),
        Molecule("c1csc(-c2ccsc2)c1"),
        Molecule("c1csc(-c2cccs2)c1"),
        Molecule("c1csc(-c2ccsc2)c1"),
        Molecule("C1=CCCC(c2ccccc2)C=C1"),
        Molecule("O=S(O)c1ccccc1"),
        Molecule("c1ccc(-c2ccsc2)cc1"),
        Molecule("c1ccc(-c2cccs2)cc1"),
        Molecule("c1ccc(-c2ccsc2)cc1"),
        Molecule("c1ccc(-c2cccs2)cc1"),
        Molecule("c1ccc(-c2ccsc2)cc1"),
        Molecule("c1ccc(-c2cccs2)cc1"),
        Molecule("c1ccc(-c2ccsc2)cc1"),
        Molecule("c1ccc(-c2ccccc2)cc1"),
        Molecule("C1=CCCC(c2ccccc2)=C1"),
        Molecule("c1csc(-c2ccc(-c3cccs3)cc2)c1"),
        Molecule("NCc1cccs1"),
        Molecule("O=[N+]([O-])c1c[nH]c(-c2ccccc2)c1"),
        Molecule("CNc1ccccc1"),
        Molecule("OCc1ccccc1"),
        Molecule("C=CC=CCc1ccccc1"),
        Molecule("NCc1ccccc1"),
        Molecule("O=[N+]([O-])c1ccc(-c2ccccc2)cc1"),
        Molecule("O=[N+]([O-])c1ccccc1-c1ccccc1"),
        Molecule("SCc1ccccc1"),
        Molecule("Oc1ccc2c(c1)S2"),
        Molecule("Oc1cccc(F)c1"),
        Molecule("Cc1cccc(O)c1"),
        Molecule("COc1ccccc1"),
        Molecule("CCOS(N)(=O)=O"),
        Molecule("Nc1nccs1"),
        Molecule("Oc1cc2ccc1-2"),
        Molecule("FCc1ccccc1"),
        Molecule("c1csc(-c2sc3cc2-3)c1"),
        Molecule("Oc1ccccn1"),
        Molecule("Oc1ccncc1"),
        Molecule("Oc1cccnc1"),
        Molecule("CCc1ccccc1"),
        Molecule("CCCCN"),
        Molecule("Fc1ccc(S)cc1"),
        Molecule("Nc1cc2c(F)cc1-2"),
        Molecule("Oc1ccc(S)cc1"),
        Molecule("Cc1ccc(S)cc1"),
        Molecule("Nc1ccc(F)cc1"),
        Molecule("Nc1cccc(F)c1"),
        Molecule("O=[N+]([O-])c1ccc(-c2cccc(-c3cccs3)c2)cc1"),
        Molecule("Sc1ccccn1"),
        Molecule("Sc1ccncc1"),
        Molecule("O=[N+]([O-])c1ccc(F)c(O)c1"),
        Molecule("c1ccc2c(c1)N2"),
        Molecule("Sc1ccc2cc1-2"),
        Molecule("Sc1cc2cc-2c1"),
        Molecule("Sc1ccc2cc1-2"),
        Molecule("c1cc2cc(c1)OC2"),
        Molecule("Nc1ccc2cc1-2"),
        Molecule("Nc1cc2ccc1-2"),
        Molecule("Sc1ccccc1"),
        Molecule("Cc1c(F)cc2cc1-2"),
        Molecule("Nc1ccccc1"),
        Molecule("Cc1ccc(F)cc1"),
        Molecule("C=CC=CCCC"),
        Molecule("c1c2n[nH]c1-2"),
        Molecule("c1[nH]c2cc1-2"),
        Molecule("c1cn[nH]c1"),
        Molecule("Fc1cc2cc-2c1"),
        Molecule("Fc1cc2ccc1-2"),
        Molecule("Fc1ccc2cc1-2"),
        Molecule("Fc1cc2cc-2c1"),
        Molecule("Fc1cc2ccc1-2"),
        Molecule("c1cc[nH]c1"),
        Molecule("c1ccc2c(c1)S2"),
        Molecule("c1sc2nc1-2"),
        Molecule("CCCO"),
        Molecule("COS(N)(=O)=O"),
        Molecule("CCCN"),
        Molecule("CC=CC=CCC"),
        Molecule("Fc1ccccc1"),
        Molecule("c1cnsc1"),
        Molecule("c1cscn1"),
        Molecule("c1cnsc1"),
        Molecule("c1cscn1"),
        Molecule("c1cnsc1"),
        Molecule("c1cscn1"),
        Molecule("c1cnsc1"),
        Molecule("c1cscn1"),
        Molecule("c1cnsc1"),
        Molecule("c1cscn1"),
        Molecule("c1cnsc1"),
        Molecule("Cc1ccccc1"),
        Molecule("c1cc2nc-2c1"),
        Molecule("c1cc2cc-2n1"),
        Molecule("c1ncc2cc1-2"),
        Molecule("c1sc2cc1-2"),
        Molecule("c1cc2sc1-2"),
        Molecule("c1sc2cc1-2"),
        Molecule("c1ccncc1"),
        Molecule("CC=CCCC"),
        Molecule("c1ccsc1"),
        Molecule("C=CC=CCC"),
        Molecule("c1cc2cc-2c1"),
        Molecule("c1cc2ccc1-2"),
        Molecule("c1cc2cc-2c1"),
        Molecule("c1cocn1"),
        Molecule("c1cnoc1"),
        Molecule("C1=CCCCCC=C1"),
        Molecule("CC1=CCCCC=C1"),
        Molecule("c1ccoc1"),
        Molecule("C=CC=CC=C"),
        Molecule("c1ccccc1"),
        Molecule("O=[N+]([O-])c1ccccc1F"),
        Molecule("O=[N+]([O-])c1cccc(F)c1"),
        Molecule("O=[N+]([O-])c1ccc(F)cc1"),
        Molecule("Cc1cccc([N+](=O)[O-])c1"),
        Molecule("C=CCCC"),
        Molecule("C1=CN=CCCC1"),
        Molecule("O=[N+]([O-])c1cc[nH]c1"),
        Molecule("O=[N+]([O-])c1ccc[nH]1"),
        Molecule("O=[N+]([O-])c1cc[nH]c1"),
        Molecule("C1=CCCCC=C1"),
        Molecule("CCCC"),
        Molecule("C1=COCCC1"),
        Molecule("C1=CC=CCC=C1"),
        Molecule("O=S(=O)(O)NS(=O)(=O)O"),
        Molecule("C1=CC2=CCC2=C1"),
        Molecule("CC=CC=CC"),
        Molecule("O=[N+]([O-])c1ccccc1"),
        Molecule("O=S(=O)(O)CO"),
        Molecule("O=S(O)CO"),
        Molecule("C1=CCCC=C1"),
        Molecule("CC1=CCC=C1"),
        Molecule("C1=CCCC=C1"),
        Molecule("NCS(=O)(=O)O"),
        Molecule("CC1C=CC=C1"),
        Molecule("CS(N)(=O)=O"),
        Molecule("C1=CNC1"),
        Molecule("CS(=O)(=O)O"),
        Molecule("CS(C)(=O)=O"),
        Molecule("[N-]=[N+]=Nc1ccccc1-c1ccncc1"),
        Molecule("CCO"),
        Molecule("C1=CCCCCC1"),
        Molecule("CCN"),
        Molecule("C1=CCC=C1"),
        Molecule("C1=NCN1"),
        Molecule("Nc1ccc(S)cc1"),
        Molecule("Nc1ccccc1S"),
        Molecule("CN1C=CC1"),
        Molecule("CS(C)=O"),
        Molecule("N"),
        Molecule("C=CC=C"),
        Molecule("O=S(O)O"),
        Molecule("C1=COC1"),
        Molecule("C1=CCCCC1"),
        Molecule("CC1C=CCC1"),
        Molecule("C1=CCC=CC1"),
        Molecule("C1=CC=C1"),
        Molecule("CCC"),
        Molecule("CO"),
        Molecule("CN"),
        Molecule("FCF"),
        Molecule("O=[N+]([O-])c1ccc(S)cc1"),
        Molecule("C1CC1"),
        Molecule("O=CO"),
        Molecule("c1c2c3cc-3c1-2"),
        Molecule("C1=CCCC1"),
        Molecule("CC1C=CC1"),
        Molecule("Nc1ccc(S)c(N)c1"),
        Molecule("C1=CC=CC=C1 "),
        Molecule("CC"),
        Molecule("C=CCCCC=C"),
        Molecule("NS(=O)(=O)O"),
        Molecule("C#N"),
    ]

    population_size = 100
    population = molecules[:population_size]

    diversities = [
        evaluator.Diversity(
            descriptor=descriptor,
            size=population_size,
        )
        for descriptor in descriptors
    ]

    for molecule in population:
        for diversity in diversities:
            diversity.evaluate_molecule(molecule)
            diversity.add_molecule(molecule)

    for diversity in diversities:
        diversity.update_entropy()

    print(diversities)

    for smiles in ["C", "CO", "CONF", "P", "PCl"]:
        mol = Molecule(smiles)
        for diversity in diversities:
            diversity.evaluate_molecule(mol)

        for diversity in diversities:
            diversity.update_entropy()


def neighborhood_tester() -> None:
    """Test the neighborhood of atoms in a molecule."""

    mg.ChangeBondMG.avoid_bond_breaking = False
    mg.RemoveGroupMG.remove_only_smallest = False
    MolecularGraph.action_space = [
        mg.AddAtomMG,
        mg.AddGroupMG,
        mg.ChangeBondMG,
        mg.CutAtomMG,
        mg.InsertCarbonMG,
        mg.MoveGroupMG,
        mg.RemoveAtomMG,
        mg.RemoveGroupMG,
        mg.SubstituteAtomMG,
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

    # import statistics

    # def gaussian(x: float, mu: float, sigma: float) -> float:
    #     # calcul d'une gaussienne
    #     pass

    # def my_fitness(molecule: Molecule) -> float:
    #     return statistics.mean(
    #         -(
    #             molecule.value("LogP") * 0.5
    #             + gaussian(molecule.value("SAScore"), 3, 1) * 0.5
    #         ),
    #         molecule.value("QED"),
    #     )

    # eval_strat = evaluator.Function("my_fitness", my_fitness)

    evaluations = [
        evaluator.CLScore(),
        evaluator.UnknownGCF(),
        evaluator.Isomer("C9H8O4"),
        evaluator.LogP,
        evaluator.NPerturbations,
        evaluator.QED,
        evaluator.RDFilters(),
        evaluator.Rediscovery("C1=CC=CC=C1"),
        evaluator.SAScore,
        evaluator.FilterUnknownECFP(),
        # eval_strat,
    ]

    smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
    smiles = "CC(=O)NC1=CC=C(C=C1)O"
    smiles = "c1ccccc1"
    smiles = "C=1=CC1"
    # smiles = "C"

    mol = Molecule(smiles)

    eval_: evaluator.Evaluation
    for eval_ in evaluations:
        print(eval_.name, eval_.evaluate(mol))
    print()


def neighbors_validity_tester():
    # smiles = "C1=CSC(=C2SC=CS2)S1"  # TTF
    # smiles = "N1=S=NC2=C1N=S=N2" # DD
    Molecule.max_heavy_atoms = 10
    # mol = Molecule(smiles)

    MolecularGraph.action_space = [
        mg.AddAtomMG,
        # mg.AddGroupMG,
        mg.ChangeBondMG,
        mg.CutAtomMG,
        mg.InsertCarbonMG,
        mg.MoveGroupMG,
        mg.RemoveAtomMG,
        # mg.RemoveGroupMG,
        mg.SubstituteAtomMG,
    ]

    evaluations = [
        evaluator.UnknownGCF(),
        evaluator.FilterUnknownGCF(threshold=0),
        evaluator.UnknownECFP(radius=2),
        evaluator.FilterUnknownECFP(threshold=0),
    ]

    print(
        "molecule,depth,nb_molecules,|molecules|"
        ",nb_valids,|valid|,filter,time,time post filter"
    )
    for smiles in [
        "C"
        # "C1=CSC(=C2SC=CS2)S1",  # TTF
        # "N1=S=NC2=C1N=S=N2",  # DD
        # "C1=CC(=CC=C1C2=C3C=CC(=C(C4=NC(=C(C5=CC=C(N5)C(=C6C=CC2=N6)C7=CC=C
        # (C=C7)C(=O)O)C8=CC=C(C=C8)C(=O)O)C=C4)C9=CC=C(C=C9)C(=O)O)N3)C(=O)O",
        # # porph
    ]:
        for apply_evaluation in [True, False]:
            for depth in range(1, 3):
                start = time.time()
                # mols = find_neighbors(
                #     Molecule(smiles), depth, evaluations, apply_evaluation
                # )
                mols = []
                duration = time.time() - start
                valid_mols = []
                start = time.time()
                for mol in mols:
                    valid = True
                    for eval_ in evaluations:
                        try:
                            eval_.evaluate(mol)
                        except evaluator.EvaluationError:
                            valid = False
                            break
                    if valid:
                        valid_mols.append(mol.id_representation.smiles)
                duration_post_filter = time.time() - start
                print(
                    ",".join(
                        map(
                            str,
                            [
                                smiles,
                                depth,
                                len(mols),
                                len(set(mols)),
                                len(valid_mols),
                                len(set(valid_mols)),
                                apply_evaluation,
                                f"{duration:.2f}",
                                f"{duration_post_filter:.2f}",
                            ],
                        )
                    )
                )
                print(set(valid_mols))

    # actions = mol.list_all_possible_actions()
    # nb_possible_actions = mol.nb_possible_actions()
    # count = 0
    # for _, action_space in actions.items():
    #     for _, actions_list in action_space.items():
    #         for action in actions_list:
    #             # print(action, "->", end=" ")
    #             new_mol = action.apply()
    #             # print(new_mol, end="")
    #             valid = True
    #             for eval_ in evaluations:
    #                 try:
    #                     eval_.evaluate(new_mol)
    #                 except evaluator.EvaluationError:
    #                     valid = False
    #                     # print(" - invalid")
    #                     break
    #             if valid:
    #                 # print(" - invalid")
    #                 count += 1
    # duration = time.time() - start
    # print(f"nb valid actions {count}/{nb_possible_actions} - {duration:.2f}s")
