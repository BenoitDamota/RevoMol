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

    smiles = "C"

    mol = Molecule(smiles)
    print(f"{mol=}")

    actions = mol.compute_possible_actions()
    pprint_action_space(actions)
    for _, action_space in actions.items():
        for _, actions_list in action_space.items():
            for action in actions_list:
                print(action, "->", end=" ")
                new_mol = action.apply()
                print(new_mol)


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
        evaluator.UnknownECFP(),
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
