"""
Convert json data (Evo10, ChEMBL, GDBChEMBL, ZINC) to csv file.
The maximum number of atoms is 12 and the atoms are limited to C, O, N, F, S.
Change default_parameters to modify the number of atoms and the atoms allowed
in convert_dataset.

Data origin :
Evo10_neutral_dict.json
https://figshare.com/articles/dataset/EVO10_dataset_dictionnary_with_SMILES_and_sillywalk_score/20054378?backTo=/collections/Sillywalks_scores_of_organic_chemistry_datasets_PC9_QM9_OD9_GDB11_/6041117

676 875 SMILES generated by EvoMol during enumeration experiments of the
chemical spaces as defined by the ECFP4 of either ChEMBL or ZINC+ChEMBL.
Only neutral singlet molecule without any atomic charges (formal or real)
composed of H, C, N, O, F and S.

Key: SMILES
value: dict {
    'HAC', # number of heavy atoms up to 10
    'swscore_ChEMBL',  # % of ECFP4 of the molecule that belong to ChEMBL
    'swscore_ZINC',  # % of ECFP4 of the molecule that belong to ZINC or ChEMBL
}

ZINC_neutral_dict.json
https://figshare.com/articles/dataset/ZINC_sub_set_dictionnary_of_CNOFS_molecules_with_SMILES_and_sillywalk_scores/20059400?backTo=/collections/Sillywalks_scores_of_organic_chemistry_datasets_PC9_QM9_OD9_GDB11_/6041117
All ZINC 20 neutral singlet molecules without any atomic charges
(formal or real) composed of only H, C, N, O, F.

Key: SMILES
value: dict {
    'nb_atoms', # number of heavy atoms
    'swscore_ChEMBL',  # % of ECFP4 of the molecule that belong to ChEMBL
    'swscore_ZINC',  # % of ECFP4 of the molecule that belong to ZINC or ChEMBL
}

ChEMBL_neutral_dict.json
https://figshare.com/articles/dataset/ChEMBL25_dataset_dictionnary_of_CNOFS_molecules_with_SMILES_and_sillywalk_score/20054381?backTo=/collections/Sillywalks_scores_of_organic_chemistry_datasets_PC9_QM9_OD9_GDB11_/6041117
Subset of 1 191 453 neutral singlet molecules without any atomic charges
(formal or real) with only H, C, N, O F and S.

Key: SMILES
value: dict {
    'HAC', # number of heavy atoms
    'swscore_ChEMBL',  # % of ECFP4 of the molecule that belong to ChEMBL
    'swscore_ZINC',  # % of ECFP4 of the molecule that belong to ZINC or ChEMBL
}

GDBChEMBL_neutral_dict.json
https://figshare.com/articles/dataset/GDBChEMBL_sub_set_dictionnary_of_CNOFS_molecules_with_SMILES_and_sillywalk_score/20769853?backTo=/collections/Sillywalks_scores_of_organic_chemistry_datasets_PC9_QM9_OD9_GDB11_/6041117
The 3 786 315 neutral singlet molecules without any atomic charges
(formal or real) composed of only H, C, N, O, F and S that belong to GDBChEMBL.

Key: SMILES
value: dict {
    'HAC', # number of heavy atoms
    'swscore_ChEMBL',  # % of ECFP4 of the molecule that belong to ChEMBL
    'swscore_ZINC',  # % of ECFP4 of the molecule that belong to ZINC or ChEMBL
}
"""

from __future__ import annotations

import json
import os
import sys
from multiprocessing import Manager, Process, Queue
from multiprocessing.managers import ListProxy

import typer

# Add the parent directory to the path to import the module evomol
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

# pylint: disable=wrong-import-position, import-error

from evomol import default_parameters as dp
from evomol import evaluation as evaluator
from evomol.representation import MolecularGraph, Molecule


def evaluate_molecule(
    input_queue: Queue,
    results: ListProxy[str],
    eval_gcf_chembl: evaluator.UnknownGCF,
    eval_ecfp_chembl: evaluator.UnknownECFP,
    eval_gcf_chembl_zinc: evaluator.UnknownGCF,
    eval_ecfp_chembl_zinc: evaluator.UnknownECFP,
) -> None:
    """Evaluate SMILES as long as there are some in the inputQueue

    Args:
        inputQueue (Queue[Optional[tuple[str, float, float]]]): Contains the
            SMILES to evaluate
        results (ListProxy[str]): List to store the results
        eval_gcf_chembl (evaluator.UnknownGCF): evaluator GCF chembl
        eval_ecfp_chembl (evaluator.UnknownECFP): evaluator ECFP chembl
        eval_gcf_chembl_zinc (evaluator.UnknownGCF): evaluator GCF chembl_zinc
        eval_ecfp_chembl_zinc (evaluator.UnknownECFP): evaluator ECFP chembl_zinc
    """
    while True:
        item: tuple[str, float, float] | None = input_queue.get()
        if item is None:
            break
        # get the SMILES and previous results from the json file
        smi, swscore_chembl, swscore_zinc = item

        mol = Molecule(smi)
        mol_graph = mol.get_representation(MolecularGraph)
        nb_atoms = mol_graph.nb_atoms

        # limit the number of atoms
        if nb_atoms > Molecule.max_heavy_atoms:
            continue

        # limit the possible atoms to C, O, N, F, S
        valid = True
        for atom in mol_graph.atoms:
            if atom not in Molecule.accepted_atoms:
                valid = False
                break
        if not valid:
            continue

        results.append(
            ",".join(
                map(
                    str,
                    (
                        smi,
                        mol_graph.canonical_smiles,
                        nb_atoms,
                        swscore_chembl,
                        swscore_zinc,
                        eval_gcf_chembl.evaluate(mol),
                        eval_ecfp_chembl.evaluate(mol),
                        eval_gcf_chembl_zinc.evaluate(mol),
                        eval_ecfp_chembl_zinc.evaluate(mol),
                    ),
                )
            )
        )


def convert_dataset(input_json: str, output_csv: str, num_proc: int) -> None:
    """Convert json data to csv file

    Args:
        input_json (str): path to the json file
        output_csv (str): path to the csv file
        num_proc (int): number of processes to use
    """

    # parameters for the Molecule class and evaluators
    dp.setup_default_parameters(
        accepted_atoms=["C", "O", "N", "F", "S"],
        max_heavy_atoms=12,
    )

    eval_gcf_chembl = evaluator.UnknownGCF(
        path_db=os.path.join("external_data", "gcf1.txt"), name="chembl"
    )
    eval_ecfp_chembl = evaluator.UnknownECFP(
        path_db=os.path.join("external_data", "ecfp4_ChEMBL.txt"),
        radius=2,
        name="chembl",
    )
    eval_gcf_chembl_zinc = evaluator.UnknownGCF(
        path_db=os.path.join("external_data", "gcf2.txt"), name="chembl_zinc"
    )
    eval_ecfp_chembl_zinc = evaluator.UnknownECFP(
        path_db=os.path.join("external_data", "ecfp4_ChEMBL_ZINC.txt"),
        radius=2,
        name="chembl_zinc",
    )

    # prepare the processes
    manager = Manager()
    results: ListProxy[str] = manager.list()
    input_queue: Queue = Queue()

    processes = [
        Process(
            target=evaluate_molecule,
            args=(
                input_queue,
                results,
                eval_gcf_chembl,
                eval_ecfp_chembl,
                eval_gcf_chembl_zinc,
                eval_ecfp_chembl_zinc,
            ),
        )
        for _ in range(num_proc)
    ]

    for p in processes:
        p.start()

    # load completely the json file
    with open(input_json, encoding="utf-8") as f:
        data = json.load(f)

    # extract data from json file and put it in the inputQueue
    for smi in data:
        input_queue.put(
            (
                smi,
                float(data[smi]["swscore_ChEMBL"]),
                float(data[smi]["swscore_ZINC"]),
            )
        )

    # put None in the inputQueue to stop the processes
    for _ in range(num_proc):
        input_queue.put(None)

    # wait for the processes to finish
    for p in processes:
        p.join()

    # write results in csv file
    with open(output_csv, "w", encoding="utf-8") as f:
        f.write(
            "smiles_aromatic,smiles_kekulized,nb_atoms,sw_filter_chembl,sw_filter_zinc,"
            "nb_unknown_gcf_chembl,nb_unknown_ecfp_chembl,"
            "nb_unknown_gcf_chembl_zinc,nb_unknown_ecfp_chembl_zinc\n"
        )
        for result in results:
            f.write(result + "\n")


if __name__ == "__main__":
    # scripts commands are too long for 1 line, I divided them in 4 lines

    # python scripts/convert_json_to_csv.py
    # external_data/smiles_datasets/Evo10_neutral_dict.json
    # external_data/smiles_datasets/Evo10.csv
    # 32

    # python scripts/convert_json_to_csv.py
    # external_data/smiles_datasets/ChEMBL_neutral_dict.json
    # external_data/smiles_datasets/ChEMBL.csv
    # 32

    # python scripts/convert_json_to_csv.py
    # external_data/smiles_datasets/GDBChEMBL_neutral_dict.json
    # external_data/smiles_datasets/GDBChEMBL.csv
    # 32

    # python scripts/convert_json_to_csv.py
    # external_data/smiles_datasets/ZINC_neutral_dict.json
    # external_data/smiles_datasets/ZINC.csv
    # 32

    typer.run(convert_dataset)

    # to run without typer :

    # database = "ChEMBL"
    # database = "GDBChEMBL"
    # database = "ZINC"
    # database = "Evo10"

    # input_json = os.path.join(
    #     "external_data",
    #     "smiles_datasets",
    #     f"{database}_neutral_dict.json",
    # )

    # output_csv = os.path.join(
    #     "external_data",
    #     "smiles_datasets",
    #     f"{database}.csv",
    # )

    # convert_dataset(input_json, output_csv, 1)
