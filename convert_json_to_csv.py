from multiprocessing import Process, Manager, Queue
import json

import typer

from evomol import evaluation as evaluator
from evomol.representation import SMILES, MolecularGraph, Molecule


def evaluateMolecule(
    inputQueue,
    results,
    eval_gcf_chembl,
    eval_ecfp_chembl,
    eval_gcf_chembl_zinc,
    eval_ecfp_chembl_zinc,
):
    while True:
        item = inputQueue.get()
        if item == "STOP":
            break
        smi, swscore_ChEMBL, swscore_ZINC = item
        mol = Molecule(smi)
        mol_graph = mol.get_representation(MolecularGraph)
        nb_atoms = mol_graph.nb_atoms
        if nb_atoms > 12:
            continue
        for atom in mol_graph.atoms:
            if atom not in Molecule.accepted_atoms:
                continue
        can_smi = mol_graph.canonical_smiles
        nb_unknown_gcf_chembl = eval_gcf_chembl.evaluate(mol)
        nb_unknown_ecfp_chembl = eval_ecfp_chembl.evaluate(mol)
        nb_unknown_gcf_chembl_zinc = eval_gcf_chembl_zinc.evaluate(mol)
        nb_unknown_ecfp_chembl_zinc = eval_ecfp_chembl_zinc.evaluate(mol)
        results.append(
            ",".join(
                map(
                    str,
                    (
                        smi,
                        can_smi,
                        nb_atoms,
                        swscore_ChEMBL,
                        swscore_ZINC,
                        nb_unknown_gcf_chembl,
                        nb_unknown_ecfp_chembl,
                        nb_unknown_gcf_chembl_zinc,
                        nb_unknown_ecfp_chembl_zinc,
                    ),
                )
            )
        )


def convert_dataset(input_json: str, output_csv: str, num_proc: int):
    Molecule.id_representation_class = SMILES
    Molecule.representations_class = [MolecularGraph]
    Molecule.accepted_atoms = ["C", "O", "N", "F", "S"]

    manager = Manager()
    results = manager.list()
    inputQueue = Queue()

    eval_gcf_chembl = evaluator.UnknownGCF(path_db="external_data/gcf1.txt")
    eval_ecfp_chembl = evaluator.UnknownECFP(
        path_db="external_data/ecfp4_ChEMBL.txt", radius=2
    )
    eval_gcf_chembl_zinc = evaluator.UnknownGCF(path_db="external_data/gcf2.txt")
    eval_ecfp_chembl_zinc = evaluator.UnknownECFP(
        path_db="external_data/ecfp4_ChEMBL_ZINC.txt", radius=2
    )

    processes = [
        Process(
            target=evaluateMolecule,
            args=(
                inputQueue,
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

    with open(input_json) as f:
        data = json.load(f)
    for smi in data:
        inputQueue.put((smi, data[smi]["swscore_ChEMBL"], data[smi]["swscore_ZINC"]))

    for _ in range(num_proc):
        inputQueue.put("STOP")

    for p in processes:
        p.join()

    # write results in csv file
    with open(output_csv, "w") as f:
        f.write(
            "smiles,smiles_canonical,nb_atoms,sw_filter_chembl,sw_filter_zinc,"
            "nb_unknown_ecfp_chembl,nb_unknown_gcf_chembl,"
            "nb_unknown_ecfp_chembl_zinc,nb_unknown_gcf_chembl_zinc\n"
        )
        for result in results:
            f.write(result + "\n")


if __name__ == "__main__":

    # python convert_json_to_csv.py output/other_datasets/Evo10_neutral_dict.json output/other_datasets/Evo10.csv 32
    # python convert_json_to_csv.py output/other_datasets/ChEMBL_neutral_dict.json output/other_datasets/ChEMBL.csv 32
    # python convert_json_to_csv.py output/other_datasets/GDBChEMBL_neutral_dict.json output/other_datasets/GDBChEMBL.csv 32
    # python convert_json_to_csv.py output/other_datasets/ZINC_neutral_dict.json output/other_datasets/ZINC.csv 32

    typer.run(convert_dataset)

    # database = "ChEMBL"
    # database = "GDBChEMBL"
    # database = "ZINC"
    # database = "Evo10"

    # input_json = f"output/other_datasets/{database}_neutral_dict.json"

    # output_csv = f"output/other_datasets/{database}.csv"

    # convert_dataset(input_json, output_csv)
