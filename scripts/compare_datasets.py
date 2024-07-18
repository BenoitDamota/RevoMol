"""
Compare data enumerated by EvoMol with other datasets

Enumerated datasets uses canonical and non aromatic smiles and are filtered
with ECFP and GCF filters (chembl or chembl_zinc depending on the name of the 
file)
The max number of atom allowed to generate the molecules is indicated in the 
name of the file.

for the other datasets, they are under csv format with the following columns:
- smiles_aromatic: SMILES of the molecule from the original json dataset
- smiles_kekulized: SMILES of the molecule after kekulization 
    (non aromatic and canonical)
- nb_atoms: number of atoms in the molecule
- sw_filter_chembl: result of chembl filter from the original json dataset
- sw_filter_zinc: result of zinc filter from the original json dataset
- nb_unknown_gcf_chembl: number of unknown GCF from the chembl database
- nb_unknown_ecfp_chembl: number of unknown ECFP4 from the chembl database
- nb_unknown_gcf_chembl_zinc: number of unknown GCF from the chembl_zinc database
- nb_unknown_ecfp_chembl_zinc: number of unknown ECFP4 from the chembl_zinc database
"""

import os
import sys

import pandas as pd

# Add the parent directory to the path to import the module evomol
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from evomol.representation import SMILES, MolecularGraph, Molecule


def compare_our_data(nb_heavy_atoms: int):
    eval_name = "chembl"  # "chembl_zinc"
    # own_data contains for each dataset, for each nb_heavy_atoms,
    # the list of valid molecules
    own_data = {}
    own_data["chembl"] = {
        nb_heavy_atom: [list() for _ in range(0, nb_heavy_atoms + 1)]
        for nb_heavy_atom in range(1, nb_heavy_atoms + 1)
    }
    own_data["chembl_zinc"] = {
        nb_heavy_atom: [list() for _ in range(0, nb_heavy_atoms + 1)]
        for nb_heavy_atom in range(1, nb_heavy_atoms + 1)
    }

    for eval_name in ["chembl", "chembl_zinc"]:
        for nb_heavy_atom in range(1, nb_heavy_atoms + 1):
            file = (
                f"output/enumeration_from_C_{eval_name}"
                f"_max_atom_{nb_heavy_atom}_depth_2.txt"
            )
            with open(file) as f:
                for line in f.readlines():
                    s = line.strip().split()[0]
                    if s == '""':
                        s = ""
                    nb_atoms = Molecule(s).get_representation(MolecularGraph).nb_atoms
                    own_data[eval_name][nb_heavy_atom][nb_atoms].append(s)

    print("nb_atoms,nb_max_atoms,nb_molecules_chembl,nb_molecules_chembl_zinc")
    for nb_atoms in range(1, nb_heavy_atoms + 1):
        for nb_heavy_atom in range(nb_atoms, nb_heavy_atoms + 1):
            # for eval_name in ["chembl", "chembl_zinc"]:
            print(
                f"{nb_atoms},{nb_heavy_atom},"
                f"{len(own_data['chembl'][nb_heavy_atom][nb_atoms])},"
                f"{len(own_data['chembl_zinc'][nb_heavy_atom][nb_atoms])}"
            )

    for nb_atoms in range(1, nb_heavy_atoms + 1):
        for nb_heavy_atom in range(nb_atoms, nb_heavy_atoms + 1):
            # diff between nb_atoms and nb_atoms - 1
            if not own_data["chembl"][nb_heavy_atom - 1][nb_atoms]:
                continue
            print(
                f"chembl {nb_heavy_atom} {nb_atoms} - "
                f"{nb_heavy_atom- 1} {nb_atoms} :",
                set(own_data["chembl"][nb_heavy_atom][nb_atoms])
                - set(own_data["chembl"][nb_heavy_atom - 1][nb_atoms]),
            )
            print(
                f"chembl_zinc {nb_heavy_atom} {nb_atoms} - "
                f"{nb_heavy_atom- 1} {nb_atoms} :",
                set(own_data["chembl_zinc"][nb_heavy_atom][nb_atoms])
                - set(own_data["chembl_zinc"][nb_heavy_atom - 1][nb_atoms]),
            )


def differences_filters(dataset: str):
    file = f"external_data/smiles_datasets/{dataset}.csv"
    data = pd.read_csv(file)

    data_ = data[
        (data["sw_filter_chembl"] == 1.0) & (data["nb_unknown_ecfp_chembl"] != 0)
    ]
    if data_.shape[0] > 0:
        data_.to_csv(
            f"output/diff_filter/new_filtered_chembl_{dataset}.csv", index=False
        )

    data_ = data[
        (data["sw_filter_zinc"] == 1.0) & (data["nb_unknown_ecfp_chembl_zinc"] != 0)
    ]
    if data_.shape[0] > 0:
        data_.to_csv(
            f"output/diff_filter/new_filtered_chembl_zinc_{dataset}.csv", index=False
        )

    data_ = data[
        (data["sw_filter_chembl"] != 1.0) & (data["nb_unknown_ecfp_chembl"] == 0)
    ]
    if data_.shape[0] > 0:
        data_.to_csv(
            f"output/diff_filter/lost_filtered_chembl_{dataset}.csv", index=False
        )

    data_ = data[
        (data["sw_filter_zinc"] != 1.0) & (data["nb_unknown_ecfp_chembl_zinc"] == 0)
    ]
    if data_.shape[0] > 0:
        data_.to_csv(
            f"output/diff_filter/lost_filtered_chembl_zinc_{dataset}.csv", index=False
        )


def differences_smiles(dataset_name: str, max_heavy_atoms: int, filter_name: str):

    print(
        "filter,nb_max_atoms,dataset_name,size_dataset,size_enumeration"
        ",size_intersection,size_diff_dataset,size_diff_enumeration"
    )

    file = f"external_data/smiles_datasets/{dataset_name}.csv"
    dataset = pd.read_csv(file)
    dataset = dataset[dataset["nb_atoms"] <= max_heavy_atoms]
    if filter_name == "chembl":
        dataset = dataset[dataset["nb_unknown_ecfp_chembl"] == 0]
        dataset = dataset[dataset["nb_unknown_gcf_chembl"] == 0]
    if filter_name == "chembl_zinc":
        dataset = dataset[dataset["nb_unknown_ecfp_chembl_zinc"] == 0]
        dataset = dataset[dataset["nb_unknown_gcf_chembl_zinc"] == 0]

    file_enum = f"output/enumeration_from_C_{filter_name}_{max_heavy_atoms}.txt"
    data_enum = pd.read_csv(file_enum, header=None)

    print(
        ",".join(
            map(
                str,
                [
                    filter_name,
                    max_heavy_atoms,
                    dataset_name,
                    dataset.shape[0],
                    data_enum.shape[0],
                    len(set(dataset["smiles_kekulized"]) & set(data_enum[0])),
                    len(set(dataset["smiles_kekulized"]) - set(data_enum[0])),
                    len(set(data_enum[0]) - set(dataset["smiles_kekulized"])),
                ],
            )
        )
    )

    print(
        f"missing from enum compared to {dataset_name} : "
        f"{set(dataset['smiles_kekulized']) - set(data_enum[0])}"
    )


if __name__ == "__main__":
    Molecule.id_representation_class = SMILES
    Molecule.representations_class = [MolecularGraph]
    Molecule.accepted_atoms = ["C", "O", "N", "F", "S"]

    database = "ChEMBL"
    database = "GDBChEMBL"
    database = "ZINC"
    database = "Evo10"
    datasets = [
        "Evo10",
        "ChEMBL",
        "GDBChEMBL",
        "ZINC",
    ]

    max_nb_atoms = 8

    # to compare the impact of filters, it will create files
    # in output/diff_filter
    for database in datasets:
        differences_filters(database)

    # compare_our_data(max_nb_atoms)
    # differences_smiles(database, max_nb_atoms, "chembl")
    # differences_smiles(database, max_nb_atoms, "chembl_zinc")
