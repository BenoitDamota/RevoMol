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

# pylint: disable=wrong-import-position, import-error

from evomol import default_parameters as dp
from evomol.representation import MolecularGraph, Molecule


def main() -> None:
    """
    Uncomment the functions to compare the data between the enumerated data and
    the other datasets.
    """

    dp.setup_default_parameters()

    datasets = [
        "Evo10",
        "ChEMBL",
        "GDBChEMBL",
        "ZINC",
    ]

    max_nb_atoms = 10

    # to compare the impact of filters, it will create files
    # in output/diff_filter
    for database in datasets:
        # differences_filters(database)
        # differences_smiles(database, max_nb_atoms, "chembl")
        differences_smiles(database, max_nb_atoms, "chembl_zinc")

    # compare_our_data(max_nb_atoms)


def compare_our_data(nb_heavy_atoms: int) -> None:
    """Compare the enumerated data between all level of max number of atoms.
    Print the results in a CSV format for the number of molecules of each size
    for each max number of atoms.
    Then print the differences of molecules not found for a number of atoms with
    less heavy atoms.
    For example, "N#CN=C1SCS1" with 7 heavy atoms is not found in the enumeration
    with 7 heavy atoms max but is found in the enumeration with 8 heavy atoms max.

    Args:
        nb_heavy_atoms (int): Maximum number of heavy atoms
    """

    # own_data contains for each filter, for each nb_heavy_atoms max,
    # for each molecule size, the list of valid molecules
    own_data: dict[str, dict[int, list[list[str]]]] = {}
    own_data["chembl"] = {
        nb_heavy_atom: [[] for _ in range(0, nb_heavy_atoms + 1)]
        for nb_heavy_atom in range(1, nb_heavy_atoms + 1)
    }
    own_data["chembl_zinc"] = {
        nb_heavy_atom: [[] for _ in range(0, nb_heavy_atoms + 1)]
        for nb_heavy_atom in range(1, nb_heavy_atoms + 1)
    }

    # read the enumerated data
    for eval_name in ["chembl", "chembl_zinc"]:
        for nb_heavy_atom in range(1, nb_heavy_atoms + 1):
            file = os.path.join(
                "output",
                "enumeration_C_depth_2",
                f"enumeration_C_{eval_name}_max_atom_{nb_heavy_atom}_depth_2.txt",
            )
            with open(file, encoding="utf-8") as f:
                for line in f.readlines():
                    s = line.strip().split()[0]
                    if s == '""':
                        s = ""
                    # add the molecule to the right list depending on the number
                    # atoms
                    nb_atoms = Molecule(s).get_representation(MolecularGraph).nb_atoms
                    own_data[eval_name][nb_heavy_atom][nb_atoms].append(s)

    # print the header of the CSV and the number of molecules for each size
    print("nb_atoms,nb_max_atoms,nb_molecules_chembl,nb_molecules_chembl_zinc")
    for nb_atoms in range(1, nb_heavy_atoms + 1):
        for nb_heavy_atom in range(nb_atoms, nb_heavy_atoms + 1):
            # for eval_name in ["chembl", "chembl_zinc"]:
            print(
                f"{nb_atoms},{nb_heavy_atom},"
                f"{len(own_data['chembl'][nb_heavy_atom][nb_atoms])},"
                f"{len(own_data['chembl_zinc'][nb_heavy_atom][nb_atoms])}"
            )

    # print the differences between the enumeration with a number of heavy atoms
    # and the enumeration with higher number of heavy atoms
    for nb_atoms in range(2, nb_heavy_atoms):
        for nb_heavy_atom in range(nb_atoms, nb_heavy_atoms + 1):
            for nb_heavy_atom_2 in range(nb_heavy_atom + 1, nb_heavy_atoms + 1):
                diff = set(own_data["chembl"][nb_heavy_atom_2][nb_atoms]) - set(
                    own_data["chembl"][nb_heavy_atom][nb_atoms]
                )
                print(
                    f"chembl {nb_atoms} {nb_heavy_atom_2} - "
                    f"{nb_atoms} {nb_heavy_atom} :",
                    diff,
                )
                diff = set(own_data["chembl_zinc"][nb_heavy_atom_2][nb_atoms]) - set(
                    own_data["chembl_zinc"][nb_heavy_atom][nb_atoms]
                )

                print(
                    f"chembl_zinc {nb_atoms} {nb_heavy_atom_2} - "
                    f"{nb_atoms} {nb_heavy_atom} :",
                    diff,
                )


def differences_filters(dataset: str) -> None:
    """Look for differences between the filters in the older version of the
    code and the new version of the code.
    Will create files in output/diff_filter with the new filtered molecules and
    the lost filtered molecules if there are any.

    Args:
        dataset (str): Name of the dataset (Evo10, ChEMBL, GDBChEMBL, ZINC)
    """
    file = os.path.join("external_data", "smiles_datasets", f"{dataset}.csv")
    data = pd.read_csv(file)

    # when the old filter is not applied but the new filter is applied (chembl)
    data_ = data[
        (data["sw_filter_chembl"] == 1.0) & (data["nb_unknown_ecfp_chembl"] != 0)
    ]
    if data_.shape[0] > 0:
        data_.to_csv(
            os.path.join(
                "output",
                "diff_filter",
                f"new_filtered_chembl_{dataset}.csv",
            ),
            index=False,
        )

    # when the old filter is applied but the new filter is not applied (chembl_zinc)
    data_ = data[
        (data["sw_filter_zinc"] == 1.0) & (data["nb_unknown_ecfp_chembl_zinc"] != 0)
    ]
    if data_.shape[0] > 0:
        data_.to_csv(
            os.path.join(
                "output",
                "diff_filter",
                f"new_filtered_chembl_zinc_{dataset}.csv",
            ),
            index=False,
        )

    # when the old filter is applied but the new filter is not applied (chembl)
    data_ = data[
        (data["sw_filter_chembl"] != 1.0) & (data["nb_unknown_ecfp_chembl"] == 0)
    ]
    if data_.shape[0] > 0:
        data_.to_csv(
            os.path.join(
                "output",
                "diff_filter",
                f"lost_filtered_chembl_{dataset}.csv",
            ),
            index=False,
        )

    # when the old filter is applied but the new filter is not applied (chembl_zinc)
    data_ = data[
        (data["sw_filter_zinc"] != 1.0) & (data["nb_unknown_ecfp_chembl_zinc"] == 0)
    ]
    if data_.shape[0] > 0:
        data_.to_csv(
            os.path.join(
                "output",
                "diff_filter",
                f"lost_filtered_chembl_zinc_{dataset}.csv",
            ),
            index=False,
        )


def differences_smiles(
    dataset_name: str, max_heavy_atoms: int, filter_name: str
) -> None:
    """Look for differences between the enumerated molecules and the dataset
    (Evo10, ChEMBL, GDBChEMBL, ZINC).
    Show the size of the dataset, the size of the enumeration, the size of the
    intersection, the size of the difference between the dataset and the
    enumeration and the size of the difference between the enumeration and the
    dataset.
    Also print the smiles that are missing from the enumeration compared to the
    dataset.

    Args:
        dataset_name (str): Name of the dataset
        max_heavy_atoms (int): Maximum number of heavy atoms
        filter_name (str): chembl or chembl_zinc as filter
    """

    # print the header of the CSV
    print(
        "filter,nb_max_atoms,dataset_name,size_dataset,size_enumeration"
        ",size_intersection,size_diff_dataset,size_diff_enumeration"
    )

    # load the dataset
    file = os.path.join("external_data", "smiles_datasets", f"{dataset_name}.csv")
    dataset = pd.read_csv(file)
    dataset = dataset[dataset["nb_atoms"] <= max_heavy_atoms]
    # filter the dataset to keep only the molecules that pass the filters
    if filter_name == "chembl":
        dataset = dataset[
            (dataset["nb_unknown_ecfp_chembl"] == 0)
            & (dataset["nb_unknown_gcf_chembl"] == 0)
        ]
    if filter_name == "chembl_zinc":
        dataset = dataset[
            (dataset["nb_unknown_ecfp_chembl_zinc"] == 0)
            & (dataset["nb_unknown_gcf_chembl_zinc"] == 0)
        ]

    file_enum = os.path.join(
        "output",
        "enumeration_C_depth_2",
        f"enumeration_C_{filter_name}_max_atom_{max_heavy_atoms}_depth_2.txt",
    )
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
        f"{filter_name} - missing from enum compared to {dataset_name} : "
        f"{set(dataset['smiles_kekulized']) - set(data_enum[0])}"
    )


if __name__ == "__main__":
    main()
