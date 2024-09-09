"""
This script generates GCF files from a list of SMILES.
"""

import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

# pylint: disable=wrong-import-position, import-error

from evomol.evaluation.generic_cyclic_features import list_gcf


def generate_gcf(list_of_smiles: list[str]) -> set[str]:
    """Generate GCF from a list of SMILES.

    Args:
        list_of_smiles (list[str]): List of SMILES.

    Returns:
        set[str]: set of GCF.
    """

    gcf = set()
    for smiles in list_of_smiles:
        gcf.update(list_gcf(smiles))

    return gcf


def save_gcf_to_file(gcf: set[str], output_file: str) -> None:
    """Save GCF to a file.

    Args:
        gcf (set[str]): GCF.
        output_file (str): Output file.
    """
    with open(output_file, "w", encoding="utf-8") as out:
        for e in gcf:
            out.write(f"{e}\n")


def main() -> None:
    """Generate GCF from a list of SMILES and save it to a file."""

    # Implement your own way to get the list of SMILES here
    smiles = [
        "O=S(=O)(O)c1ccccc1c1ccsn1",
        "O=S(=O)(O)c1ccc1c1ccsn1",
    ]

    # Choose the name of the output file
    file_name = "example_generate_gcf_file.txt"

    gcf = generate_gcf(smiles)

    print(f"Number of GCF: {len(gcf)}")

    save_gcf_to_file(gcf, file_name)


if __name__ == "__main__":
    main()
