"""
This script generates ECFP files from a list of SMILES.
"""

from rdkit import Chem
from rdkit.Chem import AllChem


def generate_ecfp(list_of_smiles: list[str]) -> set[int]:
    """Generate ECFP from a list of SMILES.

    Args:
        list_of_smiles (str): List of SMILES.

    Returns:
        set[int]: set of ECFP.
    """

    fingerprints_generator = AllChem.GetMorganGenerator(radius=2)

    ecfp = set()
    for smiles in list_of_smiles:
        mol = Chem.MolFromSmiles(smiles)

        # get the Morgan fingerprints
        ecfp.update(
            fingerprints_generator.GetSparseCountFingerprint(mol)
            .GetNonzeroElements()
            .keys()
        )

    return ecfp


def save_ecfp_to_file(ecfp: set[int], output_file: str) -> None:
    """Save ECFP to a file.

    Args:
        ecfp (set[int]): ECFP.
        output_file (str): Output file.
    """
    with open(output_file, "w", encoding="utf-8") as out:
        for e in ecfp:
            out.write(f"{e}\n")


def main() -> None:
    """Generate ECFP from a list of SMILES and save it to a file."""

    # Implement your own way to get the list of SMILES here
    smiles = [
        "O=S(=O)(O)c1ccccc1c1ccsn1",
        "O=S(=O)(O)c1ccccc1c1ccsn1",
    ]

    # Choose the name of the output file
    file_name = "example_generate_ecfp_file.txt"

    ecfp = generate_ecfp(smiles)

    print(f"Number of ECFP: {len(ecfp)}")

    save_ecfp_to_file(ecfp, file_name)


if __name__ == "__main__":
    main()
