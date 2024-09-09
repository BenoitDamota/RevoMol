"""
Script to try to get the Morgan fingerprints in a different way than the old one
that is deprecated.

Should be careful about the version of rdkit as the fingerprint could be different
between versions.
https://github.com/rdkit/rdkit/issues/2018#issuecomment-800014238
but I didn't find differences in the results between the old and the new way to
get the fingerprints between rdkit 2022 and 2024.

You can look at :
https://greglandrum.github.io/rdkit-blog/posts/2023-01-18-fingerprint-generator-tutorial.html
for the new way to get the fingerprints but it seems to always use modulo.
"""

from rdkit import Chem
from rdkit.Chem import AllChem


def main() -> None:
    """Try to get the Morgan fingerprints in a different way."""
    smiles = "CNOC#CS"

    radius = 2

    mol = Chem.MolFromSmiles(smiles)

    # get the Morgan fingerprints - old way
    fingerprints = AllChem.GetMorganFingerprint(mol, radius=radius).GetNonzeroElements()

    # get the Morgan fingerprints - new way
    fingerprints_generator = AllChem.GetMorganGenerator(radius=2)

    ecfp = fingerprints_generator.GetSparseCountFingerprint(mol).GetNonzeroElements()

    print("old :", fingerprints)
    print("new :", ecfp)


def get_fingerprints() -> None:
    """script to get the fingerprints with the new way."""
    fingerprints_generator = AllChem.GetMorganGenerator(radius=2)

    file = "external_data/smiles_datasets/ZINC.csv"
    output = "output_file_ecfp_new_zinc.txt"

    file = "external_data/smiles_datasets/ChEMBL.csv"
    output = "output_file_ecfp_new_chembl.txt"

    with open(file, encoding="utf-8") as f:
        f.readline()
        with open(output, "w", encoding="utf-8") as out:
            for line in f.readlines():
                smiles = line.split(",")[0]
                mol = Chem.MolFromSmiles(smiles)

                # get the Morgan fingerprints - new way
                ecfp = fingerprints_generator.GetSparseCountFingerprint(
                    mol
                ).GetNonzeroElements()

                out.write(f"{smiles}")
                for e in ecfp:
                    out.write(f",{e}")
                out.write("\n")


def compare_fingerprints() -> None:
    """after getting the fingerprints, compare the old and the new way.
    need to change rdkit version between the two."""

    input_file = "output_file_ecfp_new_zinc.txt"
    input_file = "output_file_ecfp_new_chembl.txt"

    with open(input_file, encoding="utf-8") as f:
        for line in f.readlines():
            smiles, *ecfp = line.strip().split(",")
            ecfp_s = set(ecfp)
            mol = Chem.MolFromSmiles(smiles)

            # get the Morgan fingerprints - old way
            fingerprints = AllChem.GetMorganFingerprint(
                mol, radius=2
            ).GetNonzeroElements()

            old_ecfp = set()
            for e in fingerprints:
                old_ecfp.add(str(e))
            if old_ecfp != ecfp_s:
                print("error ", smiles, ecfp, old_ecfp)


if __name__ == "__main__":
    main()
