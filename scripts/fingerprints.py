"""
Script to try to get the Morgan fingerprints in a different way than the old one
that is deprecated.


It doesn't work, the fingerprints are modulo 2048 so it's not what we want.

Should be careful about the version of rdkit as the fingerprint could be different
between versions.
https://github.com/rdkit/rdkit/issues/2018#issuecomment-800014238

You can look at :
https://greglandrum.github.io/rdkit-blog/posts/2023-01-18-fingerprint-generator-tutorial.html
for the new way to get the fingerprints but it seems to always use modulo.
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, rdFingerprintGenerator


def main() -> None:
    """Try to get the Morgan fingerprints in a different way."""
    smiles = "CNOC#CS"

    radius = 2

    mol = Chem.MolFromSmiles(smiles)

    # get the Morgan fingerprints - old way
    fingerprints = AllChem.GetMorganFingerprint(mol, radius=radius).GetNonzeroElements()

    # get the Morgan fingerprints - new way
    morgan_fingerprint_generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=radius, fpSize=2**10
    )

    fp = morgan_fingerprint_generator.GetFingerprint(mol)

    fp_ = np.zeros((0,), dtype=int)

    DataStructs.ConvertToNumpyArray(fp, fp_)

    print("old :", fingerprints)
    print("new :", fp_.nonzero())


if __name__ == "__main__":
    main()
