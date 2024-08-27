"""
Début de script pour tester la diversité.
Le code est à finir.
"""

import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

# pylint: disable=wrong-import-position, import-error

from evomol import default_parameters as dp
from evomol import evaluation as evaluator
from evomol.representation import Molecule


def main() -> None:
    """Test the diversity of a molecule."""

    dp.setup_default_parameters(
        accepted_atoms=["C", "O", "N", "F", "S"],
        max_heavy_atoms=38,
    )
    dp.setup_default_action_space(
        with_add_group=False,
        with_remove_group=False,
    )

    # Look at article Scalable estimator of the diversity of the size
    descriptors = [
        # evaluator.Descriptor(
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


if __name__ == "__main__":
    main()
