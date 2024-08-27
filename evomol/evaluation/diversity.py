"""
Diversity evaluation
"""

import math
import os
import tempfile
from typing import Callable

import rdkit
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Contrib.IFG import ifg
from typing_extensions import override

from evomol.evaluation.evaluation import Evaluation
from evomol.representation import MolecularGraph, Molecule


def scaffolds(molecule: Molecule) -> list[str]:
    """List of scaffolds from the given molecule.
    Warning: Does not work. Maybe the problem is the rdkit version.

    Args:
        molecule (Molecule): Molecule to extract the scaffolds.

    Returns:
        list[str]: List of scaffolds.
    """
    canonical_smiles = molecule.get_representation(
        MolecularGraph
    ).aromatic_canonical_smiles
    canonical_smiles = molecule.get_representation(MolecularGraph).canonical_smiles
    # chem_canonical_smiles = Chem.CanonSmiles(canonical_smiles)

    mol = Chem.MolFromSmiles(canonical_smiles)

    scaffold_smiles = MurckoScaffold.MurckoScaffoldSmiles(mol)

    smiles = Chem.MolToSmiles(scaffold_smiles)

    print(smiles)

    return [smiles]


def gen_scaffolds(molecule: Molecule) -> list[str]:
    """Extracting generic scaffolds from the given molecule.

    Args:
        molecule (Molecule): Molecule to extract the generic scaffolds.

    Returns:
        list[str]: List of generic scaffolds.
    """
    canonical_smiles = molecule.get_representation(
        MolecularGraph
    ).aromatic_canonical_smiles

    mol = Chem.MolFromSmiles(canonical_smiles)

    mol_scaffolds = MurckoScaffold.MakeScaffoldGeneric(mol)

    smiles = Chem.MolToSmiles(mol_scaffolds)

    return [smiles]


def compute_ifg(molecule: Molecule) -> list[str]:
    """Extracting ifg from the given molecule.

    Args:
        molecule (Molecule): Molecule to extract the ifg.

    Returns:
        list[str]: List of ifg.
    """
    curr_ifgs: list[rdkit.Contrib.IFG.ifg.IFG] = ifg.identify_functional_groups(
        Chem.MolFromSmiles(
            molecule.get_representation(MolecularGraph).aromatic_canonical_smiles
        )
    )
    results: list[str] = list(set(curr_ifg[2] for curr_ifg in curr_ifgs))
    return results


def atoms_list(molecule: Molecule) -> list[str]:
    """Extracting atoms from the given molecule.

    Args:
        molecule (Molecule): Molecule to extract the atoms.

    Returns:
        list[str]: List of atoms.
    """
    return list(set(molecule.get_representation(MolecularGraph).atoms))


def shg_1(molecule: Molecule, level: int = 1) -> list[str]:
    """
    Extracting up to the given level from the given smiles
    see https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0321-8
    """

    qry_shingles = set()

    radius_constr = level + 1

    # Reloading molecule to make it aromatic
    mol = Chem.MolFromSmiles(
        molecule.get_representation(MolecularGraph).aromatic_canonical_smiles
    )

    for atm_idx in range(mol.GetNumAtoms()):
        for r in range(1, radius_constr):
            bonds = Chem.AllChem.FindAtomEnvironmentOfRadiusN(mol, r, atm_idx)

            if not bonds:
                break

            # the reportedly faster method
            atoms = set()
            for bond_id in bonds:
                bond = mol.GetBondWithIdx(bond_id)
                atoms.add(bond.GetBeginAtomIdx())
                atoms.add(bond.GetEndAtomIdx())

            # Computed rooted shingle
            new_shingle = Chem.rdmolfiles.MolFragmentToSmiles(
                mol,
                list(atoms),
                bonds,
                0,
                0,
                False,
                False,
                atm_idx,
                True,
                False,
                False,
            )
            qry_shingles.add(new_shingle)

    return list(qry_shingles)


def checkmol(molecule: Molecule) -> list[str]:
    """
    Extracting checkmol descriptors from given molecular graph.
    see https://homepage.univie.ac.at/norbert.haider/cheminf/cmmm.html

    This code is not tested, look how to import checkmol in the code
    """

    obabel_cmd = "obabel" + ' "-:%s" -omol -O %s 2>/dev/null'

    checkmol_cmd = os.path.join("external_data", "checkmol") + " %s > %s"

    smiles = molecule.get_representation(MolecularGraph).aromatic_canonical_smiles

    with tempfile.NamedTemporaryFile(mode="w+", suffix=".mol", delete=True) as mol_fic:
        fic_name = mol_fic.name
        os.system(obabel_cmd % (smiles, fic_name))
        with tempfile.NamedTemporaryFile(
            mode="w+", suffix=".res", delete=True
        ) as mol_ficout:
            ficout_name = mol_ficout.name
            os.system(checkmol_cmd % (fic_name, ficout_name))
            lines = [line.strip() for line in mol_ficout.readlines()]
    if len(lines) == 0 or lines[0] == "unknown query file format!":
        return []
    return list(set(lines))


def ecfp4(molecule: Molecule) -> list[int]:
    """Extracting ecfp4 from given molecular graph.

    Args:
        molecule (Molecule): Molecule to extract the ecfp4.

    Returns:
        list[int]: List of ecfp4.
    """
    fp = Chem.AllChem.GetMorganFingerprint(
        Chem.MolFromSmiles(
            molecule.get_representation(MolecularGraph).aromatic_canonical_smiles
        ),
        2,
    )
    on_bits = fp.GetNonzeroElements().keys()
    return list(on_bits)


class Descriptor(Evaluation):
    """Descriptor evaluation"""

    def __init__(
        self,
        descriptor_name: str,
        function: Callable[[Molecule], list[int] | list[str]],
        nb_max_descriptors: int,
    ) -> None:
        super().__init__(f"Descriptor_{descriptor_name}")
        self.descriptor_name: str = descriptor_name
        self.nb_max_descriptors: int = nb_max_descriptors
        self.function: Callable[[Molecule], list[int] | list[str]] = function
        self.descriptor_to_index: dict[str | int, int] = {}

    @override
    def _evaluate(self, molecule: Molecule) -> list[int]:
        descriptors = self.function(molecule)

        indexes: list[int] = [self.get_index(descriptor) for descriptor in descriptors]
        return indexes

    def get_index(self, descriptor: str | int) -> int:
        """Get the index of the descriptor.

        Args:
            descriptor (str): Descriptor to get the index.

        Returns:
            int: Index of the descriptor.
        """
        if descriptor in self.descriptor_to_index:
            return self.descriptor_to_index[descriptor]
        index = len(self.descriptor_to_index)
        assert index < self.nb_max_descriptors
        self.descriptor_to_index[descriptor] = index
        return index


class Diversity:
    """Diversity evaluation"""

    def __init__(
        self,
        descriptor: Descriptor,
        size: int,
    ) -> None:
        self.size: int = size
        self.descriptor: Descriptor = descriptor
        self.count: list[int] = [0 for _ in range(descriptor.nb_max_descriptors)]
        self.entropy_descriptor: list[float] = [
            0 for _ in range(descriptor.nb_max_descriptors)
        ]
        self.entropy_total: float = 0

    def evaluate(
        self,
        new_molecule: Molecule,
        molecule_to_replace: Molecule,
    ) -> float:
        """Evaluate the diversity of the new molecule compared to the molecule
        to replace.

        Args:
            new_molecule (Molecule): New molecule to evaluate.
            molecule_to_replace (Molecule): Molecule to replace.

        Returns:
            float: Score of the diversity.
        """

        # evaluate the new molecule
        self.descriptor.evaluate(new_molecule)

        # descriptors lost if we lost the molecule to replace
        lost_descriptors: set[int] = set(
            molecule_to_replace.value(self.descriptor.name)
        )

        # descriptors gained if we add the new molecule
        new_descriptors: set[int] = set(new_molecule.value(self.descriptor.name))

        # unique descriptors lost and gained
        lost: set[int] = lost_descriptors.difference(new_descriptors)
        gain: set[int] = new_descriptors.difference(lost_descriptors)

        # entropy for lost and gained descriptors
        partial_entropy_lost: list[float] = [
            entropy_descriptor(self.count[descriptor] - 1)
            - self.entropy_descriptor[descriptor]
            for descriptor in lost
        ]

        partial_entropy_gained: list[float] = [
            entropy_descriptor(self.count[descriptor] + 1)
            - self.entropy_descriptor[descriptor]
            for descriptor in gain
        ]

        # score of entropy
        score: float = sum(partial_entropy_lost) + sum(partial_entropy_gained)

        return score

    def evaluate_molecule(self, molecule: Molecule) -> None:
        """Compute the descriptors of the molecule.

        Args:
            molecule (Molecule): Molecule to evaluate.
        """
        self.descriptor.evaluate(molecule)

    def add_molecule(self, molecule: Molecule) -> None:
        """Add the molecule to the diversity evaluation.

        Args:
            molecule (Molecule): Molecule to add.
        """
        # make it a call for observer pattern
        for descriptor_index in molecule.value(self.descriptor.name):
            self.count[descriptor_index] += 1

    def remove_molecule(self, molecule: Molecule) -> None:
        """Remove the molecule from the diversity evaluation.

        Args:
            molecule (Molecule): Molecule to remove.
        """
        # make it a call for observer pattern
        for descriptor_index in molecule.value(self.descriptor.name):
            self.count[descriptor_index] -= 1

    def update_entropy(self) -> None:
        """Update the entropy of the diversity evaluation."""
        # make it a call for observer pattern
        self.entropy_descriptor = [
            entropy_descriptor(counter) for counter in self.count
        ]
        self.entropy_total = sum(self.entropy_descriptor)


def entropy_descriptor(counter: int) -> float:
    """
    Entropy computation for one descriptor
    """
    if counter == 0:
        return 0.0
    # No need for the size as it is constant
    proportion: float = counter  # / size
    entropy_: float = -proportion * math.log(proportion)
    return entropy_


def entropy(count: list[int], size: int) -> float:
    """
    Entropy computation
    """
    proportion = [c / size for c in count if c > 0]
    entropy_ = -sum(p * math.log(p) for p in proportion)
    return entropy_
