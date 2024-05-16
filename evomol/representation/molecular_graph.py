"""
Molecular graph representation of a molecule with RDKit.
"""

from __future__ import annotations

import networkx as nx
import numpy as np
import numpy.typing as npt
from rdkit import Chem
from rdkit.Chem import rdchem
from typing_extensions import override

from evomol.representation.molecule import MoleculeRepresentation, max_valence


class MolecularGraph(MoleculeRepresentation):
    """
    Molecule representation using a molecular graph with RDKit.
    """

    def __init__(
        self,
        smiles: str,
        sanitize: bool = False,
        mutability: bool = True,
    ):
        super().__init__(smiles)
        self.mol: rdchem.RWMol = rdchem.RWMol(Chem.MolFromSmiles(smiles))
        self.sanitize: bool = sanitize
        self.mutability = mutability

        for atom in self.mol.GetAtoms():
            atom.SetBoolProp("mutability", mutability)

        self.update_representation()

        # @classmethod
        # def set_generic_action_space(cls, accepted_atoms: list[str]) -> None:
        # """Set the generic action space for the molecular graph representation."""
        # from evomol.representation.molecular_graph_actions import (
        #     ActionSpaceAddAtomMolGraph,
        #     ActionSpaceChangeBondMolGraph,
        #     ActionSpaceCutAtomMolGraph,
        #     ActionSpaceInsertCarbonMolGraph,
        #     ActionSpaceMoveFunctionalGroupMolGraph,
        #     ActionSpaceRemoveAtomMolGraph,
        #     ActionSpaceRemoveGroupMolGraph,
        #     ActionSpaceSubstituteAtomMolGraph,
        # )

        # accepted_substitutions: dict[str, list[str]] = {
        #     accepted_atom: [
        #         other_accepted_atom
        #         for other_accepted_atom in accepted_atoms
        #         if other_accepted_atom not in {accepted_atom, "H"}
        #     ]
        #     for accepted_atom in accepted_atoms
        #     if accepted_atom != "H"
        # }

        # cls.set_action_spaces(
        #     [
        #         ActionSpaceAddAtomMolGraph(
        #             accepted_atoms=accepted_atoms, allow_bonding=False
        #         ),
        #         ActionSpaceChangeBondMolGraph(avoid_break_bond=False),
        #         ActionSpaceCutAtomMolGraph(),
        #         ActionSpaceInsertCarbonMolGraph(),
        #         ActionSpaceMoveFunctionalGroupMolGraph(),
        #         ActionSpaceRemoveAtomMolGraph(),
        #         ActionSpaceRemoveGroupMolGraph(only_remove_smallest_group=True),
        #         ActionSpaceSubstituteAtomMolGraph(
        #             accepted_substitutions=accepted_substitutions
        #         ),
        #     ]
        # )

    @override
    def representation(self) -> str:
        return self.canonical_smiles

    @property
    def nb_atoms(self) -> int:
        """Return the number of atoms in the molecule."""
        nb_atoms: int = self.mol.GetNumAtoms()
        return nb_atoms

    @property
    def smiles(self) -> str:
        """Return the SMILES representation of the molecule."""
        smiles: str = Chem.MolToSmiles(self.mol)
        return smiles

    @property
    def canonical_smiles(self) -> str:
        """Return the canonical SMILES representation of the molecule."""
        canonic_smiles: str = Chem.MolToSmiles(
            Chem.MolFromSmiles(
                Chem.MolToSmiles(MolecularGraph(Chem.MolToSmiles(self.mol)).mol)
            )
        )
        return canonic_smiles

    @property
    def adjacency_matrix(self) -> list[list[int]]:
        """Return the adjacency matrix of the molecule."""
        adjacency_matrix: list[list[int]] = Chem.GetAdjacencyMatrix(self.mol)
        return adjacency_matrix

    def __copy__(self) -> MolecularGraph:
        mol = MolecularGraph(self.canonical_smiles, self.sanitize, self.mutability)
        return mol

    def atom_mutability(self, atom_idx: int) -> bool:
        """Return the mutability of the atom at the given index."""
        atom_mutability: bool = self.mol.GetAtomWithIdx(atom_idx).GetBoolProp(
            "mutability"
        )
        return atom_mutability

    def set_atom_mutability(self, atom_idx: int, mutability: bool) -> None:
        """Set the mutability of the atom at the given index."""
        self.mol.GetAtomWithIdx(atom_idx).SetBoolProp("mutability", mutability)

    def __repr__(self) -> str:
        return f"MolecularGraph({self.canonical_smiles})"

    def update_representation(self) -> None:
        """Update internal RDKit representation of the molecular graph.
        This method should be called after each action on the molecular graph.
        """
        # sanitize molecule if needed
        if self.sanitize:
            Chem.SanitizeMol(self.mol)
            # clear computed properties
            self.mol.ClearComputedProps()

        # kekulize the molecular graph
        # https://chemistry.stackexchange.com/questions/116498/what-is-kekulization-in-rdkit
        Chem.Kekulize(self.mol)

        # setting all atoms to non aromatics
        for atom in self.mol.GetAtoms():
            atom.SetIsAromatic(False)

        # setting all bonds to non aromatics
        for bond in self.mol.GetBonds():
            bond.SetIsAromatic(False)

        # update the property cache of atoms
        for i in range(self.mol.GetNumAtoms()):
            for j in range(self.mol.GetNumAtoms()):
                bond = self.mol.GetBondBetweenAtoms(i, j)
                if bond is not None:
                    bond.SetIsAromatic(False)

        # update RDKit representation
        self.mol.UpdatePropertyCache()
        Chem.FastFindRings(self.mol)

    def add_atom(self, atom_type: str) -> None:
        """Add an atom to the molecular graph.

        Args:
            atom (str):
                atom to add to the molecular graph.
        """
        atom = Chem.Atom(atom_type)
        atom.SetBoolProp("mutability", True)
        self.mol.AddAtom(atom)

    def remove_atom(self, atom_idx: int) -> None:
        """Removing the atom at the given index from the molecular graph.

        Args:
            atom_idx (int):
                atom index
        """
        self.mol.RemoveAtom(atom_idx)

    def replace_atom(self, atom_idx: int, atom_type: str) -> None:
        """Replace the atom at the given index with the given atom type.

        Args:
            atom_idx (int):
                atom index
            atom_type (str):
                atom to replace the atom at the given index.
        """
        # Changing atomic number
        self.mol.GetAtomWithIdx(atom_idx).SetAtomicNum(
            Chem.GetPeriodicTable().GetAtomicNumber(atom_type)
        )

        # Setting formal charge to 0
        self.mol.GetAtomWithIdx(atom_idx).SetFormalCharge(0)

    def set_bond(self, atom1_idx: int, atom2_idx: int, bond_type_num: int) -> None:
        """Set the bond between the two atoms to the given bond type.

        Args:
            atom1_idx (int):
                index of the first atom
            atom2_idx (int):
                index of the second atom
            bond_type_num (int):
                bond type
        """
        # Extracting current bond between given atoms
        bond = self.mol.GetBondBetweenAtoms(int(atom1_idx), int(atom2_idx))

        bond_type = None
        if bond_type_num == 0:
            bond_type = None  # no bond or Chem.BondType.ZERO ? TODO
        elif bond_type_num == 1:
            bond_type = Chem.BondType.SINGLE
        elif bond_type_num == 2:
            bond_type = Chem.BondType.DOUBLE
        elif bond_type_num == 3:
            bond_type = Chem.BondType.TRIPLE
        else:
            raise ValueError(f"Unknown bond type {bond_type_num} (use 0, 1, 2 or 3)")

        if bond is None:
            # Adding a new bond
            self.mol.AddBond(atom1_idx, atom2_idx, bond_type)
        elif bond_type is None:
            # Removing the bond
            self.mol.RemoveBond(atom1_idx, atom2_idx)
        else:
            # Changing the bond type
            bond.SetBondType(bond_type)

    def bridge_bonds_matrix(self) -> npt.NDArray[np.bool_]:
        """Return a boolean matrix of size (nb_defined_atoms, nb_defined_atoms)
        representing whether bonds of the molecule are bridge bonds.
        """
        # convert the adjacency matrix to a networkx graph
        adjacency_matrix = nx.from_numpy_array(self.adjacency_matrix)

        # initialize the output matrix
        bridge_bonds_matrix = np.full((self.nb_atoms, self.nb_atoms), False, dtype=bool)

        # find the bridge bonds and set the corresponding matrix elements to True
        for atom1_idx, atom2_idx in nx.bridges(adjacency_matrix):
            bridge_bonds_matrix[atom1_idx, atom2_idx] = True
            bridge_bonds_matrix[atom2_idx, atom1_idx] = True

        return bridge_bonds_matrix

    def atom_degree(self, atom_idx: int, as_multigraph: bool) -> int:
        """Return the degree of the atom at the given index.

        Args:
            atom_idx (int):
                atom index
            as_multigraph (bool, optional):
                if True, the degree is the number of bonds (including multiple
                bonds) connected to the atom.
                If False, the degree is the number of atoms connected to the atom.

        Returns:
            int: degree of the atom at the given index
        """
        if as_multigraph:
            return self._explicit_valence(atom_idx)
        degree: int = self.mol.GetAtomWithIdx(atom_idx).GetDegree()
        return degree

    def formal_charge_vector(self) -> list[int]:
        """Return the formal charge of each atom in the molecule.

        Returns:
            list[int]: formal charge of each atom in the molecule
        """
        return [
            self.mol.GetAtomWithIdx(atom_idx).GetFormalCharge()
            for atom_idx in range(self.nb_atoms)
        ]

    def free_electrons_vector(self) -> npt.NDArray[np.int32]:
        """Return the number of free electrons of each atom in the molecule.

        Returns:
            np.ndarray[int]: number of free electrons of each atom in the molecule
        """
        return np.array(
            [self._nb_free_electrons(atom_idx) for atom_idx in range(self.nb_atoms)]
        )

    def atom_type(self, atom_idx: int) -> str:
        """Return the atom type of the atom at the given index.

        Args:
            atom_idx (int):
                atom index

        Returns:
            str: atom type of the atom at the given index
        """
        symbol: str = self.mol.GetAtomWithIdx(atom_idx).GetSymbol()
        return symbol

    def bond_type_num(self, atom1_idx: int, atom2_idx: int) -> int:
        """Return the bond type between the two atoms.

        Args:
            atom1_idx (int):
                index of the first atom
            atom2_idx (int):
                index of the second atom

        Returns:
            int: bond type between the two atoms
        """
        bond = self.mol.GetBondBetweenAtoms(atom1_idx, atom2_idx)
        if bond is None:
            return 0
        return int(bond.GetBondType())

    def _explicit_valence(self, atom_idx: int) -> int:
        """Return the explicit valence of the atom at the given index.

        Args:
            atom_idx (int):
                atom index

        Returns:
            int: explicit valence of the atom at the given index
        """
        atom = self.mol.GetAtomWithIdx(atom_idx)
        atom.UpdatePropertyCache()
        valence: int = atom.GetExplicitValence()
        return valence

    def _nb_free_electrons(self, atom_idx: int) -> int:
        """Return the number of free electrons of the atom at the given index.

        Args:
            atom_idx (int):
                atom index

        Returns:
            int: number of free electrons of the atom at the given index
        """
        return max_valence(
            self.mol.GetAtomWithIdx(atom_idx).GetSymbol()
        ) - self._explicit_valence(atom_idx)

    def _is_new_bond_possible(self, atom1_idx: int, atom2_idx: int) -> bool:
        """Return whether a new bond between the two atoms is possible.

        Args:
            atom1_idx (int):
                index of the first atom
            atom2_idx (int):
                index of the second atom

        Returns:
            bool: whether a new bond between the two atoms is possible
        """
        return (
            min(self._nb_free_electrons(atom1_idx), self._nb_free_electrons(atom2_idx))
            > 0
        )
