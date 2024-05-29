"""
Molecular graph representation of a molecule with RDKit.
"""

from __future__ import annotations
import io
import os

from PIL import Image
import networkx as nx
from typing_extensions import override
from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem.rdDepictor import Compute2DCoords
from rdkit.Chem import Draw


from evomol.representation.molecule import MoleculeRepresentation


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
        try:
            self.mol: rdchem.RWMol = rdchem.RWMol(Chem.MolFromSmiles(smiles))
        except Exception as e:
            raise ValueError(f"Error with SMILES {smiles}: {e}")
        self.sanitize: bool = sanitize
        self._mutability = mutability

        for atom in self.mol.GetAtoms():
            atom.SetBoolProp("mutability", mutability)

        # self.update_representation() # TODO remettre

    def __copy__(self) -> MolecularGraph:
        mol = MolecularGraph("", self.sanitize, self._mutability)
        mol.mol = Chem.RWMol(self.mol, True)
        for i, atom in enumerate(self.mol.GetAtoms()):
            mol.mol.GetAtomWithIdx(i).SetBoolProp(
                "mutability", atom.GetBoolProp("mutability")
            )
        mol.update_representation()
        # mol.update_smile()
        return mol

    def __eq__(self, value: object) -> bool:
        return (
            isinstance(value, MolecularGraph)
            and self.canonical_smiles == value.canonical_smiles
        )

    @override
    def representation(self) -> str:
        return self.canonical_smiles

    def __repr__(self) -> str:
        return f"MolecularGraph({self.canonical_smiles})"

    @property
    def nb_atoms(self) -> int:
        """Return the number of atoms in the molecule."""
        nb_atoms: int = self.mol.GetNumAtoms()
        return nb_atoms

    @property
    def smiles(self) -> str:
        """Return the SMILES representation of the molecule."""
        smiles: str = Chem.MolToSmiles(self.mol, canonical=False)
        return smiles

    # def update_smile(self) -> None:
    #     """Update the SMILES representation of the molecule."""
    #     self.smiles = Chem.MolToSmiles(self.mol)

    @property
    def canonical_smiles(self) -> str:
        """Return the canonical SMILES representation of the molecule."""
        canonic_smiles: str = Chem.MolToSmiles(
            Chem.MolFromSmiles(
                Chem.MolToSmiles(MolecularGraph(Chem.MolToSmiles(self.mol)).mol)
            )
        )
        # mol = MolecularGraph(Chem.MolToSmiles(self.mol)).mol
        # print(
        #     "in canonical_smiles",
        #     Chem.MolToSmiles(mol),
        #     " mol to smiles :",
        #     Chem.MolToSmiles(self.mol),
        # )
        # for atom in mol.GetAtoms():
        #     print(
        #         atom.GetSymbol(),
        #         atom.GetFormalCharge(),
        #         atom.GetImplicitValence(),
        #         atom.GetNoImplicit(),
        #         atom.GetNumRadicalElectrons(),
        #         atom.GetBoolProp("mutability") and not atom.GetNoImplicit(),
        #     )
        return canonic_smiles

    @property
    def adjacency_matrix(self) -> list[list[int]]:
        """Return the adjacency matrix of the molecule."""
        adjacency_matrix: list[list[int]] = Chem.GetAdjacencyMatrix(self.mol)
        return adjacency_matrix

    @property
    def bonds(self) -> list[tuple[int, int]]:
        """Return the list of bonds in the molecule."""
        return [
            (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            for bond in self.mol.GetBonds()
        ]

    @property
    def bridge_bonds(self) -> list[tuple[int, int]]:
        """Return a list of the pairs of atoms in bridge bonds."""
        return list(nx.bridges(nx.Graph(self.adjacency_matrix)))

    # [
    #         (atom1, atom2)
    #         for atom1, atom2 in nx.bridges(nx.Graph(self.adjacency_matrix))
    #     ]

    @property
    def bridge_bonds_matrix(self) -> list[list[bool]]:
        """Return a boolean matrix of size (nb_atoms, nb_atoms)
        representing whether bonds of the molecule are bridge bonds.
        """
        matrix = [[False] * self.nb_atoms for _ in range(self.nb_atoms)]

        for atom1, atom2 in nx.bridges(nx.Graph(self.adjacency_matrix)):
            matrix[atom1][atom2] = True
            matrix[atom2][atom1] = True

        return matrix

    @property
    def atoms(self) -> list[str]:
        """Return the list of atoms in the molecule."""
        return [atom.GetSymbol() for atom in self.mol.GetAtoms()]

    @property
    def formal_charges(self) -> list[int]:
        """Return the formal charge of each atom in the molecule.

        Returns:
            list[int]: formal charge of each atom in the molecule
        """
        return [
            self.mol.GetAtomWithIdx(atom_idx).GetFormalCharge()
            for atom_idx in range(self.nb_atoms)
        ]

    @property
    def explicit_valences(self) -> list[int]:
        """Return the explicit valence of each atom in the molecule."""
        return [self._explicit_valence(atom_idx) for atom_idx in range(self.nb_atoms)]

    @property
    def implicit_valences(self) -> list[int]:
        """Return the implicit valence of each atom in the molecule."""
        return [self._implicit_valence(atom_idx) for atom_idx in range(self.nb_atoms)]

    def atom_mutability(self, atom_idx: int) -> bool:
        """Return the mutability of the atom at the given index."""
        atom = self.mol.GetAtomWithIdx(atom_idx)
        return atom.GetBoolProp("mutability") and not atom.GetNoImplicit()

    def set_atom_mutability(self, atom_idx: int, mutability: bool) -> None:
        """Set the mutability of the atom at the given index."""
        self.mol.GetAtomWithIdx(atom_idx).SetBoolProp("mutability", mutability)

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

    def update_representation(self, update_property_cache=False) -> None:
        """Update internal RDKit representation of the molecular graph.
        This method should be called after each action on the molecular graph.
        """
        # sanitize molecule if needed
        # print("BEFORE")
        # print("atoms", self.atoms)
        # for atom in self.mol.GetAtoms():
        #     print(
        #         atom.GetSymbol(),
        #         atom.GetFormalCharge(),
        #         atom.GetImplicitValence(),
        #         atom.GetNoImplicit(),
        #         atom.GetNumRadicalElectrons(),
        #         atom.GetBoolProp("mutability") and not atom.GetNoImplicit(),
        #     )
        if self.sanitize:
            Chem.SanitizeMol(self.mol)
            # clear computed properties
            self.mol.ClearComputedProps()

        # kekulize the molecular graph
        # https://chemistry.stackexchange.com/questions/116498/what-is-kekulization-in-rdkit
        Chem.Kekulize(self.mol, clearAromaticFlags=True)

        if update_property_cache:
            # Updating the property cache of atoms
            for atom in self.mol.GetAtoms():
                atom.UpdatePropertyCache()

            # update RDKit representation
            self.mol.UpdatePropertyCache()
        Chem.FastFindRings(self.mol)
        # print("AFTER")

        # print("atoms", self.atoms)
        # for atom in self.mol.GetAtoms():
        #     print(
        #         atom.GetSymbol(),
        #         atom.GetFormalCharge(),
        #         atom.GetImplicitValence(),
        #         atom.GetNoImplicit(),
        #         atom.GetNumRadicalElectrons(),
        #         atom.GetBoolProp("mutability") and not atom.GetNoImplicit(),
        #     )

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

    def atoms_bonded_to(self, atom_idx: int) -> list[int]:
        """Return the list of atoms bonded to the atom at the given index.

        Args:
            atom_idx (int):
                atom index

        Returns:
            list[int]: list of atoms bonded to the atom at the given index
        """
        return [
            bond.GetOtherAtomIdx(atom_idx)
            for bond in self.mol.GetAtomWithIdx(atom_idx).GetBonds()
        ]

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

    def _implicit_valence(self, atom_idx: int) -> int:
        """Return the implicit valence of the atom at the given index.

        Args:
            atom_idx (int):
                atom index

        Returns:
            int: implicit valence of the atom at the given index
        """
        atom = self.mol.GetAtomWithIdx(atom_idx)
        atom.UpdatePropertyCache()
        valence: int = atom.GetImplicitValence()
        if atom.GetSymbol() == "S":
            total_valence: int = atom.GetTotalValence()
            if total_valence in (2, 4):
                valence += 2
        return valence

    def draw(
        self,
        at_idx: bool = False,
        show: bool = True,
        size: int = 200,
        write_to_path: str = "",
    ) -> None:
        """
        Drawing the molecule
        """
        mol = self.mol.GetMol()
        atoms = mol.GetNumAtoms()

        # Setting the ids as a property if requested
        if at_idx:
            for idx in range(atoms):
                mol.GetAtomWithIdx(idx).SetProp(
                    "molAtomMapNumber", str(mol.GetAtomWithIdx(idx).GetIdx())
                )

        # Computing coordinates and making sure the properties are computed
        Compute2DCoords(mol)
        mol.UpdatePropertyCache()

        # Drawing the molecule
        dr = Draw.rdMolDraw2D.MolDraw2DCairo(size, size)
        opts = dr.drawOptions()

        # Transparent background if not writing to file
        if not write_to_path:
            opts.clearBackground = False

        dr.DrawMolecule(mol)
        dr.FinishDrawing()

        # Loading the molecule as a PIL object
        bytes_images = dr.GetDrawingText()
        image = Image.open(io.BytesIO(bytes_images))

        if show:
            image.show()

        if write_to_path:
            # Creating directories if they don't exist
            os.makedirs(os.path.dirname(write_to_path), exist_ok=True)

            # Writing image to disk
            image.save(write_to_path, "PNG")

        # return image

    ############################################################
    #                         ACTIONS                          #
    ############################################################

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
