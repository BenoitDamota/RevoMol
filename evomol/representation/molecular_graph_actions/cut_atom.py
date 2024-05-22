"""
Cut an atom from the molecular graph.
"""

from copy import copy
from typing_extensions import override

from evomol.representation.molecular_graph import MolecularGraph
from evomol.representation.molecule import Action, Molecule


class CutAtomMolGraph(Action):
    """
    Cut an atom from the molecular graph.
    """

    def __init__(
        self,
        molecule: Molecule,
        atom_to_cut: int,
        atom1_to_bond: int,
        atom2_to_bond: int,
        bond_type: int,
    ):
        super().__init__(molecule)
        # index of the atom to cut
        self.atom_to_cut: int = atom_to_cut
        # index of the atom to bond to the new atom
        self.atom1_to_bond: int = atom1_to_bond
        self.atom2_to_bond: int = atom2_to_bond
        self.bond_type: int = bond_type

    @override
    def apply(self) -> Molecule:
        mol_graph: MolecularGraph = self.molecule.get_representation(MolecularGraph)
        assert mol_graph is not None

        new_mol_graph: MolecularGraph = copy(mol_graph)

        # remove bonds
        new_mol_graph.set_bond(self.atom_to_cut, self.atom1_to_bond, 0)
        new_mol_graph.set_bond(self.atom_to_cut, self.atom2_to_bond, 0)

        # create new bond
        new_mol_graph.set_bond(self.atom1_to_bond, self.atom2_to_bond, self.bond_type)

        new_mol_graph.remove_atom(self.atom_to_cut)

        try:
            new_mol_graph.update_representation()
        except Exception as e:
            print("Cut caused an error.")
            raise e

        return Molecule(new_mol_graph.canonical_smiles)

    def __repr__(self) -> str:
        return (
            "CutAtomMolGraph("
            f"{self.molecule}, {self.atom_to_cut}, "
            f"{self.atom1_to_bond}, {self.atom2_to_bond}, {self.bond_type})"
        )

    @override
    @classmethod
    def list_actions(cls, molecule: Molecule) -> list[Action]:
        """List possible actions to cut an atom from the molecular graph."""

        mol_graph: MolecularGraph = molecule.get_representation(MolecularGraph)

        action_list: list[Action] = []

        formal_charges = mol_graph.formal_charge_vector()

        for atom_to_cut in range(mol_graph.nb_atoms):
            # we cannot cut a non mutable atom
            if not mol_graph.atom_mutability(atom_to_cut):
                # print(f"{atom_to_cut=}")
                # print(mol_graph.atoms())
                # print(mol_graph.canonical_smiles)
                # print(mol_graph.atom_mutability(atom_to_cut))
                # atom = mol_graph.mol.GetAtomWithIdx(atom_to_cut)
                # print(
                #     atom.GetSymbol(),
                #     atom.GetFormalCharge(),
                #     atom.GetImplicitValence(),
                #     atom.GetNoImplicit(),
                #     atom.GetNumRadicalElectrons(),
                #     atom.GetBoolProp("mutability") and not atom.GetNoImplicit(),
                # )
                continue

            # if there is two atoms bonded to the atom to cut
            bonds_to = [
                atom2_idx
                for atom2_idx in range(mol_graph.nb_atoms)
                if mol_graph.bond_type_num(atom_to_cut, atom2_idx) > 0
            ]
            if len(bonds_to) != 2:
                continue

            atom_1 = bonds_to[0]
            atom_2 = bonds_to[1]

            if mol_graph.bond_type_num(atom_1, atom_2) != 0:
                continue

            # atom_1 to atom_to_cut
            bond_1 = mol_graph.bond_type_num(atom_to_cut, atom_1)

            # atom_2 to atom_to_cut
            bond_2 = mol_graph.bond_type_num(atom_to_cut, atom_2)

            mutability_1 = mol_graph.atom_mutability(atom_1)
            mutability_2 = mol_graph.atom_mutability(atom_2)

            if not mutability_1 and not mutability_2:
                if bond_1 != bond_2:
                    continue
                print(
                    f"1 - {atom_to_cut=} ({mol_graph.atom_type(atom_to_cut)}) connect {atom_1=} ({mol_graph.atom_type(atom_1)}) and {atom_2=} ({mol_graph.atom_type(atom_2)}) with {bond_1}"
                )
                action_list.append(
                    CutAtomMolGraph(
                        molecule, atom_to_cut, bonds_to[0], bonds_to[1], bond_1
                    )
                )
                continue

            max_valence_1 = mol_graph.implicit_valence_vector()[atom_1] + bond_1
            max_valence_2 = mol_graph.implicit_valence_vector()[atom_2] + bond_2
            if not mutability_1 and max_valence_2 >= bond_1:
                print(
                    f"2 - {atom_to_cut=} ({mol_graph.atom_type(atom_to_cut)}) connect {atom_1=} ({mol_graph.atom_type(atom_1)}) and {atom_2=} ({mol_graph.atom_type(atom_2)}) with {bond_1}"
                )
                action_list.append(
                    CutAtomMolGraph(
                        molecule, atom_to_cut, bonds_to[0], bonds_to[1], bond_1
                    )
                )
                continue

            if not mutability_2 and max_valence_1 >= bond_2:
                print(
                    f"3 - {atom_to_cut=} ({mol_graph.atom_type(atom_to_cut)}) connect {atom_1=} ({mol_graph.atom_type(atom_1)}) and {atom_2=} ({mol_graph.atom_type(atom_2)}) with {bond_2}"
                )
                action_list.append(
                    CutAtomMolGraph(
                        molecule, atom_to_cut, bonds_to[0], bonds_to[1], bond_2
                    )
                )
                continue

            if not mutability_1 or not mutability_2:
                continue

            print(
                f"4 - {atom_to_cut=} ({mol_graph.atom_type(atom_to_cut)}) connect {atom_1=} ({mol_graph.atom_type(atom_1)}) and {atom_2=} ({mol_graph.atom_type(atom_2)}) with 1"
            )
            action_list.append(
                CutAtomMolGraph(molecule, atom_to_cut, bonds_to[0], bonds_to[1], 1)
            )

        return action_list
