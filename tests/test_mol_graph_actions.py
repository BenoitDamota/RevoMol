"""
Tests for the molecular graph actions.
"""

from typing import Type


from evomol.representation.molecular_graph import MolecularGraph
from evomol.representation.molecular_graph_actions import (
    AddAtomMolGraph,
    ChangeBondMolGraph,
    CutAtomMolGraph,
    InsertCarbonMolGraph,
    MoveGroupMolGraph,
    RemoveAtomMolGraph,
    RemoveGroupMolGraph,
    SubstituteAtomMolGraph,
)

from evomol.representation.molecule import Molecule, Action
from evomol.representation.smiles import SMILES


def set_parameters() -> None:
    """set the parameters for the tests."""
    Molecule.id_representation_class = SMILES
    Molecule.representations_class = [MolecularGraph]
    Molecule.max_heavy_atoms = 38
    Molecule.accepted_atoms = ["C", "O", "N", "F"]
    SubstituteAtomMolGraph.init_accepted_substitutions_from_accepted_atoms()


def check_actions(
    to_check: list[tuple[Molecule, list[Molecule]]],
    action_type: Type[Action],
) -> None:
    """Check the actions for the given molecules.

    Args:
        to_check (list[tuple[Molecule, list[Molecule]]]): for each molecule, the new possible molecules.
        action_type (Type[Action]): the type of action to check.
    """
    print(f"\nCheck {action_type.__name__}")
    for mol, check_mols in to_check:
        print(f"\n{action_type.__name__} for {mol}:")
        checked_mols = [False for _ in check_mols]
        actions = action_type.list_actions(mol)
        new_mols = []
        for action in actions:
            print(action, "->", end=" ")
            new_mol = action.apply()
            print(new_mol)
            new_mols.append(new_mol)
        print(f"Expected: {check_mols}")
        print(f"New mols: {new_mols}")
        assert len(actions) == len(check_mols)

        for new_mol in new_mols:
            checked = False
            for i in range(len(check_mols)):
                if check_mols[i] == new_mol and not checked_mols[i]:
                    checked_mols[i] = True
                    checked = True
                    break
            assert checked, f"{action}->{new_mol} not in {check_mols}"


def test_add_atom():
    """test the add atom action."""
    set_parameters()

    to_check = [
        (
            Molecule(""),
            [
                Molecule("C"),
                Molecule("O"),
                Molecule("N"),
                Molecule("F"),
            ],
        ),
        (
            Molecule("C"),
            [
                Molecule("CC"),
                Molecule("CO"),
                Molecule("CN"),
                Molecule("CF"),
            ],
        ),
        (
            Molecule("[CH3]"),
            [],
        ),
        (
            Molecule("[OH2]"),
            [],
        ),
        (
            Molecule("CS(=O)=O"),
            [
                Molecule("C([SH](=O)=O)C"),
                Molecule("C([SH](=O)=O)O"),
                Molecule("C([SH](=O)=O)N"),
                Molecule("C([SH](=O)=O)F"),
                Molecule("CS(=O)(=O)C"),
                Molecule("CS(=O)(=O)O"),
                Molecule("CS(=O)(=O)N"),
                Molecule("CS(=O)(=O)F"),
            ],
        ),
        (
            Molecule("[CH3]O"),
            [
                Molecule("COC"),
                Molecule("COO"),
                Molecule("CON"),
                Molecule("COF"),
            ],
        ),
        (
            Molecule("C[CH2]"),
            [
                Molecule("C([CH2])C"),
                Molecule("C([CH2])O"),
                Molecule("C([CH2])N"),
                Molecule("C([CH2])F"),
            ],
        ),
        (
            Molecule("C[C@@H](N)O"),
            [
                Molecule("C([C@@H](N)O)C"),
                Molecule("C([C@@H](N)O)O"),
                Molecule("C([C@@H](N)O)N"),
                Molecule("C([C@@H](N)O)F"),
                Molecule("C[C@@H](NC)O"),
                Molecule("C[C@@H](NO)O"),
                Molecule("C[C@@H](NN)O"),
                Molecule("C[C@@H](NF)O"),
                Molecule("C[C@@H](N)OC"),
                Molecule("C[C@@H](N)OO"),
                Molecule("C[C@@H](N)ON"),
                Molecule("C[C@@H](N)OF"),
            ],
        ),
        (
            Molecule("[N+]C"),
            [
                Molecule("[N+]CC"),
                Molecule("[N+]CO"),
                Molecule("[N+]CN"),
                Molecule("[N+]CF"),
            ],
        ),
        (
            Molecule("O=S=O"),
            [
                Molecule("O=[SH](=O)C"),
                Molecule("O=[SH](=O)O"),
                Molecule("O=[SH](=O)N"),
                Molecule("O=[SH](=O)F"),
            ],
        ),
    ]

    check_actions(to_check, AddAtomMolGraph)


def test_remove_atom():
    """test the remove atom action."""
    set_parameters()

    to_check = [
        (
            Molecule(""),
            [],
        ),
        (
            Molecule("[CH]"),
            [],
        ),
        (
            Molecule("[N+]"),
            [],
        ),
        (
            Molecule("CC"),
            [
                Molecule("C"),
                Molecule("C"),
            ],
        ),
        (
            Molecule("C[CH]"),
            [],
        ),
        (
            Molecule("CC[CH2]"),
            [
                Molecule("C[CH2]"),
            ],
        ),
        (
            Molecule("C[CH2]C"),
            [],
        ),
        (
            Molecule("c1ccc(cc1)[N+](=O)[O-]"),
            [
                Molecule("CC=C(C=C)[N+](=O)[O-]"),
                Molecule("CC=CC(=C)[N+](=O)[O-]"),
                Molecule("C(=C)C=CC[N+](=O)[O-]"),
                Molecule("C(=CC=C[N+](=O)[O-])C"),
                Molecule("C=CC=C(C)[N+](=O)[O-]"),
            ],
        ),
        (
            Molecule("C[C@@H](N)O"),
            [],
        ),
    ]

    check_actions(to_check, RemoveAtomMolGraph)


def test_change_bond_atom():
    """test the change bond action."""
    set_parameters()

    to_check = [
        (
            Molecule(""),
            [],
        ),
        (
            Molecule("[CH]"),
            [],
        ),
        (
            Molecule("C"),
            [],
        ),
        (
            Molecule("CC"),
            [
                Molecule("C=C"),
                Molecule("C#C"),
            ],
        ),
        (
            Molecule("CCN"),
            [
                Molecule("C=CN"),
                Molecule("C#CN"),
                Molecule("CC=N"),
                Molecule("CC#N"),
                Molecule("C1CN1"),
                Molecule("C1=NC1"),
            ],
        ),
        (
            Molecule("C1C[CH2]1"),
            [
                Molecule("CCC"),
                Molecule("C1=CC1"),
                Molecule("C1#CC1"),
            ],
        ),
        (
            Molecule("C1C[N+]1"),
            [
                Molecule("C[N+]C"),
                Molecule("C1=C[N+]1"),
                Molecule("C1#C[N+]1"),
            ],
        ),
        (
            Molecule("C[C@@H](N)O"),
            [
                Molecule("C1[C@H](O)N1"),
                Molecule("C1=N[C@H]1O"),
                Molecule("C1[C@@H](N)O1"),
                Molecule("C[C@H]1NO1"),
            ],
        ),
    ]

    check_actions(to_check, ChangeBondMolGraph)


# # TODO
# def test_substitute_atom():
#     """test the substitute atom action."""
#     set_parameters()

#     to_check = [
#         (
#             Molecule("C[C@@H](N)C(=O)O"),
#             [
#                 Molecule("N[C@H](O)C(=O)O"),
#                 Molecule("NC(N)C(=O)O"),
#                 Molecule("N[C@H](F)C(=O)O"),
#                 Molecule("CC(C)C(=O)O"),
#                 Molecule("C[C@@H](O)C(=O)O"),
#                 Molecule("C[C@@H](F)C(=O)O"),
#                 Molecule("C=C(O)[C@@H](C)N"),
#                 Molecule("C[C@@H](N)C(=N)O"),
#                 Molecule("CC(=O)[C@@H](C)N"),
#                 Molecule("C[C@@H](N)C(N)=O"),
#                 Molecule("C[C@@H](N)C(=O)F"),
#             ],
#         ),
#         (
#             Molecule(""),
#             [],
#         ),
#     ]

#     check_actions(to_check, SubstituteAtomMolGraph)


def test_insert_carbon():
    """test the insert carbon action."""
    set_parameters()

    to_check = [
        (
            Molecule("C"),
            [],
        ),
        (
            Molecule("CC"),
            [
                Molecule("CCC"),
            ],
        ),
        (
            Molecule("C=C"),
            [
                Molecule("CCC"),
            ],
        ),
        (
            Molecule("C#C"),
            [
                Molecule("CCC"),
            ],
        ),
        (
            Molecule("C=[N+]"),
            [
                Molecule("CC=[N+]"),
            ],
        ),
        (
            Molecule("[N+]=[N+]"),
            [
                Molecule("[N+]=C=[N+]"),
            ],
        ),
        (
            Molecule("C[NH+]([O-])O"),
            [
                Molecule("CC[NH+]([O-])O"),
                Molecule("C[NH+](O)C[O-]"),
                Molecule("C[NH+]([O-])CO"),
            ],
        ),
        (
            Molecule("C[N+](=O)[O-]"),
            [
                Molecule("CC[N+](=O)[O-]"),
                Molecule("C[N+]([O-])=CO"),
                Molecule("C[N+](=O)C[O-]"),
            ],
        ),
        (
            Molecule("[N+]#[N+]"),
            [],
        ),
        (
            Molecule("C#[N+]"),
            [
                Molecule("CC#[N+]"),
            ],
        ),
    ]

    check_actions(to_check, InsertCarbonMolGraph)


def test_cut_atom():
    """test the cut atom action."""
    set_parameters()

    to_check = [
        (
            Molecule("CCC"),
            [
                # cut middle C
                Molecule("CC"),
            ],
        ),
        (
            Molecule("C=C=C"),
            [
                # cut middle C then simple bond
                Molecule("CC"),
            ],
        ),
        (
            Molecule("[C+]#CC"),
            [
                # cut middle C then triple bond
                Molecule("[C+]#C"),
            ],
        ),
        (
            Molecule("O[CH2]O"),
            [
                # no cut (C is not mutable)
            ],
        ),
        (
            Molecule("O[CH2]OO"),
            [
                # one cut on 2nd O (lost the [CH2])
                Molecule("OCO"),
            ],
        ),
        (
            Molecule("O[CH]N=O"),
            [
                # one cut on N
                Molecule("O[CH]O"),
            ],
        ),
        (
            Molecule("O[C+]=C[C+]O"),
            [
                # not cut (C+ are not mutable and not the same bond)
            ],
        ),
        (
            Molecule("O[C+]=C=[C+]O"),
            [
                # one cut middle C with a double bond after
                Molecule("O[C+]=[C+]O"),
            ],
        ),
        (
            Molecule("O[C+]=CC"),
            [
                # one cut on middle C with a double bond after
                Molecule("O[C+]=C"),
            ],
        ),
        (
            Molecule("O[C+]C=C"),
            [
                # one cut on middle C with a simple bond after
                Molecule("O[C+]C"),
            ],
        ),
        (
            Molecule("OC=C[C+]"),
            [
                # one cut on middle C with a simple bond after
                Molecule("OC[C+]"),
                # one cut 1st C
                Molecule("OC[C+]"),
            ],
        ),
        (
            Molecule("OC=C[C]"),
            [
                Molecule("OC[C]"),
                Molecule("OC[C]"),
            ],
        ),
        (
            Molecule("OC=C[CH3]"),
            [
                Molecule("OCC"),
                Molecule("OCC"),
            ],
        ),
    ]

    check_actions(to_check, CutAtomMolGraph)


def test_move_group():
    """test the move group action."""
    set_parameters()

    to_check = [
        (
            Molecule("C1=CC=C[CH]=C1CC"),
            [
                Molecule("C1(CC)=CC=CC=C1"),
                Molecule("C1=C(CC)C=CC=C1"),
                Molecule("C1=CC(CC)=CC=C1"),
                Molecule("C1=CC=C(CC)C=C1"),
                Molecule("C1=CC=CC=C1CC"),
                Molecule("C1(C)=CC=CC=C1C"),
                Molecule("C1=C(C)C=CC=C1C"),
                Molecule("C1=CC(C)=CC=C1C"),
                Molecule("C1=CC=C(C)C=C1C"),
            ],
        ),
        (
            Molecule("C1=CC=C[CH]=C1C1=CC=C[CH]=C1"),
            [
                # 8 possibilities and not 10 as [CH] is not mutable
                Molecule("C1(C2=CC=CC=C2)=CC=CC=C1"),
                Molecule("C1=C(C2=CC=CC=C2)C=CC=C1"),
                Molecule("C1=CC(C2=CC=CC=C2)=CC=C1"),
                Molecule("C1=CC=C(C2=CC=CC=C2)C=C1"),
                Molecule("C1=CC=CC=C1C1=CC=CC=C1"),
                Molecule("C1=CC=CC=C1C1=CC=CC=C1"),
                Molecule("C1=CC=CC=C1C1=CC=CC=C1"),
                Molecule("C1=CC=CC=C1C1=CC=CC=C1"),
            ],
        ),
        (
            Molecule("c1ncncc1CCO"),
            [
                Molecule("C1(CCO)=NC=NC=C1"),
                Molecule("C1=NC(CCO)=NC=C1"),
                Molecule("C1=NC=NC(CCO)=C1"),
                Molecule("C1=NC=NC=C1C(C)O"),
                Molecule("C1=NC=NC=C1OCC"),
                Molecule("C1(CO)=NC=NC=C1C"),
                Molecule("C1=NC(CO)=NC=C1C"),
                Molecule("C1=NC=NC(CO)=C1C"),
                Molecule("C1=NC=NC=C1COC"),
                Molecule("C1(O)=NC=NC=C1CC"),
                Molecule("C1=NC(O)=NC=C1CC"),
                Molecule("C1=NC=NC(O)=C1CC"),
                Molecule("C1=NC=NC=C1C(C)O"),
            ],
        ),
        (
            Molecule("C1NC1C1SC1"),
            [
                Molecule("C1(C2SC2)NC1"),
                Molecule("C1N(C2SC2)C1"),
                Molecule("C1NC1[SH]1CC1"),
                Molecule("C1NC1C1CS1"),
            ],
        ),
        (
            Molecule("C1[N+]=C1C1SC1"),
            [
                # can't connect to N+ as it is not mutable
                Molecule("C1(C2SC2)[N+]=C1"),
                Molecule("C1[N+]=C1[SH]1CC1"),
                Molecule("C1[N+]=C1C1CS1"),
            ],
        ),
        (
            Molecule("C1C[N+]1=C1SC1"),
            [
                # can connect only N+ to the other side as it is not mutable
                Molecule("C1C[N+]1=S1CC1"),
                Molecule("C1C[N+]1=C1CS1"),
            ],
        ),
        (
            Molecule("C1C[N+]1=[N+]1SC1"),
            [
                # can't do anything as N+ are not mutable
            ],
        ),
    ]

    check_actions(to_check, MoveGroupMolGraph)


# TODO
def test_remove_group():
    """test the remove group action."""
    set_parameters()

    to_check = [
        (
            Molecule(""),
            [],
        ),
    ]

    check_actions(to_check, RemoveGroupMolGraph)
