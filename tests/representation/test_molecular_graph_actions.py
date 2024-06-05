"""
Tests for the molecular graph actions.
"""

import pytest

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
from evomol.representation.molecule import Action, Molecule
from evomol.representation.smiles import SMILES


@pytest.fixture(autouse=True, scope="module")
def setup_parameters() -> None:
    """set the parameters for the tests."""
    Molecule.id_representation_class = SMILES
    Molecule.representations_class = [MolecularGraph]
    Molecule.max_heavy_atoms = 38
    Molecule.accepted_atoms = ["C", "O", "N", "F"]
    SubstituteAtomMolGraph.init_accepted_substitutions_from_accepted_atoms()


def check_actions_smiles(
    smiles: str, check_smiles: list[str], action_type: type[Action]
) -> None:
    """Check the actions for the given smiles."""
    mol = Molecule(smiles)
    check_mols: list[Molecule] = [Molecule(smiles) for smiles in check_smiles]
    print(f"{action_type.__name__} for {mol}:")
    checked_mols = [False for _ in check_mols]
    actions = action_type.list_actions(mol)
    new_mols = []
    for action in actions:
        print(action, "->", end=" ")
        new_mol = action.apply()
        print(new_mol)
        new_mols.append(new_mol)
    print(f"Expected: {check_mols}")
    print(f"Got     : {new_mols}")
    assert len(actions) == len(check_mols)

    for new_mol in new_mols:
        checked = False
        for i, check_mol in enumerate(check_mols):
            if check_mol == new_mol and not checked_mols[i]:
                checked_mols[i] = True
                checked = True
                break
        assert (
            checked
        ), f"{smiles} - {action_type.__name__} : {new_mol} not in {check_mols}"


########################################
#                AddAtom
########################################
@pytest.mark.parametrize(
    ("smiles", "check_smiles"),
    (
        pytest.param(
            "",
            [
                "C",
                "O",
                "N",
                "F",
            ],
            id="AddAtom on empty smiles",
        ),
        pytest.param(
            "C",
            [
                "CC",
                "CO",
                "CN",
                "CF",
            ],
            id="AddAtom on C",
        ),
        pytest.param(
            "[CH3]",
            [],
            id="AddAtom on [CH3]",
        ),
        pytest.param(
            "[OH2]",
            [],
            id="AddAtom on [OH2]",
        ),
        pytest.param(
            "CS(=O)=O",
            [
                "C([SH](=O)=O)C",
                "C([SH](=O)=O)O",
                "C([SH](=O)=O)N",
                "C([SH](=O)=O)F",
                "CS(=O)(=O)C",
                "CS(=O)(=O)O",
                "CS(=O)(=O)N",
                "CS(=O)(=O)F",
            ],
            id="AddAtom on CS(=O)=O",
        ),
        pytest.param(
            "[CH3]O",
            [
                "COC",
                "COO",
                "CON",
                "COF",
            ],
            id="AddAtom on [CH3]O",
        ),
        pytest.param(
            "C[CH2]",
            [
                "C([CH2])C",
                "C([CH2])O",
                "C([CH2])N",
                "C([CH2])F",
            ],
            id="AddAtom on C[CH2]",
        ),
        pytest.param(
            "C[C@@H](N)O",
            [
                "C([C@@H](N)O)C",
                "C([C@@H](N)O)O",
                "C([C@@H](N)O)N",
                "C([C@@H](N)O)F",
                "C[C@@H](NC)O",
                "C[C@@H](NO)O",
                "C[C@@H](NN)O",
                "C[C@@H](NF)O",
                "C[C@@H](N)OC",
                "C[C@@H](N)OO",
                "C[C@@H](N)ON",
                "C[C@@H](N)OF",
            ],
            id="AddAtom on C[C@@H](N)O",
        ),
        pytest.param(
            "[N+]C",
            [
                "[N+]CC",
                "[N+]CO",
                "[N+]CN",
                "[N+]CF",
            ],
            id="AddAtom on [N+]C",
        ),
        pytest.param(
            "O=S=O",
            [
                "O=[SH](=O)C",
                "O=[SH](=O)O",
                "O=[SH](=O)N",
                "O=[SH](=O)F",
            ],
            id="AddAtom on O=S=O",
        ),
    ),
)
def test_actions_add_atom(smiles: str, check_smiles: list[str]) -> None:
    """Test the AddAtom action."""
    check_actions_smiles(smiles, check_smiles, AddAtomMolGraph)


########################################
#               RemoveAtom
########################################
@pytest.mark.parametrize(
    ("smiles", "check_smiles"),
    (
        pytest.param(
            "",
            [],
            id="RemoveAtom on empty smiles",
        ),
        pytest.param(
            "C",
            [""],
            id="RemoveAtom on C",
        ),
        pytest.param(
            "CC",
            [
                "C",
                "C",
            ],
            id="RemoveAtom on CC",
        ),
        pytest.param(
            "[CH]",
            [],
            id="RemoveAtom on [CH]",
        ),
        pytest.param(
            "[N+]",
            [],
            id="RemoveAtom on [N+]",
        ),
        pytest.param(
            "CC[CH3]",
            [
                "CC",
                "CC",
            ],
            id="RemoveAtom on CC[CH3]",
        ),
        # pytest.param(
        #     "CC[CH2]",
        #     [
        #         "C[CH2]",
        #         "CC",
        #     ],
        #     id="RemoveAtom on CC[CH2]",
        # ),
        # pytest.param(
        #     "C[CH2]C",
        #     [],
        #     id="RemoveAtom on C[CH2]C",
        # ),
        pytest.param(
            "c1ccc(cc1)[N+](=O)[O-]",
            [
                "CC=C(C=C)[N+](=O)[O-]",
                "CC=CC(=C)[N+](=O)[O-]",
                "C(=C)C=CC[N+](=O)[O-]",
                "C(=CC=C[N+](=O)[O-])C",
                "C=CC=C(C)[N+](=O)[O-]",
            ],
            id="RemoveAtom on c1ccc(cc1)[N+](=O)[O-]",
        ),
        # pytest.param( # TODO
        #     "C[C@@H](N)O",
        #     [],
        #     id="RemoveAtom on C[C@@H](N)O",
        # ),
    ),
)
def test_actions_remove_atom(smiles: str, check_smiles: list[str]) -> None:
    """Test the RemoveAtom action."""
    check_actions_smiles(smiles, check_smiles, RemoveAtomMolGraph)


########################################
#               ChangeBond
########################################
@pytest.mark.parametrize(
    ("smiles", "check_smiles"),
    (
        pytest.param(
            "",
            [],
            id="ChangeBond on empty smiles",
        ),
        pytest.param(
            "[CH]",
            [],
            id="ChangeBond on [CH]",
        ),
        pytest.param(
            "C",
            [],
            id="ChangeBond on C",
        ),
        pytest.param(
            "CC",
            [
                "C=C",
                "C#C",
            ],
            id="ChangeBond on CC",
        ),
        pytest.param(
            "CCN",
            [
                "C=CN",
                "C#CN",
                "CC=N",
                "CC#N",
                "C1CN1",
                "C1=NC1",
            ],
            id="ChangeBond on CCN",
        ),
        pytest.param(
            "C1C[CH2]1",
            [
                "CCC",
                "C1=CC1",
                "C1#CC1",
                "CC[CH2]",
                "C(C)[CH2]",
            ],
            id="ChangeBond on C1C[CH2]1",
        ),
        pytest.param(
            "C1C[N+]1",
            [
                "C[N+]C",
                "C1=C[N+]1",
                "C1#C[N+]1",
            ],
            id="ChangeBond on C1C[N+]1",
        ),
        pytest.param(
            "C[C@@H](N)O",
            [
                "C1[C@H](O)N1",
                "C1=N[C@H]1O",
                "C1[C@@H](N)O1",
                "C[C@H]1NO1",
            ],
            id="ChangeBond on C[C@@H](N)O",
        ),
    ),
)
def test_actions_change_bond(smiles: str, check_smiles: list[str]) -> None:
    """Test the ChangeBond action."""
    check_actions_smiles(smiles, check_smiles, ChangeBondMolGraph)


########################################
#             SubstituteAtom
########################################
@pytest.mark.parametrize(
    ("smiles", "check_smiles"),
    (
        pytest.param(
            "",
            [],
            id="SubstituteAtom on empty smiles",
        ),
        pytest.param(
            "C",
            [
                "O",
                "N",
                "F",
            ],
            id="SubstituteAtom on C",
        ),
        pytest.param(
            "[CH3]",
            [],
            id="SubstituteAtom on [CH3]",
        ),
        pytest.param(
            "[N+]=C",
            [
                "[N+]=O",
                "[N+]=N",
            ],
            id="SubstituteAtom on [N+]=C",
        ),
        pytest.param(
            "[N+]=[N+]",
            [],
            id="SubstituteAtom on [N+]=[N+]",
        ),
        pytest.param(
            "C[C@@H](N)C(=O)O",
            [
                "O[C@@H](N)C(=O)O",
                "NC(N)C(=O)O",
                "F[C@@H](N)C(=O)O",
                "CC(C)C(=O)O",
                "C[C@@H](O)C(=O)O",
                "C[C@@H](F)C(=O)O",
                "C[C@@H](N)C(=C)O",
                "C[C@@H](N)C(=N)O",
                "C[C@@H](N)C(=O)C",
                "C[C@@H](N)C(=O)N",
                "C[C@@H](N)C(=O)F",
            ],
            id="SubstituteAtom on C[C@@H](N)C(=O)O",
        ),
    ),
)
def test_actions_substitute(smiles: str, check_smiles: list[str]) -> None:
    """Test the SubstituteAtom action."""
    check_actions_smiles(smiles, check_smiles, SubstituteAtomMolGraph)


########################################
#              InsertCarbon
########################################
@pytest.mark.parametrize(
    ("smiles", "check_smiles"),
    (
        pytest.param(
            "C",
            [],
            id="InsertCarbon on C",
        ),
        pytest.param(
            "CC",
            [
                "CCC",
            ],
            id="InsertCarbon on CC",
        ),
        pytest.param(
            "C=C",
            [
                "CCC",
            ],
            id="InsertCarbon on C=C",
        ),
        pytest.param(
            "C#C",
            [
                "CCC",
            ],
            id="InsertCarbon on C#C",
        ),
        pytest.param(
            "C=[N+]",
            [
                "CC=[N+]",
            ],
            id="InsertCarbon on C=[N+]",
        ),
        pytest.param(
            "[N+]=[N+]",
            [
                "[N+]=C=[N+]",
            ],
            id="InsertCarbon on [N+]=[N+]",
        ),
        pytest.param(
            "C[NH+]([O-])O",
            [
                "CC[NH+]([O-])O",
                "C[NH+](O)C[O-]",
                "C[NH+]([O-])CO",
            ],
            id="InsertCarbon on C[NH+]([O-])O",
        ),
        pytest.param(
            "C[N+](=O)[O-]",
            [
                "CC[N+](=O)[O-]",
                "C[N+]([O-])=CO",
                "C[N+](=O)C[O-]",
            ],
            id="InsertCarbon on C[N+](=O)[O-]",
        ),
        pytest.param(
            "[N+]#[N+]",
            [],
            id="InsertCarbon on [N+]#[N+]",
        ),
        pytest.param(
            "C#[N+]",
            [
                "CC#[N+]",
            ],
            id="InsertCarbon on C#[N+]",
        ),
    ),
)
def test_actions_insert_carbon(smiles: str, check_smiles: list[str]) -> None:
    """Test the InsertCarbon action."""
    check_actions_smiles(smiles, check_smiles, InsertCarbonMolGraph)


########################################
#               CutAtom
########################################
@pytest.mark.parametrize(
    ("smiles", "check_smiles"),
    (
        pytest.param(
            "CCC",
            [
                "CC",
            ],
            id="CutAtom on CCC",
        ),
        pytest.param(
            "C=C=C",
            [
                "CC",
            ],
            id="CutAtom on C=C=C",
        ),
        pytest.param(
            "[C+]#CC",
            [
                "[C+]#C",
            ],
            id="CutAtom on [C+]#CC",
        ),
        pytest.param(
            "O[CH2]OO",
            [
                "OOO",
                "OCO",
            ],
            id="CutAtom on O[CH2]OO",
        ),
        pytest.param(
            "O[CH]N=O",
            [
                "O[CH]O",
            ],
            id="CutAtom on O[CH]N=O",
        ),
        pytest.param(
            "O[C+]=C[C+]O",
            [],
            id="CutAtom on O[C+]=C[C+]O",
        ),
        pytest.param(
            "O[C+]=C=[C+]O",
            [
                "O[C+]=[C+]O",
            ],
            id="CutAtom on O[C+]=C=[C+]O",
        ),
        pytest.param(
            "O[C+]=CC",
            [
                "O[C+]=C",
            ],
            id="CutAtom on O[C+]=CC",
        ),
        pytest.param(
            "O[C+]C=C",
            [
                "O[C+]C",
            ],
            id="CutAtom on O[C+]C=C",
        ),
        pytest.param(
            "OC=C[C+]",
            [
                "OC[C+]",
                "OC[C+]",
            ],
            id="CutAtom on OC=C[C+]",
        ),
        pytest.param(
            "OC=C[C]",
            [
                "OC[C]",
                "OC[C]",
            ],
            id="CutAtom on OC=C[C]",
        ),
        pytest.param(
            "OC=C[CH3]",
            [
                "OCC",
                "OCC",
            ],
            id="CutAtom on OC=C[CH3]",
        ),
    ),
)
def test_actions_cut_atom(smiles: str, check_smiles: list[str]) -> None:
    """Test the CutAtom action."""
    check_actions_smiles(smiles, check_smiles, CutAtomMolGraph)


########################################
#               MoveGroup
########################################
@pytest.mark.parametrize(
    ("smiles", "check_smiles"),
    (
        pytest.param(
            "C1=CC=C[CH]=C1CC",
            [
                "C1(CC)=CC=CC=C1",
                "C1=C(CC)C=CC=C1",
                "C1=CC(CC)=CC=C1",
                "C1=CC=C(CC)C=C1",
                "C1=CC=CC=C1CC",
                "C1(C)=CC=CC=C1C",
                "C1=C(C)C=CC=C1C",
                "C1=CC(C)=CC=C1C",
                "C1=CC=C(C)C=C1C",
            ],
            id="MoveGroup on C1=CC=C[CH]=C1CC",
        ),
        pytest.param(
            "C1=CC=C[CH]=C1C1=CC=C[CH]=C1",
            [
                "C1(C2=CC=CC=C2)=CC=CC=C1",
                "C1=C(C2=CC=CC=C2)C=CC=C1",
                "C1=CC(C2=CC=CC=C2)=CC=C1",
                "C1=CC=C(C2=CC=CC=C2)C=C1",
                "C1=CC=CC=C1C1=CC=CC=C1",
                "C1=CC=CC=C1C1=CC=CC=C1",
                "C1=CC=CC=C1C1=CC=CC=C1",
                "C1=CC=CC=C1C1=CC=CC=C1",
            ],
            id="MoveGroup on C1=CC=C[CH]=C1C1=CC=C[CH]=C1",
        ),
        pytest.param(
            "c1ncncc1CCO",
            [
                "C1(CCO)=NC=NC=C1",
                "C1=NC(CCO)=NC=C1",
                "C1=NC=NC(CCO)=C1",
                "C1=NC=NC=C1C(C)O",
                "C1=NC=NC=C1OCC",
                "C1(CO)=NC=NC=C1C",
                "C1=NC(CO)=NC=C1C",
                "C1=NC=NC(CO)=C1C",
                "C1=NC=NC=C1COC",
                "C1(O)=NC=NC=C1CC",
                "C1=NC(O)=NC=C1CC",
                "C1=NC=NC(O)=C1CC",
                "C1=NC=NC=C1C(C)O",
            ],
            id="MoveGroup on c1ncncc1CCO",
        ),
        pytest.param(
            "C1NC1C1SC1",
            [
                "C1(C2SC2)NC1",
                "C1N(C2SC2)C1",
                "C1NC1[SH]1CC1",
                "C1NC1C1CS1",
            ],
            id="MoveGroup on C1NC1C1SC1",
        ),
        pytest.param(
            "C1[N+]=C1C1SC1",
            [
                "C1(C2SC2)[N+]=C1",
                "C1[N+]=C1[SH]1CC1",
                "C1[N+]=C1C1CS1",
            ],
            id="MoveGroup on C1[N+]=C1C1SC1",
        ),
        pytest.param(
            "C1C[N+]1=C1SC1",
            [
                "C1C[N+]1=S1CC1",
                "C1C[N+]1=C1CS1",
            ],
            id="MoveGroup on C1C[N+]1=C1SC1",
        ),
        pytest.param(
            "C1C[N+]1=[N+]1SC1",
            [],
            id="MoveGroup on C1C[N+]1=[N+]1SC1",
        ),
    ),
)
def test_actions_move_group(smiles: str, check_smiles: list[str]) -> None:
    """Test the MoveGroup action."""
    check_actions_smiles(smiles, check_smiles, MoveGroupMolGraph)


########################################
#               RemoveGroup
########################################
@pytest.mark.parametrize(
    ("smiles", "check_smiles"),
    (
        pytest.param(
            "",
            [],
            id="RemoveGroup on empty smiles",
        ),
        pytest.param(
            "C1=CC=C[CH]=C1CC",
            ["C1=CC=CC=C1", "C1=CC=CC=C1C", "CC", "C"],
            id="RemoveGroup on C1=CC=C[CH]=C1CC",
        ),
        pytest.param(
            "c1ccc[S+]1[N+]1[CH]ccocc1",
            [],
            id="RemoveGroup on c1ccc[S+]1[N+]1[CH]ccocc1",
        ),
    ),
)
def test_actions_remove_group(smiles: str, check_smiles: list[str]) -> None:
    """Test the RemoveGroup action."""
    check_actions_smiles(smiles, check_smiles, RemoveGroupMolGraph)


# @pytest.mark.usefixtures("setup_parameters")
# @pytest.mark.parametrize(
#     ("input_a", "input_b", "input_c"),
#     (
#         ("a", "b", "c"),
#         ([], None, 0),
#     ),
# )
# def test_error(input_a, input_b, input_c):
#     with pytest.raises(AttributeError):
#         check_actions_mol(input_a, input_b, input_c)
