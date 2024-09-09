"""
Tests for the molecular graph actions.
"""

import pytest

from evomol import default_parameters as dp
from evomol.action import Action
from evomol.action import molecular_graph as mg
from evomol.representation import Molecule


@pytest.fixture(autouse=True, scope="module")
def setup_parameters() -> None:
    """set the parameters for the tests."""
    dp.setup_default_parameters(
        accepted_atoms=["C", "O", "N", "F"],
        max_heavy_atoms=38,
    )
    dp.setup_default_action_space()


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
    assert len(actions) == len(
        check_mols
    ), f"Expected {len(check_mols)} actions, got {len(actions)}"

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
            [
                "CO",
                "OO",
                "NO",
                "OF",
            ],
            id="AddAtom on [OH2]",
        ),
        pytest.param(
            "CS(=O)=O",
            [
                "CC[SH](=O)=O",
                "O=[SH](=O)CO",
                "NC[SH](=O)=O",
                "O=[SH](=O)CF",
                "CS(C)(=O)=O",
                "CS(=O)(=O)O",
                "CS(N)(=O)=O",
                "CS(=O)(=O)F",
            ],
            id="AddAtom on CS(=O)=O",
        ),
        pytest.param(
            "[CH3]O",
            [
                "CCO",
                "OCO",
                "NCO",
                "OCF",
                "COC",
                "COO",
                "CON",
                "COF",
            ],
            id="AddAtom on [CH3]O",
        ),
        pytest.param(
            "[SH]O",
            [
                "CSO",
                "OSO",
                "NSO",
                "OSF",
                "COS",
                "OOS",
                "NOS",
                "FOS",
            ],
            id="AddAtom on [SH]O",
        ),
        pytest.param(
            "C[CH2]",
            [
                "[CH2]CC",
                "[CH2]CO",
                "[CH2]CN",
                "[CH2]CF",
            ],
            id="AddAtom on C[CH2]",
        ),
        pytest.param(
            "C[C@@H](N)O",
            [
                "CC[C@@H](N)O",
                "N[C@@H](O)CO",
                "NC[C@@H](N)O",
                "N[C@@H](O)CF",
                "CC(C)(N)O",
                "CC(N)(O)O",
                "CC(N)(N)O",
                "C[C@@](N)(O)F",
                "CN[C@H](C)O",
                "C[C@H](O)NO",
                "C[C@H](O)NN",
                "C[C@H](O)NF",
                "CO[C@@H](C)N",
                "C[C@@H](N)OO",
                "C[C@@H](N)ON",
                "C[C@@H](N)OF",
            ],
            id="AddAtom on C[C@@H](N)O",
        ),
        pytest.param(
            "[N+]C",
            [
                "CC[N+]",
                "[N+]CO",
                "[N+]CN",
                "[N+]CF",
            ],
            id="AddAtom on [N+]C",
        ),
        pytest.param(
            "O=S=O",
            [
                "C[SH](=O)=O",
                "O=[SH](=O)O",
                "N[SH](=O)=O",
                "O=[SH](=O)F",
            ],
            id="AddAtom on O=S=O",
        ),
    ),
)
def test_actions_add_atom(smiles: str, check_smiles: list[str]) -> None:
    """Test the AddAtom action."""
    check_actions_smiles(smiles, check_smiles, mg.AddAtomMG)


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
        pytest.param(
            "CC[CH2]",
            [
                "[CH2]C",
            ],
            id="RemoveAtom on CC[CH2]",
        ),
        pytest.param(
            "C[CH2]C",
            [
                "CC",
                "CC",
            ],
            id="RemoveAtom on C[CH2]C",
        ),
        pytest.param(
            "c1ccc(cc1)[N+](=O)[O-]",
            [
                "C=CC(=CC)[N+](=O)[O-]",
                "C=C(C=CC)[N+](=O)[O-]",
                "C=CC=CC[N+](=O)[O-]",
                "CC=CC=C[N+](=O)[O-]",
                "C=CC=C(C)[N+](=O)[O-]",
            ],
            id="RemoveAtom on c1ccc(cc1)[N+](=O)[O-]",
        ),
        pytest.param(
            "C[C@@H](N)O",
            [
                "NCO",
                "CCO",
                "CCN",
            ],
            id="RemoveAtom on C[C@@H](N)O",
        ),
    ),
)
def test_actions_remove_atom(smiles: str, check_smiles: list[str]) -> None:
    """Test the RemoveAtom action."""
    check_actions_smiles(smiles, check_smiles, mg.RemoveAtomMG)


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
                "CCC",
                "C1=CC1",
                "C1#CC1",
                "CCC",
                "C1=CC1",
                "C1#CC1",
            ],
            id="ChangeBond on C1C[CH2]1",
        ),
        pytest.param(
            "CC[CH2]",
            [
                "[CH2]C=C",
                "C#C[CH2]",
            ],
            id="ChangeBond on CC[CH2]",
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
                "C=C(N)O",
                "O[C@H]1CN1",
                "O[C@H]1C=N1",
                "N[C@@H]1CO1",
                "CC(=N)O",
                "CC(N)=O",
                "C[C@H]1NO1",
            ],
            id="ChangeBond on C[C@@H](N)O",
        ),
    ),
)
def test_actions_change_bond(smiles: str, check_smiles: list[str]) -> None:
    """Test the ChangeBond action."""
    check_actions_smiles(smiles, check_smiles, mg.ChangeBondMG)


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
            "CN",
            [
                "NO",
                "NN",
                "NF",
                "CC",
                "CO",
                "CF",
            ],
            id="SubstituteAtom on CN",
        ),
        pytest.param(
            "C=N",
            [
                "N=O",
                "N=N",
                "C=C",
                "C=O",
            ],
            id="SubstituteAtom on C=N",
        ),
        pytest.param(
            "CC=[CH]",
            [
                "[CH]=CO",
                "[CH]=CN",
                "[CH]=CF",
                "[CH]=NC",
            ],
            id="SubstituteAtom on CC=[CH]",
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
            "C[C@@H](N)O",
            [
                "NC(O)O",
                "NC(N)O",
                "N[C@@H](O)F",
                "CN(N)O",
                "CC(C)O",
                "CC(O)O",
                "C[C@H](O)F",
                "CC(C)N",
                "CC(N)N",
                "C[C@@H](N)F",
            ],
            id="SubstituteAtom on C[C@@H](N)O",
        ),
    ),
)
def test_actions_substitute(smiles: str, check_smiles: list[str]) -> None:
    """Test the SubstituteAtom action."""
    check_actions_smiles(smiles, check_smiles, mg.SubstituteAtomMG)


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
            "C=N",
            [
                "CCN",
            ],
            id="InsertCarbon on C=N",
        ),
        pytest.param(
            "[N+]=[CH]",
            [
                "[CH]=C=[N+]",
            ],
            id="InsertCarbon on [N+]=[CH]",
        ),
        pytest.param(
            "C=NC",
            [
                "CCNC",
                "C=NCC",
            ],
            id="InsertCarbon on C=NC",
        ),
        pytest.param(
            "[CH]=[CH]",
            [
                "[CH]=C=[CH]",
            ],
            id="InsertCarbon on [CH]=[CH]",
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
        pytest.param(
            "C[C@@H](N)O",
            [
                "CC[C@@H](N)O",
                "C[C@@H](O)CN",
                "C[C@@H](N)CO",
            ],
            id="InsertCarbon on C[C@@H](N)O",
        ),
        pytest.param(
            "C=[S+]=C",
            [
                "C=[S+]=CC",
                "C=[S+]=CC",
            ],
            id="InsertCarbon on C=[S+]=C",
        ),
        pytest.param(
            "[S+]=[S+]",
            [
                "[S+]=C=[S+]",
            ],
            id="InsertCarbon on [S+]=[S+]",
        ),
    ),
)
def test_actions_insert_carbon(smiles: str, check_smiles: list[str]) -> None:
    """Test the InsertCarbon action."""
    check_actions_smiles(smiles, check_smiles, mg.InsertCarbonMG)


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
            "C=NC",
            [
                "CC",
            ],
            id="CutAtom on C=NC",
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
                "C=[C+]O",
            ],
            id="CutAtom on O[C+]=CC",
        ),
        pytest.param(
            "O[C+]C=C",
            [
                "C[C+]O",
            ],
            id="CutAtom on O[C+]C=C",
        ),
        pytest.param(
            "OC=C[C+]",
            [
                "[C+]CO",
                "[C+]CO",
            ],
            id="CutAtom on OC=C[C+]",
        ),
        pytest.param(
            "OC=C[C]",
            [
                "[C]CO",
                "[C]CO",
            ],
            id="CutAtom on OC=C[C]",
        ),
        pytest.param(
            "OC=C[CH3]",
            [
                "CCO",
                "CCO",
            ],
            id="CutAtom on OC=C[CH3]",
        ),
        pytest.param(
            "[N+]1=C=[N+]1",
            [
                "[N+]#[N+]",
            ],
            id="CutAtom on [N+]1=C=[N+]1",
        ),
        pytest.param(
            "[N+]1=C[N+]1",
            [],
            id="CutAtom on [N+]1=C[N+]1",
        ),
        pytest.param(
            "N1=C=C1",
            [
                "C#C",
                "C#N",
                "C#N",
            ],
            id="CutAtom on N1=C=C1",
        ),
        pytest.param(
            "N1C=[N+]1",
            [
                "C#[N+]",
                "N#[N+]",
            ],
            id="CutAtom on N1C=[N+]1",
        ),
        pytest.param(
            "[N+]1=CN1",
            [
                "N#[N+]",
                "C#[N+]",
            ],
            id="CutAtom on [N+]1=CN1",
        ),
        pytest.param(
            "[N+]1C=N1",
            [
                "[N+]=N",
                "C=[N+]",
            ],
            id="CutAtom on [N+]1C=N1",
        ),
        pytest.param(
            "N1C=N1",
            [
                "C#N",
                "N=N",
                "C=N",
            ],
            id="CutAtom on N1C=N1",
        ),
        pytest.param(
            "[N+]=C=[CH]",
            [
                "[CH]=[N+]",
            ],
            id="CutAtom on [N+]=C=[CH]",
        ),
        pytest.param(
            "C[C@@H](N)O",
            [],
            id="CutAtom on C[C@@H](N)O",
        ),
    ),
)
def test_actions_cut_atom(smiles: str, check_smiles: list[str]) -> None:
    """Test the CutAtom action."""
    check_actions_smiles(smiles, check_smiles, mg.CutAtomMG)


########################################
#               MoveGroup
########################################
@pytest.mark.parametrize(
    ("smiles", "check_smiles"),
    (
        pytest.param(
            "C1=CC=C[CH]=C1CC",
            [
                "CCC1=CC=CC=C1",
                "CCC1=CC=CC=C1",
                "CCC1=CC=CC=C1",
                "CCC1=CC=CC=C1",
                "CCC1=CC=CC=C1",
                "CCC1=CC=CC=C1",
                "CC1=CC=CC=C1C",
                "CC1=CC=CC(C)=C1",
                "CC1=CC=C(C)C=C1",
                "CC1=CC=CC(C)=C1",
                "CC1=C(C)C=CC=C1",
            ],
            id="MoveGroup on C1=CC=C[CH]=C1CC",
        ),
        pytest.param(
            "C1=CC=C[CH]=C1C1=CC=C[CH]=C1",
            [
                "C1=CC=C(C2=CC=CC=C2)C=C1",
                "C1=CC=C(C2=CC=CC=C2)C=C1",
                "C1=CC=C(C2=CC=CC=C2)C=C1",
                "C1=CC=C(C2=CC=CC=C2)C=C1",
                "C1=CC=C(C2=CC=CC=C2)C=C1",
                "C1=CC=C(C2=CC=CC=C2)C=C1",
                "C1=CC=C(C2=CC=CC=C2)C=C1",
                "C1=CC=C(C2=CC=CC=C2)C=C1",
                "C1=CC=C(C2=CC=CC=C2)C=C1",
                "C1=CC=C(C2=CC=CC=C2)C=C1",
            ],
            id="MoveGroup on C1=CC=C[CH]=C1C1=CC=C[CH]=C1",
        ),
        pytest.param(
            "c1ncncc1CCO",
            [
                "OCCC1=NC=NC=C1",
                "OCCC1=NC=CC=N1",
                "OCCC1=CC=NC=N1",
                "CC(O)C1=CN=CN=C1",
                "CCOC1=CN=CN=C1",
                "CC1=CN=CN=C1CO",
                "CC1=CN=C(CO)N=C1",
                "CC1=C(CO)N=CN=C1",
                "COCC1=CN=CN=C1",
                "CCC1=CN=CN=C1O",
                "CCC1=CN=C(O)N=C1",
                "CCC1=C(O)N=CN=C1",
                "CC(O)C1=CN=CN=C1",
            ],
            id="MoveGroup on c1ncncc1CCO",
        ),
        pytest.param(
            "C1NC1C1SC1",
            [
                "C1NC1C1CS1",
                "C1SC1N1CC1",
                "C1NC1[SH]1CC1",
                "C1NC1C1CS1",
            ],
            id="MoveGroup on C1NC1C1SC1",
        ),
        pytest.param(
            "C1[N+]=C1C1SC1",
            [
                "C1=[N+]C1C1CS1",
                "C1[N+]=C1[SH]1CC1",
                "C1[N+]=C1C1CS1",
            ],
            id="MoveGroup on C1[N+]=C1C1SC1",
        ),
        pytest.param(
            "C1C[N+]1=C1SC1",
            [
                "C1C[N+]1=S1CC1",
                "C1SC1=[N+]1CC1",
            ],
            id="MoveGroup on C1C[N+]1=C1SC1",
        ),
        pytest.param(
            "C1C[N+]1=[N+]1SC1",
            [],
            id="MoveGroup on C1C[N+]1=[N+]1SC1",
        ),
        pytest.param(
            "C=NC",
            [
                "CCN",
                "CC=N",
            ],
            id="MoveGroup on C=NC",
        ),
        pytest.param(
            "[N+]=CN",
            [
                "CN=[N+]",
            ],
            id="MoveGroup on [N+]=CN",
        ),
        pytest.param(
            "C[C@@H](N)O",
            [
                "CNCO",
                "COCN",
                "NCCO",
                "CCON",
                "NCCO",
                "CCNO",
            ],
            id="MoveGroup on C[C@@H](N)O",
        ),
        pytest.param(
            "C[C+](N)O",
            [],
            id="MoveGroup on C[C+](N)O",
        ),
        pytest.param(
            "S1C([CH2])=CC=C1",
            [
                "[CH2][SH]1C=CC=C1",
                "[CH2]C1=CSC=C1",
                "[CH2]C1=CSC=C1",
                "[CH2]C1=CC=CS1",
            ],
            id="MoveGroup on S1C([CH2])=CC=C1",
        ),
    ),
)
def test_actions_move_group(smiles: str, check_smiles: list[str]) -> None:
    """Test the MoveGroup action."""
    check_actions_smiles(smiles, check_smiles, mg.MoveGroupMG)


########################################
#                 AddGroup
########################################
@pytest.mark.parametrize(
    ("smiles", "check_smiles"),
    (
        pytest.param(
            "",
            [
                "C1=CSC=C1",
                "C1=CC=CC=C1",
                "O=[N+][O-]",
                "[N-]=[N+]=N",
                "O=[SH](=O)O",
            ],
            id="AddGroup on empty smiles",
        ),
        pytest.param(
            "[CH2]",
            [],
            id="AddGroup on [CH2]",
        ),
        pytest.param(
            "C",
            [
                "CC1=CC=CS1",
                "CC1=CC=CC=C1",
                "C[N+](=O)[O-]",
                "CN=[N+]=[N-]",
                "CS(=O)(=O)O",
            ],
            id="AddGroup on C",
        ),
        pytest.param(
            "N[CH2]",
            [
                "[CH2]NC1=CC=CS1",
                "[CH2]NC1=CC=CC=C1",
                "[CH2]N[N+](=O)[O-]",
                "[CH2]NN=[N+]=[N-]",
                "[CH2]NS(=O)(=O)O",
            ],
            id="AddGroup on N[CH2]",
        ),
    ),
)
def test_actions_add_group(smiles: str, check_smiles: list[str]) -> None:
    """Test the RemoveGroup action."""
    check_actions_smiles(smiles, check_smiles, mg.AddGroupMG)


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
            "C1=CC=CC=C1CC",
            [
                "C1=CC=CC=C1",
                "CC1=CC=CC=C1",
            ],
            id="RemoveGroup on C1=CC=CC=C1CC",
        ),
        pytest.param(
            "c1ccc[S+]1[N+]1[CH]ccocc1",
            [],
            id="RemoveGroup on c1ccc[S+]1[N+]1[CH]ccocc1",
        ),
        pytest.param(
            "C=C",
            [
                "C",
                "C",
            ],
            id="RemoveGroup on C=C",
        ),
        pytest.param(
            "C1=CC=CC=C1[CH]C",
            [],
            id="RemoveGroup on C1=CC=CC=C1[CH]C",
        ),
        pytest.param(
            "C1=CC=CC=C1C[CH]",
            [
                "C1=CC=CC=C1",
            ],
            id="RemoveGroup on C1=CC=CC=C1C[CH]",
        ),
    ),
)
def test_actions_remove_group(smiles: str, check_smiles: list[str]) -> None:
    """Test the RemoveGroup action."""
    check_actions_smiles(smiles, check_smiles, mg.RemoveGroupMG)


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
