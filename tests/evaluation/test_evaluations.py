"""
Tests for the molecular graph actions.
"""

import pytest

import evomol.evaluation as evaluator
from evomol import default_parameters as dp
from evomol.evaluation.generic_cyclic_features import (
    convert_to_carbon,
    divide_subgraphs,
    list_gcf,
)
from evomol.evaluation.unknown_ecfp import list_ecfp
from evomol.representation import Molecule


@pytest.fixture(autouse=True, scope="module")
def setup_parameters() -> None:
    """set the parameters for the tests."""
    dp.setup_default_parameters(
        accepted_atoms=["C", "O", "N", "F"],
        max_heavy_atoms=38,
    )
    dp.setup_default_action_space()


########################################
#              CycleScore
########################################
@pytest.mark.parametrize(
    ("smiles", "value"),
    (
        pytest.param(
            "",
            0,
            id="CycleScore on empty smiles",
        ),
        pytest.param(
            "C1CCCCC1",
            0,
            id="CycleScore on C1CCCCC1",
        ),
        pytest.param(
            "c1ccccccc1",
            -2,
            id="CycleScore on c1ccccccc1",
        ),
        pytest.param(
            "c1ccccccccc1",
            -4,
            id="CycleScore on c1ccccccccc1",
        ),
    ),
)
def test_cycle_score(smiles: str, value: float) -> None:
    """Test the Cycle Score."""
    assert evaluator.CycleScore.evaluate(Molecule(smiles)) == value


########################################
#                  GCF
########################################
@pytest.mark.parametrize(
    ("smiles", "features"),
    (
        pytest.param(
            "",
            set(),
            id="GCF : list_gcf on empty smiles",
        ),
        pytest.param(
            "C",
            set(),
            id="GCF : list_gcf on C",
        ),
        pytest.param(
            "C1=CC=CC=C1",
            {"C1=CC=CC=C1"},
            id="GCF : list_gcf on C1=CC=CC=C1",
        ),
        pytest.param(
            "C1=CC2CCCCC2C=C1",
            {"C1=CC2CCCCC2C=C1"},
            id="GCF : list_gcf on C1=CC2CCCCC2C=C1",
        ),
        pytest.param(
            "C1=CC=CC=C1C1=CC=CC=C1",
            {"C1=CC=CC=C1"},
            id="GCF : list_gcf on C1=CC=CC=C1C1=CC=CC=C1",
        ),
        pytest.param(
            "C1=CN=CN=C1C1=CN=CN=C1",
            {"C1=CC=CC=C1"},
            id="GCF : list_gcf on C1=CN=CN=C1C1=CN=CN=C1",
        ),
        pytest.param(
            "S1=S=S#SS=S1S1=S=S#S#S=S1",
            {"C1=CC=CC#[SH]=1", "C1=CC=[SH]#S#[SH]=1"},
            id="GCF : list_gcf on S1=S=S#SS=S1S1=S=S#S#S=S1",
        ),
    ),
)
def test_gcf_conversion(smiles: str, features: set[str]) -> None:
    """Test the list_gcf function."""
    assert list_gcf(smiles) == features


@pytest.mark.parametrize(
    ("smiles", "subgraphs"),
    (
        pytest.param(
            "",
            {""},
            id="GCF : divide_subgraph on empty smiles",
        ),
        pytest.param(
            "C",
            {"C"},
            id="GCF : divide_subgraph on C",
        ),
        pytest.param(
            "C1=CC=CC=C1",
            {"C1=CC=CC=C1"},
            id="GCF : divide_subgraph on C1=CC=CC=C1",
        ),
        pytest.param(
            "C1=CC2CCCCC2C=C1",
            {"C1=CC2CCCCC2C=C1"},
            id="GCF : divide_subgraph on C1=CC2CCCCC2C=C1",
        ),
        pytest.param(
            "C1=CC=CC=C1C1=CC=CC=C1",
            {"C1=CC=CC=C1"},
            id="GCF : divide_subgraph on C1=CC=CC=C1C1=CC=CC=C1",
        ),
        pytest.param(
            "C1=CN=CN=C1C1=CN=CN=C1",
            {"C1=CN=CN=C1"},
            id="GCF : divide_subgraph on C1=CN=CN=C1C1=CN=CN=C1",
        ),
    ),
)
def test_divide_subgraphs(smiles: str, subgraphs: set[str]) -> None:
    """Test the divide_subgraphs function."""
    assert divide_subgraphs(smiles) == subgraphs


@pytest.mark.parametrize(
    ("smiles", "carbon_smiles"),
    (
        pytest.param(
            "",
            "",
            id="GCF : convert_to_carbon on empty smiles",
        ),
        pytest.param(
            "C",
            "",
            id="GCF : convert_to_carbon on C",
        ),
        pytest.param(
            "C1=CC=CC=C1",
            "C1=CC=CC=C1",
            id="GCF : convert_to_carbon on C1=CC=CC=C1",
        ),
        pytest.param(
            "C1=CC2CCCCC2C=C1",
            "C1=CC2CCCCC2C=C1",
            id="GCF : convert_to_carbon on C1=CC2CCCCC2C=C1",
        ),
        pytest.param(
            "C1=CC=CC=C1C1=CC=CC=C1",
            "C1=CC=C(C2=CC=CC=C2)C=C1",
            id="GCF : convert_to_carbon on C1=CC=CC=C1C1=CC=CC=C1",
        ),
        pytest.param(
            "C1=CN=CN=C1C1=CN=CN=C1",
            "C1=CC=C(C2=CC=CC=C2)C=C1",
            id="GCF : convert_to_carbon on C1=CN=CN=C1C1=CN=CN=C1",
        ),
        pytest.param(
            "C=S1=C=CC=CC=1",
            "C=S1=C=CC=CC=1",
            id="GCF : convert_to_carbon on C=S1=C=CC=CC=1",
        ),
        pytest.param(
            "S1=SS=SS=S1",
            "C1=CC=CC=C1",
            id="GCF : convert_to_carbon on S1=SS=SS=S1",
        ),
        pytest.param(
            "S1=S=S#SS=S1",
            "C1=CC=CC#[SH]=1",
            id="GCF : convert_to_carbon on S1=S=S#SS=S1",
        ),
        pytest.param(
            "S1=S=S#S#S=S1",
            "C1=CC=[SH]#S#[SH]=1",
            id="GCF : convert_to_carbon on S1=S=S#S#S=S1",
        ),
    ),
)
def test_convert_to_carbon(smiles: str, carbon_smiles: str) -> None:
    """Test the convert_to_carbon function."""
    assert convert_to_carbon(smiles) == carbon_smiles


########################################
#                  ECFP
########################################
@pytest.mark.parametrize(
    ("smiles", "ecfp"),
    (
        pytest.param(
            "",
            [],
            id="ECFP : list_ecfp on empty smiles",
        ),
        pytest.param(
            "C",
            [2246733040],
            id="ECFP : list_ecfp on C",
        ),
        pytest.param(
            "C1=CC=CC=C1",
            [98513984, 2763854213, 3218693969],
            id="ECFP : list_ecfp on C1=CC=CC=C1",
        ),
        pytest.param(
            "C1=CC2CCCCC2C=C1",
            [
                98466654,
                1721979728,
                1959335747,
                2117068077,
                2142032900,
                2308821956,
                2447748155,
                2832976762,
                2968968094,
                2976033787,
                3218693969,
                3570805596,
                4247667222,
            ],
            id="ECFP : list_ecfp on C1=CC2CCCCC2C=C1",
        ),
        pytest.param(
            "C1=CC=CC=C1C1=CC=CC=C1",
            [
                98513984,
                669222828,
                951226070,
                2353112200,
                2763854213,
                2984966880,
                3217380708,
                3218693969,
                3999906991,
            ],
            id="ECFP : list_ecfp on C1=CC=CC=C1C1=CC=CC=C1",
        ),
        pytest.param(
            "C1=CN=CN=C1C1=CN=CN=C1",
            [
                64070612,
                725322217,
                951226070,
                1100037548,
                1103832125,
                1178619885,
                1542020373,
                1717044408,
                1755931940,
                2041434490,
                2065735937,
                2173788520,
                2353112200,
                2519390547,
                3118255683,
                3217380708,
                3218693969,
                3229069614,
                3229334257,
                3776905034,
                3777168895,
            ],
            id="ECFP : list_ecfp on C1=CN=CN=C1C1=CN=CN=C1",
        ),
        pytest.param(
            "S1=S=S#SS=S1S1=S=S#S#S=S1",
            [
                487731247,
                803858804,
                830972580,
                867306036,
                872258609,
                1101358677,
                1471721975,
                1876084888,
                2344124339,
                2728471623,
                2781990650,
                2812620228,
                2849135538,
                2850007887,
                2856568483,
                2980543007,
                3193912799,
                3201391669,
                3201416467,
                3201428904,
                3254496537,
                3351556771,
                3818747402,
                3834138505,
                3990869196,
                4181607724,
            ],
            id="ECFP : list_ecfp on S1=S=S#SS=S1S1=S=S#S#S=S1",
        ),
    ),
)
def test_list_ecfp(smiles: str, ecfp: list[int]) -> None:
    """Test the list_ecfp function."""
    assert list_ecfp(Molecule(smiles)) == ecfp


########################################
#                  logP
########################################
@pytest.mark.parametrize(
    ("smiles", "value"),
    (
        pytest.param(
            "",
            0,
            id="logP on empty smiles",
        ),
        pytest.param(
            "C1CCCCC1",
            2.3406,
            id="logP on C1CCCCC1",
        ),
        pytest.param(
            "c1ccccccc1",
            2.2248,
            id="logP on c1ccccccc1",
        ),
        pytest.param(
            "c1ccccccccc1",
            2.7810,
            id="logP on c1ccccccccc1",
        ),
        pytest.param(
            "C1=CC=CC=C1C1=CC=CC=C1",
            3.0912,
            id="logP on C1=CC=CC=C1C1=CC=CC=C1",
        ),
        pytest.param(
            "C1=CN=CN=C1C1=CN=CN=C1",
            0.9796,
            id="logP on C1=CN=CN=C1C1=CN=CN=C1",
        ),
        pytest.param(
            "S1=S=S#SS=S1S1=S=S#S#S=S1",
            -0.3337999999999999,
            id="logP on S1=S=S#SS=S1S1=S=S#S#S=S1",
        ),
    ),
)
def test_logp(smiles: str, value: float) -> None:
    """Test the logP evaluation."""
    assert evaluator.LogP.evaluate(Molecule(smiles)) == pytest.approx(value)


########################################
#                  QED
########################################
@pytest.mark.parametrize(
    ("smiles", "value"),
    (
        pytest.param(
            "",
            0.33942358984550913,
            id="QED on empty smiles",
        ),
        pytest.param(
            "C1CCCCC1",
            0.42231618686094674,
            id="QED on C1CCCCC1",
        ),
        pytest.param(
            "c1ccccccc1",
            0.4417581656415751,
            id="QED on c1ccccccc1",
        ),
        pytest.param(
            "c1ccccccccc1",
            0.5062560331071598,
            id="QED on c1ccccccccc1",
        ),
        pytest.param(
            "C1=CC=CC=C1C1=CC=CC=C1",
            0.5905018880213778,
            id="QED on C1=CC=CC=C1C1=CC=CC=C1",
        ),
        pytest.param(
            "C1=CN=CN=C1C1=CN=CN=C1",
            0.6191775203459348,
            id="QED on C1=CN=CN=C1C1=CN=CN=C1",
        ),
        pytest.param(
            "S1=S=S#SS=S1S1=S=S#S#S=S1",
            0.2562939068049502,
            id="QED on S1=S=S#SS=S1S1=S=S#S#S=S1",
        ),
    ),
)
def test_qed(smiles: str, value: float) -> None:
    """Test the QED evaluation."""
    assert evaluator.QED.evaluate(Molecule(smiles)) == pytest.approx(value)


########################################
#               SA Score
########################################
@pytest.mark.parametrize(
    ("smiles", "value"),
    (
        pytest.param(
            "",
            1,
            id="SA Score on empty smiles",
        ),
        pytest.param(
            "C1CCCCC1",
            1.0,
            id="SA Score on C1CCCCC1",
        ),
        pytest.param(
            "c1ccccccc1",
            2.3639650013068003,
            id="SA Score on c1ccccccc1",
        ),
        pytest.param(
            "c1ccccccccc1",
            2.670852803830762,
            id="SA Score on c1ccccccccc1",
        ),
        pytest.param(
            "C1=CC=CC=C1C1=CC=CC=C1",
            3.467864638330914,
            id="SA Score on C1=CC=CC=C1C1=CC=CC=C1",
        ),
        pytest.param(
            "C1=CN=CN=C1C1=CN=CN=C1",
            4.7614646383309145,
            id="SA Score on C1=CN=CN=C1C1=CN=CN=C1",
        ),
        pytest.param(
            "S1=S=S#SS=S1S1=S=S#S#S=S1",
            8.232101638846137,
            id="SA Score on S1=S=S#SS=S1S1=S=S#S#S=S1",
        ),
    ),
)
def test_sa_score(smiles: str, value: float) -> None:
    """Test the SA Score evaluation."""
    assert evaluator.SAScore.evaluate(Molecule(smiles)) == pytest.approx(value)
