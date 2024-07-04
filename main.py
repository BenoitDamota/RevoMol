import time

import typer

from evomol import evaluation as evaluator
from evomol.action import molecular_graph as mg
from evomol.representation import SMILES, MolecularGraph, Molecule

# from evomol.evomol import main


def explore_smiles(
    smiles: str, depth: int, evaluations: evaluator.Evaluation, apply_evaluation: bool
) -> dict[str, int]:
    to_try = set(smiles)
    tried = set()
    found = {}
    while to_try:
        smiles = to_try.pop()
        tried.add(smiles)
        mols = find_neighbors(
            molecule=Molecule(smiles),
            depth=depth,
            evaluations=evaluations,
            apply_evaluation=apply_evaluation,
        )
        for mol in mols:
            valid = True
            for eval_ in evaluations:
                try:
                    eval_.evaluate(mol)
                except evaluator.EvaluationError:
                    valid = False
                    break
            if valid:
                smiles = mol.get_representation(MolecularGraph).canonical_smiles
                if smiles not in tried:
                    to_try.add(smiles)
                if smiles not in found:
                    found[smiles] = 1
                else:
                    found[smiles] += 1
    return found


def find_neighbors(
    molecule: Molecule, depth: int, evaluations, apply_evaluation: bool
) -> list[Molecule]:
    molecule.list_all_possible_actions()
    new_mols = []
    while molecule.nb_possible_actions() != 0:
        action_space = list(molecule.possible_actions.keys())[0]
        while (
            molecule.possible_actions is not None
            and molecule.possible_actions.get(action_space) is not None
        ):
            action = list(molecule.possible_actions[action_space].keys())[0]
            new_mol = molecule.possible_actions[action_space][action][0].apply()

            if apply_evaluation:
                valid = True
                for eval_ in evaluations:
                    try:
                        eval_.evaluate(new_mol)
                    except evaluator.EvaluationError:
                        valid = False
                        break
                if not valid:
                    continue

            new_mols.append(new_mol)
            if depth > 1:
                new_mols.extend(
                    find_neighbors(new_mol, depth - 1, evaluations, apply_evaluation)
                )

    return new_mols


def from_c_change_heavy_atom(eval_name: str, nb_heavy_atoms: int):
    Molecule.id_representation_class = SMILES
    Molecule.representations_class = [MolecularGraph]
    Molecule.max_heavy_atoms = 10

    Molecule.accepted_atoms = ["C", "O", "N", "F", "S"]

    MolecularGraph.action_space = [
        mg.AddAtomMG,
        # mg.AddGroupMG,
        mg.ChangeBondMG,
        mg.CutAtomMG,
        mg.InsertCarbonMG,
        mg.MoveGroupMG,
        mg.RemoveAtomMG,
        # mg.RemoveGroupMG,
        mg.SubstituteAtomMG,
    ]

    evaluations = {
        "chembl": [
            evaluator.UnknownGCF(path_db="external_data/gcf1.txt"),
            evaluator.FilterUnknownGCF(threshold=0),
            evaluator.UnknownECFP(path_db="external_data/ecfp4_ChEMBL.txt", radius=2),
            evaluator.FilterUnknownECFP(threshold=0),
        ],
        "chembl_zinc": [
            evaluator.UnknownGCF(path_db="external_data/gcf2.txt"),
            evaluator.FilterUnknownGCF(threshold=0),
            evaluator.UnknownECFP(
                path_db="external_data/ecfp4_ChEMBL_ZINC.txt", radius=2
            ),
            evaluator.FilterUnknownECFP(threshold=0),
        ],
    }

    smiles = "C"
    apply_evaluation = False
    depth = 2

    print("molecule,heavy_atom_limit,depth,eval,nb_molecules,time")

    # for nb_heavy_atoms in range(1, 11):
    Molecule.max_heavy_atoms = nb_heavy_atoms
    # for eval_name, evaluations in evaluationss.items():
    # for apply_evaluation in [True, False]:

    start = time.time()

    found = explore_smiles(
        smiles=smiles,
        depth=depth,
        evaluations=evaluations[eval_name],
        apply_evaluation=apply_evaluation,
    )

    duration = time.time() - start
    print(
        ",".join(
            map(
                str,
                [
                    "C",
                    apply_evaluation,
                    nb_heavy_atoms,
                    2,
                    eval_name,
                    len(found),
                    # len(mols),
                    # len(set(mols)),
                    # len(valid_mols),
                    # len(set(valid_mols)),
                    f"{duration:.2f}",
                    # f"{duration_post_filter:.2f}",
                ],
            )
        )
    )
    # print(set(valid_mols))
    with open(f"output/{eval_name}_{nb_heavy_atoms}.txt", "w") as f:
        for mol, nb in found.items():
            if mol == "":
                f.write(f'"" {nb}\n')
            else:
                f.write(f"{mol} {nb}\n")


if __name__ == "__main__":
    typer.run(from_c_change_heavy_atom)
    # main()
