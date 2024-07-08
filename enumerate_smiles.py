import time
import typer
from evomol import evaluation as evaluator
from evomol.action import molecular_graph as mg
from evomol.representation import SMILES, MolecularGraph, Molecule


def explore_smiles(
    smiles: str, depth: int, evaluations: list[evaluator.Evaluation]
) -> dict[str, int]:
    to_try = set([smiles])
    found = {smiles: 1}
    while to_try:
        current_smiles = to_try.pop()
        mols = find_neighbors(Molecule(current_smiles), depth)
        for mol in mols:
            smiles = mol.get_representation(MolecularGraph).canonical_smiles
            if smiles in to_try:
                continue
            if smiles in found:
                found[smiles] += 1
                continue

            valid = True
            for eval_ in evaluations:
                try:
                    eval_.evaluate(mol)
                except evaluator.EvaluationError:
                    valid = False
                    break

            if valid:
                to_try.add(smiles)
                if smiles not in found:
                    found[smiles] = 1
    return found


def find_neighbors(molecule: Molecule, depth: int) -> list[Molecule]:
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

            new_mols.append(new_mol)
            if depth > 1:
                new_mols.extend(find_neighbors(new_mol, depth - 1))

    return new_mols


def from_c_change_heavy_atom(
    start_smiles: str, eval_name: str, nb_heavy_atoms: int, depth: int
):
    Molecule.id_representation_class = SMILES
    Molecule.representations_class = [MolecularGraph]
    Molecule.max_heavy_atoms = nb_heavy_atoms

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

    evaluations = []
    if eval_name == "chembl":
        evaluations = [
            evaluator.UnknownGCF(path_db="external_data/gcf1.txt"),
            evaluator.FilterUnknownGCF(threshold=0),
            evaluator.UnknownECFP(path_db="external_data/ecfp4_ChEMBL.txt", radius=2),
            evaluator.FilterUnknownECFP(threshold=0),
        ]
    elif eval_name == "chembl_zinc":
        evaluations = [
            evaluator.UnknownGCF(path_db="external_data/gcf2.txt"),
            evaluator.FilterUnknownGCF(threshold=0),
            evaluator.UnknownECFP(
                path_db="external_data/ecfp4_ChEMBL_ZINC.txt", radius=2
            ),
            evaluator.FilterUnknownECFP(threshold=0),
        ]

    print("molecule,heavy_atom_limit,depth,eval,nb_molecules,time")

    start = time.time()

    found = explore_smiles(start_smiles, depth, evaluations)

    duration = time.time() - start
    print(
        ",".join(
            map(
                str,
                [
                    start_smiles,
                    nb_heavy_atoms,
                    depth,
                    eval_name,
                    len(found),
                    f"{duration:.2f}",
                ],
            )
        )
    )

    file = f"output/enumerate_{start_smiles}_{eval_name}_{nb_heavy_atoms}.txt"

    with open(file, "w") as f:
        for mol, nb in found.items():
            if mol == "":
                f.write(f'"" {nb}\n')
            else:
                f.write(f"{mol} {nb}\n")


if __name__ == "__main__":
    # parallel -j 20 python enumerate_smiles.py C {1} {2} 2 ::: chembl chembl_zinc ::: {1..10}
    typer.run(from_c_change_heavy_atom)
