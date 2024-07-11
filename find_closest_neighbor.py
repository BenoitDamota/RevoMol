import time
from dataclasses import dataclass

from evomol import evaluation as evaluator
from evomol.action import molecular_graph as mg
from evomol.representation import SMILES, MolecularGraph, Molecule


def find_neighbors(molecule: Molecule, depth: int) -> list[Molecule]:
    """Recursively find the neighbors of a molecule up to a certain depth

    Args:
        molecule (Molecule): Molecule to explore
        depth (int): Maximum depth to explore

    Returns:
        list[Molecule]: List of molecules found
    """
    # list all possible actions in the molecule representation
    molecule.compute_possible_actions()

    new_mols = []

    # while there are possible actions
    while molecule.nb_remaining_actions() != 0:
        # get the first action space
        action_space = list(molecule.possible_actions.keys())[0]
        # while there are possible actions in the action space
        while (
            molecule.possible_actions is not None
            and molecule.possible_actions.get(action_space) is not None
        ):
            # get the first action
            action = list(molecule.possible_actions[action_space].keys())[0]
            # apply the action to get a new molecule
            new_mol = molecule.possible_actions[action_space][action][0].apply()

            # add the new molecule to the list
            new_mols.append(new_mol)

            # recursively find the neighbors of the new molecule
            if depth > 1:
                new_mols.extend(find_neighbors(new_mol, depth - 1))

    return new_mols


def find_closest_neighbors(start_smiles: str, eval_name: str, nb_heavy_atoms: int):
    """Enumerate the neighbors of a SMILES

    Args:
        start_smiles (str): SMILES to start from
        eval_name (str): chembl or chembl_zinc as filter
        nb_heavy_atoms (int): maximum number of heavy atoms
        depth (int): maximum depth to explore
    """
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
            evaluator.UnknownGCF(path_db="external_data/gcf1.txt", name="chembl"),
            evaluator.FilterUnknownGCF(threshold=0, name="chembl"),
            evaluator.UnknownECFP(
                path_db="external_data/ecfp4_ChEMBL.txt", radius=2, name="chembl"
            ),
            evaluator.FilterUnknownECFP(threshold=0, name="chembl"),
        ]
    elif eval_name == "chembl_zinc":
        evaluations = [
            evaluator.UnknownGCF(path_db="external_data/gcf2.txt", name="chembl_zinc"),
            evaluator.FilterUnknownGCF(threshold=0, name="chembl_zinc"),
            evaluator.UnknownECFP(
                path_db="external_data/ecfp4_ChEMBL_ZINC.txt",
                radius=2,
                name="chembl_zinc",
            ),
            evaluator.FilterUnknownECFP(threshold=0, name="chembl_zinc"),
        ]

    valid_mols: dict[str, int] = {}
    depth = 0
    nb_mols = 0

    print(start_smiles, f"- filter: {eval_name}")

    while not valid_mols:
        depth += 1

        time_start = time.time()
        mols = find_neighbors(Molecule(start_smiles), depth)

        nb_mols = len(mols)

        for mol in mols:
            can_smi = mol.get_representation(MolecularGraph).canonical_smiles
            if can_smi in valid_mols:
                valid_mols[can_smi] += 1
                continue
            valid = True
            for eval_ in evaluations:
                try:
                    eval_.evaluate(mol)
                except evaluator.EvaluationError:
                    valid = False
                    break
            if valid:
                valid_mols[can_smi] = 1
        duration = time.time() - time_start
        if not valid_mols:
            print(
                f"\tno valid molecules at depth {depth} "
                f"(among {nb_mols} molecules)"
                f" in {duration:.2f}s"
            )
        else:
            print(
                f"\tvalid molecules at depth {depth} "
                f"({len(valid_mols)}/{nb_mols} molecules)"
                f" in {duration:.2f}s"
            )

    for mol, nb_found in valid_mols.items():
        print("\t\t", mol, nb_found)


@dataclass
class Counter:
    nb_found: int
    min_depth: int


def find_unique_neighbors(
    molecule: Molecule, max_depth: int
) -> dict[Molecule, Counter]:
    found_molecules: dict[Molecule, Counter] = {molecule: Counter(1, 0)}

    to_explore: set[Molecule] = {molecule}

    while to_explore:
        current_mol = to_explore.pop()
        current_depth = found_molecules[current_mol].min_depth
        if current_depth >= max_depth:
            continue
        current_mol.compute_possible_actions()
        while current_mol.nb_remaining_actions() != 0:
            action_space = list(current_mol.possible_actions.keys())[0]
            while (
                current_mol.possible_actions is not None
                and current_mol.possible_actions.get(action_space) is not None
            ):
                action = list(current_mol.possible_actions[action_space].keys())[0]
                new_mol = current_mol.possible_actions[action_space][action][0].apply()
                if new_mol in found_molecules:
                    found_molecules[new_mol].nb_found += 1
                    found_molecules[new_mol].min_depth = min(
                        found_molecules[new_mol].min_depth, current_depth + 1
                    )
                else:
                    found_molecules[new_mol] = Counter(1, current_depth + 1)
                    to_explore.add(new_mol)

    return found_molecules


def find_closest_neighbors_unique(
    start_smiles: str, eval_name: str, nb_heavy_atoms: int
):
    """Enumerate the neighbors of a SMILES

    Args:
        start_smiles (str): SMILES to start from
        eval_name (str): chembl or chembl_zinc as filter
        nb_heavy_atoms (int): maximum number of heavy atoms
        depth (int): maximum depth to explore
    """
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
            evaluator.UnknownGCF(path_db="external_data/gcf1.txt", name="chembl"),
            evaluator.FilterUnknownGCF(threshold=0, name="chembl"),
            evaluator.UnknownECFP(
                path_db="external_data/ecfp4_ChEMBL.txt", radius=2, name="chembl"
            ),
            evaluator.FilterUnknownECFP(threshold=0, name="chembl"),
        ]
    elif eval_name == "chembl_zinc":
        evaluations = [
            evaluator.UnknownGCF(path_db="external_data/gcf2.txt", name="chembl_zinc"),
            evaluator.FilterUnknownGCF(threshold=0, name="chembl_zinc"),
            evaluator.UnknownECFP(
                path_db="external_data/ecfp4_ChEMBL_ZINC.txt",
                radius=2,
                name="chembl_zinc",
            ),
            evaluator.FilterUnknownECFP(threshold=0, name="chembl_zinc"),
        ]

    valid_mols: dict[str, int] = {}
    depth = 0
    nb_mols = 0

    can_smi_start = (
        Molecule(start_smiles).get_representation(MolecularGraph).canonical_smiles
    )

    print(start_smiles, f" can : {can_smi_start} - filter: {eval_name}")

    while not valid_mols:
        depth += 1

        time_start = time.time()

        mols: dict[Molecule, Counter] = find_unique_neighbors(
            Molecule(can_smi_start), depth
        )

        nb_mols = len(mols) - 1

        for mol in mols:
            can_smi = mol.get_representation(MolecularGraph).canonical_smiles
            if can_smi == can_smi_start:
                continue
            valid = True
            for eval_ in evaluations:
                try:
                    eval_.evaluate(mol)
                except evaluator.EvaluationError:
                    valid = False
                    break
            if valid:
                valid_mols[can_smi] = mols[mol].nb_found
        duration = time.time() - time_start
        print(
            "\tno" if not valid_mols else "\t",
            f"valid molecules at depth {depth} -",
            f"{nb_mols} unique / {sum([mols[mol].nb_found for mol in mols]) -1} molecules"
            f" in {duration:.2f}s",
        )

    for mol, nb_found in valid_mols.items():
        print("\t\t", mol, nb_found)


if __name__ == "__main__":

    # eval_name = "chembl"
    # eval_name = "chembl_zinc"

    smiles = [
        # "C1=CSC(=C2SC=CS2)S1",  # TTF
        # "N1=S=NC2=C1N=S=N2",  # DD
        # "CC1=NOC(C)(O)C1=NO",
        # "CN1ON(C)ON(C)O1",
        # "c1coc2occoc=2o1",
        "CSC1N=NC(SC)N=N1",
        # "S=c1[nH]ssc2nnc1=2",
        # "N1=NC(=C2N=NN=N2)N=N1",
        # "N1=S=NC2=C1N=S=N2",
    ]

    for smi in smiles:
        for eval_name in [
            # "chembl",
            "chembl_zinc",
        ]:
            # find_closest_neighbors(smi, eval_name, 10)
            find_closest_neighbors_unique(smi, eval_name, 10)
        print()

    # examples from CLI:
    # python find_closest_neighbor.py C chembl 10

    # TTF
    # python find_closest_neighbor.py "C1=CSC(=C2SC=CS2)S1" chembl 10

    # DD
    # python find_closest_neighbor.py "N1=S=NC2=C1N=S=N2" chembl 10

    # typer.run(find_closest_neighbors)
