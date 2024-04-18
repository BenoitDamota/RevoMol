from collections import namedtuple
from itertools import product as product
from rdkit import Chem as Chem
from typing import Any, Optional

class Default:
    timeout: Any = ...
    timeoutString: str = ...
    maximize: str = ...
    atomCompare: str = ...
    bondCompare: str = ...
    matchValences: bool = ...
    ringMatchesRingOnly: bool = ...
    completeRingsOnly: bool = ...

class AtomSmartsNoAromaticity(dict):
    def __missing__(self, eleno: Any): ...

def atom_typer_any(atoms: Any): ...
def atom_typer_elements(atoms: Any): ...
def atom_typer_isotopes(atoms: Any): ...
def bond_typer_any(bonds: Any): ...
def bond_typer_bondtypes(bonds: Any): ...

atom_typers: Any
bond_typers: Any
default_atom_typer: Any
default_bond_typer: Any

def get_isotopes(mol: Any): ...
def set_isotopes(mol: Any, isotopes: Any) -> None: ...
def save_isotopes(mol: Any, isotopes: Any) -> None: ...
def save_atom_classes(mol: Any, atom_classes: Any) -> None: ...
def get_selected_atom_classes(mol: Any, atom_indices: Any): ...
def restore_isotopes(mol: Any) -> None: ...
def assign_isotopes_from_class_tag(mol: Any, atom_class_tag: Any) -> None: ...

class TypedMolecule:
    rdmol: Any = ...
    rdmol_atoms: Any = ...
    rdmol_bonds: Any = ...
    atom_smarts_types: Any = ...
    bond_smarts_types: Any = ...
    canonical_bondtypes: Any = ...
    def __init__(
        self,
        rdmol: Any,
        rdmol_atoms: Any,
        rdmol_bonds: Any,
        atom_smarts_types: Any,
        bond_smarts_types: Any,
        canonical_bondtypes: Any,
    ) -> None: ...

class FragmentedTypedMolecule:
    rdmol: Any = ...
    rdmol_atoms: Any = ...
    orig_atoms: Any = ...
    orig_bonds: Any = ...
    atom_smarts_types: Any = ...
    bond_smarts_types: Any = ...
    canonical_bondtypes: Any = ...
    def __init__(
        self,
        rdmol: Any,
        rdmol_atoms: Any,
        orig_atoms: Any,
        orig_bonds: Any,
        atom_smarts_types: Any,
        bond_smarts_types: Any,
        canonical_bondtypes: Any,
    ) -> None: ...

class TypedFragment:
    rdmol: Any = ...
    orig_atoms: Any = ...
    orig_bonds: Any = ...
    atom_smarts_types: Any = ...
    bond_smarts_types: Any = ...
    canonical_bondtypes: Any = ...
    def __init__(
        self,
        rdmol: Any,
        orig_atoms: Any,
        orig_bonds: Any,
        atom_smarts_types: Any,
        bond_smarts_types: Any,
        canonical_bondtypes: Any,
    ) -> None: ...

def get_canonical_bondtypes(
    rdmol: Any, bonds: Any, atom_smarts_types: Any, bond_smarts_types: Any
): ...
def get_typed_molecule(
    rdmol: Any,
    atom_typer: Any,
    bond_typer: Any,
    matchValences: Any = ...,
    ringMatchesRingOnly: Any = ...,
): ...
def get_specified_types(rdmol: Any, atom_types: Any, ringMatchesRingOnly: Any): ...
def convert_input_to_typed_molecules(
    mols: Any,
    atom_typer: Any,
    bond_typer: Any,
    matchValences: Any,
    ringMatchesRingOnly: Any,
): ...
def get_counts(it: Any): ...
def intersect_counts(counts1: Any, counts2: Any): ...
def get_canonical_bondtype_counts(typed_mols: Any): ...
def remove_unknown_bondtypes(typed_mol: Any, supported_canonical_bondtypes: Any): ...
def find_upper_fragment_size_limits(rdmol: Any, atoms: Any): ...

EnumerationMolecule = namedtuple("Molecule", "rdmol atoms bonds directed_edges")

Atom = namedtuple("Atom", "real_atom atom_smarts bond_indices is_in_ring")

Bond = namedtuple(
    "Bond", "real_bond bond_smarts canonical_bondtype atom_indices is_in_ring"
)

DirectedEdge = namedtuple("DirectedEdge", "bond_index end_atom_index")

Subgraph = namedtuple("Subgraph", "atom_indices bond_indices")

def get_typed_fragment(typed_mol: Any, atom_indices: Any): ...
def fragmented_mol_to_enumeration_mols(typed_mol: Any, minNumAtoms: int = ...): ...
def gen_primes() -> None: ...

class CangenNode:
    index: Any = ...
    atom_smarts: Any = ...
    value: int = ...
    neighbors: Any = ...
    rank: int = ...
    outgoing_edges: Any = ...
    def __init__(self, index: Any, atom_smarts: Any) -> None: ...

OutgoingEdge = namedtuple(
    "OutgoingEdge", "from_atom_index bond_index bond_smarts other_node_idx other_node"
)

def get_initial_cangen_nodes(
    subgraph: Any,
    enumeration_mol: Any,
    atom_assignment: Any,
    do_initial_assignment: bool = ...,
): ...
def rerank(cangen_nodes: Any) -> None: ...
def find_duplicates(cangen_nodes: Any, start: Any, end: Any): ...
def canon(cangen_nodes: Any): ...
def get_closure_label(bond_smarts: Any, closure: Any): ...
def generate_smarts(cangen_nodes: Any): ...
def make_canonical_smarts(
    subgraph: Any, enumeration_mol: Any, atom_assignment: Any
): ...
def make_arbitrary_smarts(
    subgraph: Any, enumeration_mol: Any, atom_assignment: Any
): ...
def powerset(iterable: Any): ...
def nonempty_powerset(iterable: Any): ...

tiebreaker: Any

def find_extensions(
    atom_indices: Any, visited_bond_indices: Any, directed_edges: Any
): ...
def all_subgraph_extensions(
    enumeration_mol: Any,
    subgraph: Any,
    visited_bond_indices: Any,
    internal_bonds: Any,
    external_edges: Any,
) -> None: ...
def find_extension_size(
    enumeration_mol: Any, known_atoms: Any, exclude_bonds: Any, directed_edges: Any
): ...

class CachingTargetsMatcher(dict):
    targets: Any = ...
    required_match_count: Any = ...
    def __init__(
        self, targets: Any, required_match_count: Optional[Any] = ...
    ) -> None: ...
    def shift_targets(self) -> None: ...
    def __missing__(self, smarts: Any): ...

class VerboseCachingTargetsMatcher:
    targets: Any = ...
    cache: Any = ...
    required_match_count: Any = ...
    num_lookups: int = ...
    num_search_true: int = ...
    def __init__(
        self, targets: Any, required_match_count: Optional[Any] = ...
    ) -> None: ...
    def shift_targets(self) -> None: ...
    def __getitem__(self, smarts: Any, missing: Any = ...): ...
    def report(self) -> None: ...

def prune_maximize_bonds(
    subgraph: Any,
    mol: Any,
    num_remaining_atoms: Any,
    num_remaining_bonds: Any,
    best_sizes: Any,
): ...
def prune_maximize_atoms(
    subgraph: Any,
    mol: Any,
    num_remaining_atoms: Any,
    num_remaining_bonds: Any,
    best_sizes: Any,
): ...

class _SingleBest:
    best_num_atoms: int = ...
    best_smarts: Any = ...
    sizes: Any = ...
    timer: Any = ...
    verbose: Any = ...
    def __init__(self, timer: Any, verbose: Any) -> None: ...
    def get_result(self, completed: Any): ...

class MCSResult:
    num_atoms: Any = ...
    num_bonds: Any = ...
    smarts: Any = ...
    completed: Any = ...
    def __init__(
        self, num_atoms: Any, num_bonds: Any, smarts: Any, completed: Any
    ) -> None: ...
    def __nonzero__(self): ...

class SingleBestAtoms(_SingleBest):
    def add_new_match(self, subgraph: Any, mol: Any, smarts: Any): ...

class SingleBestBonds(_SingleBest):
    def add_new_match(self, subgraph: Any, mol: Any, smarts: Any): ...

def check_completeRingsOnly(smarts: Any, subgraph: Any, enumeration_mol: Any): ...

class SingleBestAtomsCompleteRingsOnly(_SingleBest):
    def add_new_match(self, subgraph: Any, mol: Any, smarts: Any): ...

class SingleBestBondsCompleteRingsOnly(_SingleBest):
    def add_new_match(self, subgraph: Any, mol: Any, smarts: Any): ...

def enumerate_subgraphs(
    enumeration_mols: Any,
    prune: Any,
    atom_assignment: Any,
    matches_all_targets: Any,
    hits: Any,
    timeout: Any,
    heappush: Any,
    heappop: Any,
): ...

class Uniquer(dict):
    counter: Any = ...
    def __init__(self) -> None: ...
    def __missing__(self, key: Any): ...

def MATCH(mol: Any, pat: Any): ...

class VerboseHeapOps:
    num_seeds_added: int = ...
    num_seeds_processed: int = ...
    verboseDelay: Any = ...
    trigger: Any = ...
    def __init__(self, trigger: Any, verboseDelay: Any) -> None: ...
    def heappush(self, seeds: Any, item: Any): ...
    def heappop(self, seeds: Any): ...
    def trigger_report(self) -> None: ...
    def report(self) -> None: ...

def compute_mcs(
    fragmented_mols: Any,
    typed_mols: Any,
    minNumAtoms: Any,
    threshold_count: Optional[Any] = ...,
    maximize: Any = ...,
    completeRingsOnly: Any = ...,
    timeout: Any = ...,
    timer: Optional[Any] = ...,
    verbose: bool = ...,
    verboseDelay: float = ...,
): ...

class Timer:
    mark_times: Any = ...
    def __init__(self) -> None: ...
    def mark(self, name: Any) -> None: ...

def fmcs(
    mols: Any,
    minNumAtoms: int = ...,
    maximize: Any = ...,
    atomCompare: Any = ...,
    bondCompare: Any = ...,
    threshold: float = ...,
    matchValences: Any = ...,
    ringMatchesRingOnly: bool = ...,
    completeRingsOnly: bool = ...,
    timeout: Any = ...,
    times: Optional[Any] = ...,
    verbose: bool = ...,
    verboseDelay: float = ...,
): ...
def subgraph_to_fragment(mol: Any, subgraph: Any): ...
def make_fragment_smiles(
    mcs: Any, mol: Any, subgraph: Any, args: Optional[Any] = ...
): ...
def make_fragment_sdf(mcs: Any, mol: Any, subgraph: Any, args: Any): ...
def make_complete_sdf(mcs: Any, mol: Any, subgraph: Any, args: Any): ...

structure_format_functions: Any

def make_structure_format(
    format_name: Any, mcs: Any, mol: Any, subgraph: Any, args: Any
): ...
def parse_num_atoms(s: Any): ...
def parse_threshold(s: Any): ...
def parse_timeout(s: Any): ...

class starting_from:
    left: Any = ...
    def __init__(self, left: Any) -> None: ...
    def __contains__(self, value: Any): ...

range_pat: Any
value_pat: Any

def parse_select(s: Any): ...

compare_shortcuts: Any

def main(args: Optional[Any] = ...) -> None: ...
