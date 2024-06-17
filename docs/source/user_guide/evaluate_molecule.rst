Evaluate a molecule
===================


>>> from evomol.representation.molecule import Molecule
>>> from evomol.representation.smiles import SMILES
>>> from evomol.representation.molecular_graph import MolecularGraph
>>> Molecule.id_representation_class = SMILES
>>> Molecule.representations_class = [MolecularGraph]
>>> mol_graph = Molecule('CCO').get_representation(MolecularGraph)
>>> print(mol_graph.smiles)
CCO