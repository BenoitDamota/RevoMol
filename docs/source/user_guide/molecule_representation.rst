Molecule Representation
=======================

In EvoMol, a molecule is represented by a :class:`.Molecule` object that contains different representations.

Each representation is a subclass of the :class:`.MoleculeRepresentation` class.

That class that contains the information of the molecule in a specific format.
It can also propose actions to modify the molecule, such as adding an atom or changing a bond.
For example, the :class:`.MolecularGraph` representation contains a graph representation of the molecule with the help of the rdkit library, while the :class:`.SMILES` representation contains the SMILES string of the molecule.


Create a new representation
---------------------------

To create a new representation, you need to create a new class that inherits from the :class:`.MoleculeRepresentation` class.
The class must implement the `representation` method that returns the representation of the molecule in a string format.
The `__eq__` method must be implemented to compare two molecules.

For example, the following code shows a new representation called `ECFP4` that uses the Extended Connectivity Fingerprint 4 (ECFP4) to represent a molecule, the complete code is available in `scripts/new_representation_example.py`.

.. code-block:: python


    from rdkit import Chem
    from rdkit.Chem import AllChem
    from typing_extensions import override

    from evomol.representation import MoleculeRepresentation, SMILES, Molecule


    class ECFP4(MoleculeRepresentation):
        """
        Molecule representation using the Extended Connectivity Fingerprint 4 (ECFP4).
        """

        def __init__(self, smiles: str):
            super().__init__(smiles)

            mol = Chem.MolFromSmiles(smiles)

            # create a Morgan fingerprints generator
            self.fingerprints_generator = AllChem.GetMorganGenerator(radius=2)
            fingerprints = self.fingerprints_generator.GetSparseCountFingerprint(
                mol
            ).GetNonzeroElements()

            self.fingerprints = set(fingerprints.keys())

        @override
        def representation(self) -> str:
            return str(self.fingerprints)

        def __eq__(self, value: object) -> bool:
            if not isinstance(value, ECFP4):
                return False
            return self.fingerprints == value.fingerprints

Once the new representation is created, it can be used to create a new molecule object.
You can also create actions on this representation, see :ref:`the mutation page<mutation>` for more informations.

You can then use the new representation in the molecule object by setting the `representations_class` attribute of the :class:`.Molecule` class.
When the molecule object is created, it will automatically create each representation and it can be accessed by the `get_representation` method that require the class of the representation as argument and return the representation object.

.. code-block:: python

    # initialize the molecule representation
    Molecule.id_representation_class = SMILES
    Molecule.representations_class = [ECFP4]
    Molecule.max_heavy_atoms = 6

    Molecule.accepted_atoms = ["C", "O", "N", "S"]

    smiles = "CCO"

    molecule = Molecule(smiles)

    print(molecule.get_representation(ECFP4).representation())
