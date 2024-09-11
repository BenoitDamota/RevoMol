======
EvoMol
======


EvoMol is a Python package for evolutionary molecular design.
It is designed to be a flexible and extensible platform for developing and testing new algorithms for molecular design.
Reimplementation of EvoMol from Jules Leguy : https://github.com/jules-leguy/EvoMol.

Installation
------------

After cloning the repository, you can install EvoMol with the following command:

.. code-block:: bash

    python3.11 -m venv venv
    source venv/bin/activate
    pip install --upgrade pip
    pip install -r requirements_dev.txt

The code should work with Python 3.10, 3.11 and 3.12.

You should now be able to use EvoMol.

Tests
-----

To run the tests, you can use the following command:

.. code-block:: bash

    # to perform all the tests without details
    pytest
    # to list the tests
    pytest -v
    # to run only the tests related to the GCF
    pytest -v -k gcf

You can also use `tox <https://tox.wiki>`_ to run the tests with different versions of Python and also on different platforms (Linux, Windows) :

.. code-block:: bash

    # to run the tests with all the versions of Python
    tox
    # to run the tests with a specific version of Python
    tox -e py310

The available versions are : py310, py311, py312 and also code analysis with flake8, ruff, pylint.
Check the `tox.ini` file for more information.
There is also a mypy environment but rdkit cause some issues with mypy so I removed it from tox.

Documentation
-------------

To generate the documentation, you can use the following commands:

.. code-block:: bash

    # to generate the documentation
    make docs

You can also run the doctests with the following command:

.. code-block:: bash

    # to run the doctests
    cd docs
    make doctest

There is currently no doctests as EvoMol require a lot of configuration just to use the Molecule object.

Code Analysis
-------------

To run the code analysis, you can use the following commands:

.. code-block:: bash

    # to run the import sorter
    isort evomol/ scripts/ tests/
    # to run the formatter
    black evomol/ scripts/ tests/
    # to run ruff
    ruff check --fix evomol/ scripts/ tests/
    # to run the linter
    pylint evomol/ scripts/ tests/
    # to run the linter
    flake8 evomol/ scripts/ tests/
    # to run the type checker
    mypy evomol/ scripts/ tests/

mypy as some troubles with rdkit so it shows a lot of errors, there is no typings and some problems in the stubs.
There is also an error with the use of Queue, mypy wants type hints for the Queue and pylint don't so I let the error in mypy.

Scripts and Examples
--------------------

You can find some scripts and examples in the `scripts` and `notebooks` directories.

`scripts` :

- `cache_test.py` : short example on the use of the cache decorator with python and joblib. No parameters.

    .. code-block:: bash

        python scripts/cache_test.py

- `draw_mol.py` : example on how to draw a molecule, add annotations and draw multiple molecules in a grid. No parameters.

    .. code-block:: bash

        python scripts/draw_mol.py

- `explore_neighborhood.py` : For a molecule given in parameters, explore the neighborhood with molecular graph neighborhood. The molecule is given in SMILES format. It also check the validity of the neighborhood and the cardinality of the neighborhood.

    .. code-block:: bash

        python scripts/explore_neighborhood.py "COC(C)C(O)"

- `molecule_validity.py` : Given a SMILES string, check if the molecule is valid and print the number of unknown ECFP4 fingerprints and GCF (Generic Cyclic Features).

    .. code-block:: bash

        python scripts/molecule_validity.py "COC(C)C(O)"

- `compare_datasets.py` : Compare the dataset against the enumerations from EvoMol. Read comments in the script for more information.

- `enumerate_smiles.py` : Enumerate all reachable molecules from a given SMILES string. Read comments in the script for more information.

- `find_closest_neighbor.py` : Find the closest neighbor of a molecule in the neighborhood. The script takes 3 parameters: the molecule in SMILES format, the dataset to use to filter (chembl or chembl_zinc), and the maximum number of atoms in the molecule.

    .. code-block:: bash

        python scripts/find_closest_neighbor.py "N1=S=NC2=C1N=S=N2" chembl 10

- `neighborhood_validity_cardinality.py` : enumerate neighborhood of a molecule and check the validity and cardinality of the neighborhood. The script takes in parameter the molecule in SMILES format and prints the results in the console in csv format.

    .. code-block:: bash

        python scripts/neighborhood_validity_cardinality.py "CS(=O)(=O)C1=C(O)C(CC2=CC=C(F)S2)=CC(F)=C1"

- `convert_json_to_csv.py` : convert dataset (zinc, chembl, evo10, ...) from json to csv.

- `evaluation.py` : short example on the use of the evaluation function. No parameters. The script will print the results of the evaluation and show an example of an error in the order of the evaluations.

    .. code-block:: bash

        python scripts/evaluation.py

- `old_new_fingerprints.py` : compare the old and new way of doing fingerprints with rdkit.

- `generate_ecfp_file.py` : generate a file with the ECFP4 fingerprints from a list of SMILES. You have to implement the extraction of the SMILES from your dataset.

    .. code-block:: bash

        python scripts/generate_ecfp_file.py

- `generate_gcf_file.py` : generate a file with the GCF fingerprints from a list of SMILES. You have to implement the extraction of the SMILES from your dataset.

    .. code-block:: bash

        python scripts/generate_gcf_file.py

- `solution_to_remove_error_rdkit.py` : proposition to hide error messages from rdkit.

- `evolutionary_algorithm.py` : example of an evolutionary algorithm to optimize QED. No parameters.

    .. code-block:: bash

        python scripts/evolutionary_algorithm.py

- `logging_tests.py` : shows how the logging works in python.

    .. code-block:: bash

        python scripts/logging_tests.py
- `tanimoto.py` : begining of an attempt to replace the tanimoto function from guacamol. You will need to install guacamol to run this script. See the warning in the installation page of the docs on how to install guacamol.

- `new_representation_example.py` : example of the creation of a new representation. No parameters.

    .. code-block:: bash

        python scripts/new_representation_example.py

- `new_evaluation_example.py` : example of the creation of a new evaluation. No parameters.

    .. code-block:: bash

        python scripts/new_evaluation_example.py

- `new_mutation_example.py` : example of the creation of a new mutation. No parameters.

    .. code-block:: bash

        python scripts/new_mutation_example.py
