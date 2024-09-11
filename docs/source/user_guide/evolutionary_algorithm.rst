Optimization with Evolutionary Algorithms
=========================================

EvoMol uses evolutionary algorithms to optimize molecules.
The evolutionary algorithm is a population-based optimization algorithm that mimics the process of natural selection.
In the case of EvoMol, the population is a set of molecules where at each iteration, the algorithm select molecules to be replaced, molecules to be mutated, and try for each molecule to be replaced to find a better molecule starting from a molecule to mutate.

You can use the script `scripts/evolutionary_algorithm.py` to run the evolutionary algorithm on a molecule.

.. code-block:: bash

    python scripts/evolutionary_algorithm.py

The different steps in the script are:

1. Initialization of the defaults parameters :
    - the accepted atoms, the max number of heavy atoms and the representation class (SMILES and MolecularGraph by default)
    - the default actions spaces (add atom, change bond, cut atom, insert carbon, move group, remove atom, substitute atom)
    - the input and output files for logging
    - random seed

2. Setting the evaluations :
    - list of evaluations and filters to apply on the molecules
    - selection of the evaluation to use for the optimization with `evomol.search.parameters.search_parameters.fitness_criteria` from :class:`.SearchParameters`.

3. Initialization of the population :
    - list the molecules to start with, if not provided or not completed, the algorithm will fill the population with empty molecules
    - create the :class:`.EvolutionaryAlgorithm` object with the parameters, the population, the parameters (population size, max number of generations, max time, number of molecules to replace per generations), the selection strategy, the evaluations and the filters, and the neighborhood strategy.

4. Run the optimization with the `run` method of the :class:`.EvolutionaryAlgorithm` object.
