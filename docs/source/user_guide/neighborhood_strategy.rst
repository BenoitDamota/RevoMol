Neighborhood Strategy
=====================

The neighborhood strategy is used to select a neighbor of a molecule among the different possible actions to apply to the molecule.

You can create your own neighborhood strategy by inheriting from the :class:`.NeighborhoodStrategy` class and implementing the :meth:`.NeighborhoodStrategy.mutate` method. :meth:`.NeighborhoodStrategy.mutate` takes in input a molecule to mutate, a molecule to replace, the evaluations to use, and the set of tabu SMILES of molecules that have already been in the population, and returns either a new molecule or None if no neighbor have been found.

You can try the script `scripts/new_neighborhood_example.py` to see how to create a new neighborhood strategy.


.. code-block:: bash

    python scripts/new_neighborhood_example.py

The mutate method is called in a loop until a neighbor is found. When the neighborhood strategy return None, the molecule to replace will be sent again to the neighborhood strategy with the next molecule to mutate in the order of the selection.
