"""
This module contains the ActionMolGraph class, which is a subclass of the
Action class. It is an abstract class that defines the interface for actions
that modify a molecular graph.
It is used to load the molecular graph from the molecule, call the apply_action
method from the subclass, and update the molecule with the new molecular graph.

Warning: the new molecular graph must be modified in place in the apply_action
method.
"""

from abc import abstractmethod
from copy import copy

from typing_extensions import override

from evomol.action import Action, ActionError
from evomol.representation import MolecularGraph, Molecule


class ActionMolGraph(Action):
    """
    Abstract class for actions that modify a molecular graph.
    """

    @override
    def _apply(self) -> Molecule:
        """Apply the action to the molecule.

        Raises:
            ActionError: Error if the action cannot be applied.

        Returns:
            Molecule: Molecule after the action.
        """
        # get the molecular graph representation
        mol_graph: MolecularGraph = self.molecule.get_representation(MolecularGraph)
        assert mol_graph is not None
        # make a copy of the molecular graph
        new_mol_graph: MolecularGraph = copy(mol_graph)

        # edit the new molecular graph in place
        self.apply_action(new_mol_graph)

        # update the representation of the molecule
        try:
            new_mol_graph.update_representation()
        except Exception as e:
            # raise an error if the new molecular graph cannot be converted to
            # a molecule should not happen
            # if it does, it is a bug in the code of the actions, either in the
            # apply_action method or in the list_actions method
            # it would mean that the molecular graph is not a valid molecule
            # in terms of valence, etc.
            raise ActionError(self, new_mol_graph.canonical_smiles, repr(e)) from e

        return Molecule(new_mol_graph.canonical_smiles)

    @abstractmethod
    def apply_action(self, new_mol_graph: MolecularGraph) -> None:
        """Subclass must implement this method to apply the action to the
        molecular graph.

        Args:
            new_mol_graph (MolecularGraph): Molecular graph to which the action
            is applied. The method should modify this object in place.
        """

    @override
    def representation_name(self) -> str:
        return MolecularGraph.class_name()
