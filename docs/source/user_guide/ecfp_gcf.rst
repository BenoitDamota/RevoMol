Filters Unknown ECFP and GCF
============================

EvoMol will propose molecules that are valid in terms of valence and other "simple" chemical properties but are impossible to find in the nature or to synthesize.
To remove those impossible molecules, we propose filters that will check if parts of the molecule are valid based on ECFP and GCP (Generic Cyclic Features).

To do that, we use the classes :class:`.UnknownGCF`, :class:`.FilterUnknownGCF`, :class:`.UnknownECFP`, and :class:`.FilterUnknownECFP` to count the number of unknown ECFP and GCF then filter the molecules that have more than a certain threshold of unknown ECFP or GCF.

The classes :class:`.UnknownGCF` and :class:`.UnknownECFP` ask for a file containing a list of ECFP or GCF that are known and valid.
For each molecule, the class will extract the ECFP or GCF and check if they are in the list of known ECFP or GCF.
The classes :class:`.FilterUnknownGCF` and :class:`.FilterUnknownECFP` ask for a threshold of unknown ECFP or GCF.
It will raise an :class:`.EvaluationError` if the molecule has more unknown ECFP or GCF than the threshold.

Create your own filters on ECFP and GCF
---------------------------------------

To create your own filters on ECFP and GCF, you need to create new files of ECFP and GCF that are allowed for the molecules you want to generate.
To help you with that you can use the scripts `scripts/generate_ecfp.py` and `scripts/generate_gcf.py` that will generate the files for you, you just need to provide the SMILES of molecules that are allowed.

Once you have the files, when creating the :class:`.UnknownGCF` or :class:`.UnknownECFP` class, you need to provide the path to the file.

For example, the following code shows how to create a filter that will check with the files `my_ecfp1.txt` and `my_gcf1.txt` that the molecule is valid:

.. code-block:: python

    from evomol.evaluation import UnknownECFP, FilterUnknownECFP, UnknownGCF, FilterUnknownGCF, is_valid_molecule


    # make sure the name is the same for those two classes
    evaluations = [
        UnknownECFP(path_db="my_ecfp1.txt", name="my_ecfp1"),
        FilterUnknownECFP(threshold=0, name="my_ecfp1"),
        UnknownGCF(path_db="my_gcf1.txt", name="my_gcf1"),
        FilterUnknownGCF(threshold=0, name="my_gcf1"),
    ]

    molecule = Molecule("C1=C=C=1")

    print(is_valid_molecule(molecule, evaluations))

The name given to the classes :class:`.UnknownECFP` and :class:`.FilterUnknownECFP` must be the same as it is used in the filter (and the same for GCF).
