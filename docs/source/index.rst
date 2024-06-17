Welcome to EvoMol's documentation!
==================================


EvoMol is a Python package for evolutionary molecular design.
It is designed to be a flexible and extensible platform for developing and testing new algorithms for molecular design.

Check out the :doc:`quick_start/quick_start` section for a quick start guide, including the :ref:`installation instructions<installation>`.

EvoMol is originally proposed by `Jules Leguy <https://github.com/jules-leguy/EvoMol>`_ and related to the publications [Leguy2020]_ and [Leguy2021]_, please cite them if you use EvoMol in your research.

.. [Leguy2020] Leguy, J., Cauchy, T., Glavatskikh, M., Duval, B., Da Mota, B. EvoMol: a flexible and interpretable evolutionary algorithm for unbiased de novo molecular generation. J Cheminform 12, 55 (2020). https://doi.org/10.1186/s13321-020-00458-z

.. [Leguy2021] Leguy, J., Glavatskikh, M., Cauchy, T. et al. Scalable estimator of the diversity for de novo molecular generation resulting in a more robust QM dataset (OD9) and a more efficient molecular optimization. J Cheminform 13, 76 (2021). https://doi.org/10.1186/s13321-021-00554-8

.. note::

   This project is under active development.


.. _introduction:

Introduction
------------

EvoMol is a Python package for evolutionary molecular design.
It provides a framework for the design of molecules with desired properties using evolutionary algorithms.
EvoMol is built on top of RDKit, a widely used open-source cheminformatics toolkit.
EvoMol provides a set of tools for molecular design, including:

- Molecular representation
- Fitness functions
- Evolutionary algorithms
- Mutation operators
- Selection strategies



.. toctree::
    :maxdepth: 2
    :hidden:

    installation
    quick_start/quick_start


.. toctree::
    :maxdepth: 2
    :hidden:
    :caption: User Guide

    user_guide


.. toctree::
    :maxdepth: 2
    :hidden:
    :caption: API Reference

    generated/modules

.. toctree::
    :maxdepth: 2
    :hidden:
    :caption: Development Informations

    authors
    contributing
    history




Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


