.. _installation:

Installation
------------

To install EvoMol, you can use the following command:

.. code-block:: bash

    pip install evomol

.. warning::
    `GuacaMol <https://github.com/BenevolentAI/guacamol>`_ currently use an old version of scipy that had a `histogram` function that is now deprecated.
    To fix this, once installed in your python environment, in the file `venv/lib/pythonx.y/site-packages/guacamol/utils/chemistry.py` you need to replace the line 10 `from scipy import histogram` with `from numpy import histogram`.

Evaluation functions
""""""""""""""""""""

To use some of the Evaluation functions, you will need external data. You can create a directory called `external_data` in the root directory of the repository and download the data there.

.. code-block:: bash

    mkdir external_data
    cd external_data

To use CLScore, you need to download the external data from the `GDBChEMBL repository <https://github.com/reymond-group/GDBChEMBL>`_ :

.. code-block:: bash

    wget https://github.com/reymond-group/GDBChEMBL/raw/master/shingle_libs/chembl_24_1_shingle_scores_log10_rooted_nchir_min_freq_100.pkl



To use RDFilter, you need to download the external data from the `rd_filters repository <https://github.com/PatWalters/rd_filters>`_ :

.. code-block:: bash

    wget https://github.com/PatWalters/rd_filters/raw/master/rd_filters/data/alert_collection.csv
    wget https://github.com/PatWalters/rd_filters/raw/master/rd_filters/data/rules.json

To use SillyWalks, you need to download the external data from the there :


.. code-block:: bash

    # TODO add link to ecfp4_CHEMBL_ZINC.txt
    # TODO add link to ecfp4_CHEMBL.txt

To use generic cyclic features or scaffolds, you need to download the external data from the there :

.. code-block:: bash

    # TODO add link to gcf_CHEMBL_ZINC.txt
    # TODO add link to gcf_CHEMBL.txt

