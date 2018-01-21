.. _prsice-quick-start:

PRsice Quick Start
******************

===========
Preparation
===========
To perform PRSice, you will need

------
Input
------
- PRSice.R file: A wrapper for the PRSice binary and for plotting
- PRSice executable file: Perform all analysis except plotting
- Base data set: GWAS summary results, which the PRS is based on
- Target data set: Raw genotype data of **target phenotype**. Can be in the form of  `PLINK binary <https://www.cog-genomics.org/plink2/formats#bed>`_ or `BGEN <http://www.well.ox.ac.uk/~gav/bgen_format/>`_

Note: You should first perform quality control on your genotype file before passing it to PRSice
You can do that using PLINK. A good starting point is (assume **($target)** is the prefix of your target binary file)

.. code-block:: bash

    plink --bfile ($target) \
      --maf 0.05 \
      --mind 0.1 \
      --geno 0.1 \
      --hwe 1e-6 \
      --make-just-bim \
      --make-just-fam \
      --out ($target).qc

Then you can add :code:`--keep ($target).qc.fam --extract ($target).qc.bim` to PRSice command to filter out
the samples and SNPs

.. note::
   This is only a short example. For detail information and tutorial, please refer to our [step by step tutorial](step-by-step-tutorial)*

===============
Running PRSice
===============
In most case, you can simply run PRSice using the following command, assuming your
PRSice executable is located in :code:`($HOME)/PRSice/` and you are working in :code:`($HOME)/PRSice`

.. note::
   For window users, please use **Rscript.exe** instead of **Rscript**

--------------
Binary Traits
--------------
For binary traits, you can use the following command

.. code-block:: bash

  Rscript PRSice.R --dir . \
    --prsice ./PRSice \
    --base TOY_BASE_GWAS.assoc \
    --target TOY_TARGET_DATA \
    --thread 1 \
    --stat OR \
    --binary-target T


--------------------
Quantitative Traits
--------------------

For quantitative traits, you can use

.. code-block:: bash

  Rscript PRSice.R --dir . \
    --prsice ./PRSice \
    --base TOY_BASE_GWAS.assoc \
    --target TOY_TARGET_DATA \
    --thread 1 \
    --stat BETA \
    --binary-target F


.. note::
  The default of PRSice is based on the header of your base file:
  1. When *BETA* (case insensitive) is found in the header and :code:`--stat` was not provided, :code:`--beta` will be added to your command, and if :code:`--binary-target` was not provided, :code:`--binary-target F` will be added to your command
  2. When *OR* (case insensitive) is found in the header and :code:`--binary-target` was not provided, :code:`--binary-target T` will be added to your command

  Most importantly, PRSice will detail all effective option in its log file where you can simply copy and paste it to get the same output

  You can also disable this behaviour by using the :code:`--no-default` option
