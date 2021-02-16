CLI
===

CLI instructions for running spi-dmg.


Examples
--------

The examples covers a simple case where we have output from MethylKit (DMCs) and MethylSig (DMRs).
DMR regions (significant q <= 0.1) with at least 60% of DMCs
(q < 0.1) agreeing with the DMR change in methylation direction are kept.
Genes with multiple DMRs are removed if the DMRs are not in agreement (meth. Diff. direction).
If the DMRs are in agreement, the CpG with the highest DNA methylation difference in the direction of change is
assigned as the methylation value (change and padj) for that gene i.e. as the driver CpG behind the gene’s change in
DNA methylation. Note the cutoff values are all adjustable.

Assign DMR and DMC
------------------

.. code-block:: bash

    scidmg --dmr "data/methylSig_prom.csv" --dmc "data/methyKit_DMC.csv"

Arguments
---------

.. argparse::
   :module: scidmg
   :func: gen_parser
   :prog: scidmg
   :nodefaultconst:
