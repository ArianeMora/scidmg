*******
sci-dmg
*******


Sci-dmg aims to assign a change in DNA methylation (as calculated by an external tool) to genes in a consistent and
 unbiased manner.

The user provides a DMR file, a file with the percentage of DNA Methylation, and also the DMCs. Using these, sci-DMG
consolidates the DMR's and DMC's that are consistent. DMR regions (significant q <= 0.1) with at least 60% of DMCs
(q < 0.1) agreeing with the DMR change in methylation direction were kept.
Genes with multiple DMRs were removed if the DMRs were not in agreement (meth. Diff. direction).
If the DMRs were in agreement, the CpG with the highest DNA methylation difference in the direction of change is
assigned as the methylation value (change and padj) for that gene i.e. as the driver CpG behind the gene’s change in
DNA methylation. Note the cutoff values are all adjustable.

Any tool can be used to produce the DMC's and DMR's, two such tools are
<MethylKit https://bioconductor.org/packages/release/bioc/html/methylKit.html>_ and
<MethylSig https://pubmed.ncbi.nlm.nih.gov/24836530/>_ many others exist.

Running sci-dmg
===============

1. Install sci-dmg (:ref:`Installing <installing>`)

2. View a python examples in (:ref:`examples <examples/python_notebook>`)

3. Look at CLI examples in (:ref:`cli <examples/cli>`)

Extending sci-dmg
=================

1. Make a pull request on github.

Citing sci-dmg
==============
Sci-DMG can be cited as in :ref:`references`, where we also provide citations for the used tools (e.g. numpy).

.. toctree::
   :caption: Getting started
   :maxdepth: 1

   about
   installing/index


.. toctree::
   :caption: Running sci-dmg
   :maxdepth: 1

   examples/python_notebook
   examples/cli


.. toctree::
   :caption: About
   :maxdepth: 1

   faq
   changelog
   references
