*******
sci-dmg
*******


Sci-dmg aims to assign a change in DNA methylation (as calculated by an external tool) to genes in a consistent and
 unbiased manner.

The user provides a DMR file, a file with the percentage of DNA Methylation, and also the DMCs. Using these, sci-DMG
consolidates the DMR's and DMC's that are consistent. DMR regions (significant q <= 0.1) with at least 60% of DMCs
(q < 0.1) agreeing with the DMR change in methylation direction are kept.
Genes with multiple DMRs are removed if the DMRs are not in agreement (meth. Diff. direction).
If the DMRs are in agreement, the CpG with the highest DNA methylation difference in the direction of change is
assigned as the methylation value (change and padj) for that gene i.e. as the driver CpG behind the gene’s change in
DNA methylation. Note the cutoff values are all adjustable.
