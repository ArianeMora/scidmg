# sci-dmg
[![codecov.io](https://codecov.io/github/ArianeMora/scidmg/coverage.svg?branch=main)](https://codecov.io/github/ArianeMora/scidmg?branch=main)
[![PyPI](https://img.shields.io/pypi/v/scidmg)](https://pypi.org/project/scidmg/)
[![DOI](https://zenodo.org/badge/316410924.svg)](https://zenodo.org/badge/latestdoi/316410924)

[Link to docs](https://arianemora.github.io/scidmg/)

## Differentially Methylated Genes
This package aims to consolidate DMRs (differentially methylated regions) and DMCs (differentially methylated cytosines)
and most importantly, develop a consistent, **unbiased** method for assigning a change in DNA methylation to a gene.
See the [docs](https://arianemora.github.io/scidmg/) for detail.

The user provides a DMR file, a file with the percentage of DNA Methylation, and also the DMCs. Using these, sci-DMG
consolidates the DMR's and DMC's that are consistent. DMR regions (significant q <= 0.1) with at least 60% of DMCs
(q < 0.1) agreeing with the DMR change in methylation direction were kept.
Genes with multiple DMRs were removed if the DMRs were not in agreement (meth. Diff. direction).
If the DMRs were in agreement, the CpG with the highest DNA methylation difference in the direction of change is
assigned as the methylation value (change and padj) for that gene i.e. as the driver CpG behind the gene’s change in
DNA methylation. Note the cutoff values are all adjustable.

Any tool can be used to produce the DMC's and DMR's, two such tools are
[MethylKit](https://bioconductor.org/packages/release/bioc/html/methylKit.html) and
[MethylSig](https://pubmed.ncbi.nlm.nih.gov/24836530/) many others exist.

**Note:** The CpGs and DMRs must have already been assigned to genes, this tool consolidates the DMRs and DMCs after
they have been assigned to genes. 

### Example format for Methylsig
``` 
idx,seqnames,start,end,gene_idx,meth_diff,uid,pvalue,fdr,ensembl_gene_id,external_gene_name,chromosome_name,start_position,end_position,strand
1,chr1,10,100,1,-0.6,dmr_1,0.001,0.01,ENSG00000278267,AC114488.2,chr1,1,30,-1
``` 

### Example format for MethylKit
``` 
idx,chr,start,end,gene_idx,meth.diff,uid,pvalue,qvalue,ensembl_gene_id,external_gene_name,chromosome_name,start_position,end_position,strand
1,chr1,1,1,1,-0.6,dmc_1,0.001,0.01,ENSG00000278267,AC114488.2,chr1,1,1,-1
``` 

### Example output

```
--------------------------------------------------------------------------------
                  Length of all merged methylation data: 	16	                   
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
            Length of merged methylation data grouped by region: 	6	            
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
                Number of CpGs to keep based on the regions: 	4	                
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
Length of filtered methylation dataframe: 	7	
Number of genes with Methylation:	3	
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
                 Dropping any genes with disagreeing DMRs: 	1	                  
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
     Length of dataframe filtered to only keep top DMC mapped to genes:	2	      
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
                               Printing stats: 	                                
--------------------------------------------------------------------------------
Length grouped by DMRs 6
Number of CpGs to keep from grouped DMRs 4
Length of merged DMR and DMC 16
Length grouped by Genes 3
Number of Genes with DMRs that disagree in direction 1
```

### Example output file
```
,idx_dmc,chr,start_dmc,end_dmc,gene_idx_dmc,meth.diff,uid_dmc,pvalue_dmc,qvalue,ensembl_gene_id_dmc,external_gene_name,chromosome_name_dmc,start_position_dmc,end_position_dmc,strand_dmc,gene_dmc-uid_dmc,WT_1,WT_2,WT_3,KO_1,KO_2,KO_3,idx_dmr,seqnames,start_dmr,end_dmr,gene_idx_dmr,meth_diff,uid_dmr,pvalue_dmr,fdr,ensembl_gene_id_dmr,chromosome_name_dmr,start_position_dmr,end_position_dmr,strand_dmr,gene_dmc-uid_dmr,CpGsToKeep,abs_logfc
14,7,chr1,1,1,7,-0.1,dmc_8,0.01,0.001,ENSG00000116273,HOXA1,chr4,2,12,1,HOXA1 dmc_8,0.98,0.9,0.64,0.18,0.06,0.51,5,chr1,123,190,2,-0.1,dmr_5,0.01,0.001,ENSG00000116273,chr4,2,12,1,HOXA1 dmr_5,1,0.1
```

### Other
Please post questions and issues related to sci-dmg on the `Issues <https://github.com/ArianeMora/scidmg/issues>`_  section of the GitHub repository.

