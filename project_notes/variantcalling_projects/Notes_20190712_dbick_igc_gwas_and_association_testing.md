# IGC GWAS and association testing
---
*7/12/2019*

These are my notes on performing association testing and curating Kiranmayee's GWAS results for subsequent publication.

## Table of Contents


## Where we left off

I need to organize the details so that we can put together a coherent publication. My thoughts are to perform association testing (via plink using Cochran-Armitage after pruning SNPs with extreme deviations from HW) with just the custom genotype SNPs, then expand the analysis to a GWAS with the HD chip SNPs and custom SNPs that survived HW pruning. I will prune the larger HD + custom dataset for LD and use in subsequent GWAS.

Kiranmayee has already preformed the logistic regression analysis and identified a lambda value that indicates that population stratification does not influence the analysis of the data. 

**Here are the locations of the data:**

* /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/neogen_results
* /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/neogen_results/round_2  <- plink files for association testing and GWAS
* /mnt/nfs/nfs2/bickhart-users/natdb_sequencing/prelim_gwas/round2/gmmat  <- kiranmayee's gmmat data
* /beegfs/project/rumen_longread_metagenome_assembly/kiranmayee/IGC/gmmat <- gmmat data on CERES

And Kiranmayee's note file for [GMMAT (best model)](https://github.com/bkiranmayee/My_Labnotes/blob/master/IGC/second_round_GMMAT_GWAS.md) and [GEMMA (faulty assumption)](https://github.com/bkiranmayee/My_Labnotes/blob/master/IGC/second_round_GEMMA_gwas.md) are at the links. 