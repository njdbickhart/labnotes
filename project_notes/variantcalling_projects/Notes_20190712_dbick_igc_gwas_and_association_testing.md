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
* /project/rumen_longread_metagenome_assembly/kiranmayee/IGC/gmmat <- gmmat data on CERES
* /project/rumen_longread_metagenome_assembly/kiranmayee/IGC/prelim_gwas/round2/data_files/ <- data files for GWAS

And Kiranmayee's note file for [GMMAT (best model)](https://github.com/bkiranmayee/My_Labnotes/blob/master/IGC/second_round_GMMAT_GWAS.md) and [GEMMA (faulty assumption)](https://github.com/bkiranmayee/My_Labnotes/blob/master/IGC/second_round_GEMMA_gwas.md) are at the links. 

Let's see if I can combine all of the markers from our 1st and second rounds and analyze them using plink. First, the case-control analysis. 

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/kiranmayee/IGC/case_control

```bash
module load plink/1.9

# Moving and merging files
cp ../neogen_results/round_1/neogen_updated.map ./
cp ../neogen_results/round_1/neogen_updated.ped ./
cp ../prelim_gwas/round2/data_files/pheno_cov_master_file.txt ./

plink --file neogen_updated --merge ../neogen_results/round_2/neogen_2ndRound --out neogen_merged_unfiltered --allow-extra-chr

# Calculating summary statistics
plink --bfile neogen_merged_unfiltered --missing --out neogen_merged_unfiltered.miss --allow-extra-chr

# Filtering variants
plink --bfile neogen_merged_unfiltered --allow-extra-chr --geno 0.20 --make-bed --out neogen_merged_genopop_filt
plink --bfile neogen_merged_genopop_filt --hardy --allow-extra-chr --out neogen_merged_hardy
plink --bfile neogen_merged_genopop_filt --hwe 1E-10 --make-bed --allow-extra-chr --out neogen_merged_hwe_filt
#plink --bfile neogen_merged_hwe_filt --maf 0.05 --make-bed  --allow-extra-chr --out neogen_merged_hwe_maf_filt

# running the association using logistic regression
plink --bfile neogen_merged_hwe_filt --logistic --genotypic --allow-extra-chr --covar pheno_cov_master_file.txt --hide-covar --out neogen_merged_logistic 

# Didn't work well. Going to reformat the covariate file because I think that was screwed up
perl -lane 'if($F[0] eq "id"){next;} print "$F[0]\t$F[1]\t$F[5]\t$F[6]\t$F[8]";' < pheno_cov_master_file.txt > pheno_cov_master_file.crop.tab

plink --bfile neogen_merged_hwe_filt --logistic --genotypic --allow-extra-chr --covar pheno_cov_master_file.crop.tab --hide-covar --out neogen_merged_logistic 
```

OK, the logistic regression wasn't fantastic and I think that it was primarily due to how the files were formatted. Before I start getting fancy with the data, let's start from the beginning and go back to basics. 

## Replicating the other studies

Let's try to replicate [Bermingham et al 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3998787/) and [Wilkinson et al. 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5483290/) first. I think that Bermingham is far too narrow and does not fit ALL of the data. Instead, the approach by Wilkinson et al using regional heritability estimates may be the best approach. Let's see how it deals with SNPs from our new set that have reasonable heterozygosity and minor allele frequency.

### Generating the new SNPs for use in the model

```bash

```

