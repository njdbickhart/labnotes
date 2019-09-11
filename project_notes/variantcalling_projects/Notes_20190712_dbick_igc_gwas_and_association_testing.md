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

Actually, Wilkinson used a proprietary program that isn't easy to install. At all. Very few directions and hard-coded makefile instructions based on system ID names. I'll try to replicate Kiranmayee's data instead using different filtered variant sites.

I am going to re-filter the custom markers to remove those with excessive heterozygosity, low MAF, and keep the remainder. 

```bash
module load plink/1.9
module load r

plink --bfile neogen_merged_genopop_filt --freqx --allow-extra-chr --out neogen_merged_genopop_stats

# Now I processed the file manually to check for distributions
perl -lane 'if($F[0] eq "CHR"){next;} $maf = $F[6] / ($F[4] + $F[5] + $F[6]); $het = $F[5] / ($F[4] + $F[5] + $F[6]); print "$F[1]\t$maf\t$het";' < neogen_merged_genopop_stats.frqx | perl -lane 'print $F[-1]' | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   115
Sum:    26.9967031417674
Minimum 0
Maximum 1
Average 0.234754
Median  0.179294389820706
Standard Deviation      0.237575

# I will select a cutoff of 50% for the het rate and 95% for the MAF
perl -lane 'if($F[0] eq "CHR"){next;} $maf = $F[6] / ($F[4] + $F[5] + $F[6]); $het = $F[5] / ($F[4] + $F[5] + $F[6]); print "$F[1]\t$maf\t$het";' < neogen_merged_genopop_stats.frqx | perl -lane 'if($F[1] <= 0.95 && $F[2] < 0.50){print $F[0];}' > non_het_snps.txt

plink --bfile neogen_merged_genopop_filt --extract non_het_snps.txt --allow-extra-chr --make-bed --out neogen_merged_nonhet_filt

# Just testing for curiosity's sake
plink --bfile neogen_merged_nonhet_filt --logistic --genotypic --allow-extra-chr --covar pheno_cov_master_file.crop.tab --hide-covar --out neogen_merged_nonhet_logistic

perl -lane 'if($F[-1] ne "NA" && $F[-1] < 0.005){print $_;}' < neogen_merged_nonhet_logistic.assoc.logistic
 CHR                    SNP         BP   A1       TEST    NMISS         OR         STAT            P
   5             5_99190989   99190989    T   GENO_2DF     1601         NA         21.6    2.039e-05
  18            18_62774779   62774779    A     DOMDEV     1763     0.5293       -3.694    0.0002204
  18            18_62774779   62774779    A   GENO_2DF     1763         NA        18.73    8.557e-05
  18            18_63417698   63417698    A   GENO_2DF     1678         NA        18.45    9.872e-05
 MHC             MHC_260913     260913    G        ADD     1654      1.207        2.916     0.003545
LIB14427_MHC     LIB14427_MHC_73766      73766    A     DOMDEV     1736      2.027         3.85    0.0001179
LIB14427_MHC     LIB14427_MHC_73766      73766    A   GENO_2DF     1736         NA        15.13    0.0005195

# OK, now to trace Kiranmayee's steps to rerun the GWAS
# I need to convert the coordinates so that the SNPs can be merged properly
grep 'ARS_PIR' ../prelim_gwas/round2/data_files/merge1.bim > previously_reconciled_snps.tab

# Running a union analysis to see how many overlap with Kiranmayee's previous coordinates
perl -lane 'print $F[1]' < neogen_merged_nonhet_filt.bim > new_filt.list
perl -lane '@b = split(/_/, $F[1]); if(scalar(@b) == 4 || $F[1] =~ /LRC/ || $F[1] =~ /KIR/){print "$b[-2]\_$b[-1]";}else{print "$b[-3]\_$b[-2]\_$b[-1]";}' < previously_reconciled_snps.tab > previously_reconciled_snps.list

perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl new_filt.list previously_reconciled_snps.list
File Number 1: new_filt.list
File Number 2: previously_reconciled_snps.list
Set     Count
1       17
1;2     55
2       25

# Damn! Well, in order to convert the coordinates, I will use a regression based on other members of the group
perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o new_filt.list previously_reconciled_snps.list

perl -lane '@b = split(/_/, $F[1]); if(scalar(@b) == 4 || $F[1] =~ /LRC/ || $F[1] =~ /KIR/){$F[1] = "$b[-2]\_$b[-1]"; print join("\t", @F);}else{$F[1] = "$b[-3]\_$b[-2]\_$b[-1]"; print join("\t", @F);}' < previously_reconciled_snps.tab | perl -lane '@b = split(/_/, $F[1]); $old = pop(@b); $val = join("_", @b); print "$F[0]\t$val\t$F[1]\t$old\t$F[3]\t$F[4]\t$F[5]";' > previously_reconciled_snps.training.tab
```

Now to do the conversion in R

```R
library(dplyr)
traindf <- read.delim("previously_reconciled_snps.training.tab", header=F)
names(traindf) <- c("chr", "set", "ID", "old", "pos", "ref", "alt")
preddf <- read.delim("neogen_merged_nonhet_filt.bim", header=F)
names(preddf) <- c("set", "ID", "sex", "old", "ref", "alt")

preddf$set <- as.character(preddf$set)
preddf$set[preddf$set == "MHC"] <- "MHCclassI_MHC"
preddf$set <- as.factor(preddf$set)

# Creating the models
models = list(); for(x in unique(preddf$set)){fit <- lm(pos ~ old, data=traindf[traindf$set == x,]); models[[x]] <- fit}

# now to convert coords
for (x in names(models)){v <- predict(models[[x]], preddf[preddf$set == x,]); preddf[preddf$set == x,"pos"] <- v;}

# dammit! It didn't work very well
# I found out why:
traindf %>% mutate(diff = pos - old) %>% group_by(set) %>% summarize(mean = mean(diff), stdev = sd(diff))
# A tibble: 7 x 3
  set                 mean   stdev
  <fct>              <dbl>   <dbl>
1 18               255796.  65849.
2 23               -59833.   9038.
3 5                484491   18554.
4 KIR            62932209      NA
5 LIB14427_MHC   28456035.  13267.
6 LRC            63107830. 156684.
7 MHCclassI_MHC  28259399. 120393.

# It's gotta be because of gaps in the region and Kiranmayee's liftover. The linear regression is probably fine as is for most values?
preddf %>% mutate(diff = pos - old) %>% group_by(set) %>% summarize(mean = mean(diff), stdev = sd(diff))
# A tibble: 7 x 3
  set                 mean  stdev
  <fct>              <dbl>  <dbl>
1 18               256757. 23365.
2 23               -64878.   271.
3 5                483395.  2788.
4 KIR            62932209     NA
5 LIB14427_MHC   28457215. 11947.
6 LRC            63103197. 38638.
7 MHCclassI_MHC  28274157. 46393.

# Yeah, not the best, but less dispersion than the original values and less than using the mean in a simple subtraction.
write.table(x=preddf, file="neogen_merged_nonhet_filt.revised.bim",row.names=FALSE, col.names=FALSE, quote=FALSE)
```

I then manually edited the file to change the chromosome names and positions.

```bash
perl -lane '$F[-1] = int($F[-1]); if($F[0] =~ /KIR/){$F[0] = 18;}elsif($F[0] =~ /MHC/){$F[0] = 23;}elsif($F[0] =~ /LRC/){$F[0] = 18;} print "$F[0]\t$F[1]\t$F[2]\t$F[6]\t$F[4]\t$F[5]";' < neogen_merged_nonhet_filt.revised.bim > neogen_merged_nonhet_alter.bim
cp neogen_merged_nonhet_filt.bed neogen_merged_nonhet_alter.bed
cp neogen_merged_nonhet_filt.fam neogen_merged_nonhet_alter.fam

# Merging together
plink --bfile ../prelim_gwas/round1/trial_1/set1_ld_filtered --bmerge neogen_merged_nonhet_alter --cow --out hd_neogen_merged_snpset

# Ah, crap! Kiranmayee must have done something to this file. I'll have to replicate from her later results
grep 'ARS_PIR' ../prelim_gwas/round2/data_files/merge1.bim | perl -lane 'print $F[1];' > markers_to_remove.txt
plink --bfile ../prelim_gwas/round2/data_files/merge1 --cow --exclude markers_to_remove.txt --make-bed --out hd_array_ld_filt
plink --bfile hd_array_ld_filt --bmerge neogen_merged_nonhet_alter --cow --out hd_neogen_merged_snpset

# OK, I have the SNPs, let's now try to replicate her phenotype data sorting
mkdir pheno_groups
plink --bfile hd_neogen_merged_snpset --cow --keep-fam ../prelim_gwas/round2/data_files/pheno1.fam --make-bed --out pheno_groups/phenotype1
plink --bfile hd_neogen_merged_snpset --cow --keep-fam ../prelim_gwas/round2/data_files/pheno2.fam --make-bed --out pheno_groups/phenotype2
plink --bfile hd_neogen_merged_snpset --cow --keep-fam ../prelim_gwas/round2/data_files/pheno3.fam --make-bed --out pheno_groups/phenotype3
plink --bfile hd_neogen_merged_snpset --cow --keep-fam ../prelim_gwas/round2/data_files/pheno4.fam --make-bed --out pheno_groups/phenotype4

# Changing her phenotype codes to fit with plink
for i in `seq 1 4`; do echo $i; perl -lane 'if($F[-1] == -9){$F[-1] = 0;} print join(" ", @F);' < pheno_groups/phenotype${i}.fam > temp; mv temp pheno_groups/phenotype${i}.fam; done

# And, running a quick logistic analysis
for i in `seq 1 4`; do echo $i; plink --bfile pheno_groups/phenotype${i} --cow --1 --logistic --genotypic --covar pheno_cov_master_file.crop.tab --hide-covar --out pheno_groups/logistic_test_pheno${i}; done

# Making a separate folder for analysis
mkdir gmmat_output
mkdir gmmat_data
```

OK, I need to create GDS files for each dataset for use in GMMAT.

```R
library(SeqArray)
seqBED2GDS("pheno_groups/phenotype1.bed", "pheno_groups/phenotype1.fam", "pheno_groups/phenotype1.bim", "gmmat_data/phenotype1.gds")
seqBED2GDS("pheno_groups/phenotype2.bed", "pheno_groups/phenotype2.fam", "pheno_groups/phenotype2.bim", "gmmat_data/phenotype2.gds")
seqBED2GDS("pheno_groups/phenotype3.bed", "pheno_groups/phenotype3.fam", "pheno_groups/phenotype3.bim", "gmmat_data/phenotype3.gds")
seqBED2GDS("pheno_groups/phenotype4.bed", "pheno_groups/phenotype4.fam", "pheno_groups/phenotype4.bim", "gmmat_data/phenotype4.gds")
```

Now to adjust the covariate file for scripting and create separate snpsets for the slurm job array. 

```bash
module load r/3.5.2
for i in `seq 1 4`; do echo $i; perl -e 'chomp(@ARGV); %h; open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\s+/); $h{$s[0]} = $s[5];} close IN; open(IN, "< $ARGV[1]"); $t = <IN>; print "$t"; while(<IN>){chomp; @s = split(/\t/); if(exists($h{$s[0]})){$s[2] = $h{$s[0]}; print join("\t", @s); print "\n";}} close IN;' pheno_groups/phenotype${i}.fam pheno_cov_master_file.txt > gmmat_data/pheno${i}_covariates.tab; done

# Creating 1000 SNP subsets for processing
mkdir snp_sets
perl -e '$c = 0; $i = 1; open(OUT, "> snp_sets/snp_list_$i.txt"); while(<>){chomp; @s = split(/\t/); print OUT "$s[1]\n"; $c++; if($c >= 1000){close OUT; $i++; open(OUT, "> snp_sets/snp_list_$i.txt"); $c = 0;}} close OUT;' < pheno_groups/phenotype1.bim

# I wrote a job array script that processes each SNP set separately. Now to run it!
sbatch queue_gwas_job_array.sh gmmat_data/pheno1_covariates.tab /project/rumen_longread_metagenome_assembly/kiranmayee/IGC/prelim_gwas/round2/gemma/gemma_output/pheno1.cXX.txt gmmat_data/phenotype1.gds phenotype1

# Now to queue it all up to run overnight
for i in `seq 1 4`; do mkdir gmmat_output/phenotype${i}; done
mv gmmat_output/phenotype1*.tab gmmat_output/phenotype1/

# Damn! Kiranmayee mislabelled the covariate files because of her previous conventions!
cp /project/rumen_longread_metagenome_assembly/kiranmayee/IGC/prelim_gwas/round2/gemma/gemma_output/merge1.cXX.txt ./pheno4.cXX.txt
cp /project/rumen_longread_metagenome_assembly/kiranmayee/IGC/prelim_gwas/round2/gemma/gemma_output/pheno2.cXX.txt ./
cp /project/rumen_longread_metagenome_assembly/kiranmayee/IGC/prelim_gwas/round2/gemma/gemma_output/pheno3.cXX.txt ./

for i in `seq 2 4`; do echo $i; sbatch queue_gwas_job_array.sh gmmat_data/pheno${i}_covariates.tab pheno${i}.cXX.txt gmmat_data/phenotype${i}.gds phenotype${i}; done
```

And here are the scripts used to run the GWAS:

#### run_gmmat_gwas.R

```R
#!/software/7/apps/r/3.5.2/bin/Rscript
# A R script to automate GMMAT GWAS on subsets of SNPs
# ARGS = 1 : phenotype covariates, 2 : snp file, 3 : kinship file, 4 : input gds file, 5 : output filename

library(GMMAT)

args = commandArgs(trailingOnly=TRUE)

# test if there are at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one arguments must be supplied (snp ids file).\n", call.=FALSE)
}

pheno <- read.table(args[1], header=TRUE, stringsAsFactors=F)

snps <- scan(args[2], what="character")

kins <- as.matrix(read.table(args[3]))

colnames(kins)<-pheno$id

row.names(kins)<-pheno$id

results <- glmm.wald(fixed = bTB_status ~ breed + age + yr + season, data = pheno, id = "id", kins = kins, family = binomial(link = "logit"), infile = args[4], snps = snps)

write.table(results, file=args[5], row.names=F, sep = "\t", quote=FALSE)
```

#### queue_gwas_job_array.sh

```bash
#!/usr/bin/bash
# This is a script to queue up a job array of GWAS calculations on SNPs
# Input: $1 : Phenotype covariates file, $2 : kinship matrix file, $3 : input GDS file, $4 : output basename
#SBATCH --job-name=gmmatGWAS
#SBATCH --nodes=1
#SBATCH --time=3-0
#SBATCH --output=gmmat_stdout/gmmat_%A_%a.out
#SBATCH -p msn
#SBATCH --array=1-188
#SBATCH --mem=10000
#SBATCH --ntasks-per-node=1

N=$SLURM_ARRAY_TASK_ID

module load r/3.5.2
echo "Working on Slurm job array: $N"

time /software/7/apps/r/3.5.2/bin/Rscript --vanilla run_gmmat_gwas.R $1 snp_sets/snp_list_${N}.txt $2 $3 gmmat_output/${4}_${N}.tab

echo "Finished with Slurm task: $N"
```