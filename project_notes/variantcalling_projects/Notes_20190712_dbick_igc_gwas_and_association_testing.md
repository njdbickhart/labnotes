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

# Pheno groups 3 and 4 failed. I need to calculate the kinship matrix for these datasets. Let's do that now
module load gemma/0.96

# OK, so the goofballs managing the packages in CERES didn't finish the installation of gemma. I need to do a workaround
cp /software/7/apps/gemma/0.96/gemma.linux ./
chmod +x gemma.linux
./gemma.linux -bfile pheno_groups/phenotype3 -gk 1 -o pheno3.correct
./gemma.linux -bfile pheno_groups/phenotype4 -gk 1 -o pheno4.correct

# This goes by fast, so lets queue up the rest
./gemma.linux -bfile pheno_groups/phenotype2 -gk 1 -o pheno2.correct
./gemma.linux -bfile pheno_groups/phenotype1 -gk 1 -o pheno1.correct

# Now, let's resubmit the scripts and try again
for i in 1 2 3 4; do echo $i; sbatch queue_gwas_job_array.sh gmmat_data/pheno${i}_covariates.tab output/pheno${i}.correct.cXX.txt gmmat_data/phenotype${i}.gds phenotype${i}; done

# Damn Ceres storage issues! I just rewrote the script to avoid regenerating already existing files
for i in 1 2 3 4; do echo $i; sbatch queue_gwas_job_array.sh gmmat_data/pheno${i}_covariates.tab output/pheno${i}.correct.cXX.txt gmmat_data/phenotype${i}.gds phenotype${i}; done
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

if [ -f gmmat_output/${4}/${4}_${N}.tab ]; then
        echo "Skipping gmmat_output/${4}/${4}_${N}.tab as it already exists!"
else
        time /software/7/apps/r/3.5.2/bin/Rscript --vanilla run_gmmat_gwas.R $1 snp_sets/snp_list_${N}.txt $2 $3 gmmat_output/${4}/${4}_${N}.tab
fi

echo "Finished with Slurm task: $N"
```

## Printing the manhattan plots

OK, now I have replicated the GWAS but with a different set of criteria SNPs. Now to produce the plots for the data and check them all out.

Kiranmayee used an R vector to store and highlight the SNPs of interest and I will adopt her pipeline to save time. 

> Ceres: /project/rumen_longread_metagenome_assembly/kiranmayee/IGC/case_control

```bash
# Merging the data in order
for i in gmmat_output/phenotype*; do echo $i; perl -e 'chomp(@ARGV); @f = `ls $ARGV[0]/*.tab`; @f = sort{@s = split(/[_\.]/, $a); @j = split(/[_\.]/, $b); $s[-2] <=> $j[-2]}@f; chomp(@f); open(OUT, "> $ARGV[1]"); print {OUT} "SNP\tCHR\tPOS\tREF\tALT\tN\tAF\tBETA\tSE\tPVAL\tconverged\n"; foreach $k (@f){open(IN, "< $k"); <IN>; while($l = <IN>){print {OUT} $l;} close IN;}' $i ${i}_merged.tab; done

# Now to copy and modify Kiranmayees scripts
cp /project/rumen_longread_metagenome_assembly/kiranmayee/IGC/prelim_gwas/round2/gemma/plot.sh ./manhattan_plot.sh

# SNPs of interest to add to the plot
perl -ne 'chomp; @s = split(/\s+/); print "\"$s[1]\", ";' < neogen_merged_nonhet_filt.revised.bim ; echo

# Now to generate the plots
sbatch manhattan_plot.sh

# Not satisfying, but I can work with this. I'm also going to pull the effect sizes and the standard errors for the SNPs of interest:
for i in `seq 1 4`; do echo $i; python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f gmmat_output/phenotype${i}_merged.tab -c 0 -l new_filt.list -d '\t' > gmmat_output/phenotype${i}_onlycustom.tab; done

grep '18_63417698' neogen_merged_hardy.hwe
  18            18_63417698      ALL    A    T          74/503/1101   0.2998   0.3127       0.1007
  18            18_63417698      AFF    A    T           52/338/847   0.2732   0.2935      0.01961
  18            18_63417698    UNAFF    A    T           22/165/254   0.3741   0.3616       0.5126
grep '5_99190989' neogen_merged_hardy.hwe
   5             5_99190989      ALL    T    A           5/272/1324   0.1699   0.1606      0.01865
   5             5_99190989      AFF    T    A           2/170/1001   0.1449   0.1373      0.05646
   5             5_99190989    UNAFF    T    A            3/102/323   0.2383   0.2205       0.1231

grep '18_63417698' gmmat_output/*onlycustom.tab
gmmat_output/phenotype1_onlycustom.tab:18_63417698      18      63640110        A       T       931     0.804511278195489       0.4852857470566790.12678233520238 0.000129342118284114    TRUE
gmmat_output/phenotype2_onlycustom.tab:18_63417698      18      63640110        A       T       1434    0.808228730822873       0.4314844652376850.107676779156407        6.14378180734818e-05    TRUE
gmmat_output/phenotype3_onlycustom.tab:18_63417698      18      63640110        A       T       1468    0.806198910081744       0.4262955850530670.106223275965894        5.9900320331582e-05     TRUE
gmmat_output/phenotype4_onlycustom.tab:18_63417698      18      63640110        A       T       1672    0.806519138755981       0.4052558522973330.105021658954359        0.00011395279421179     TRUE

grep '5_99190989' gmmat_output/*onlycustom.tab
gmmat_output/phenotype1_onlycustom.tab:5_99190989       5       99676855        T       A       886     0.899548532731377       0.5943330751059990.180746260762023        0.00100823800625793     TRUE
gmmat_output/phenotype2_onlycustom.tab:5_99190989       5       99676855        T       A       1367    0.909656181419166       0.6349225238747220.155864545342061        4.63012607020109e-05    TRUE
gmmat_output/phenotype3_onlycustom.tab:5_99190989       5       99676855        T       A       1398    0.910228898426323       0.6610244380176640.154400077159731        1.85851256334299e-05    TRUE
gmmat_output/phenotype4_onlycustom.tab:5_99190989       5       99676855        T       A       1595    0.912225705329154       0.6181864161234680.151326691057893        4.40556430047843e-05    TRUE
```

All beta values were relatively small, suggesting that the effects are very tiny for each variant on the regression. This is reflected by the hodge-podge of genotype states among affected and unaffected individuals. Unfortunately, there isn't much of a story here, but we can make the case that this is new and will have impacts downstream. 

## Calculating heritability

I am going to attempt to do this using the SOMMER package

```R
library(sommer)

```

## Generating final table

```bash
perl -e 'chomp(@ARGV); %d; %m; foreach $n (@ARGV){open(IN, "< $n"); while(<IN>){chomp; @s = split(/\t/); $m{$s[0]} = [$s[1], $s[2], $s[3], $s[4]]; $d{$s[0]}->{$n} = [sprintf("%0.3f", $s[7]), sprintf("%0.3f",$s[9])];} close IN;} foreach $snp (keys(%m)){ print "$snp\t" . join("\t", @{$m{$snp}}); foreach $f (@ARGV){print "\t"; print join("\t", @{$d{$snp}->{$f}});} print "\n";}' gmmat_output/phenotype1_onlycustom.tab gmmat_output/phenotype2_onlycustom.tab gmmat_output/phenotype3_onlycustom.tab gmmat_output/phenotype4_onlycustom.tab
MHC_86084       23      28405622        A       G       -0.001  0.994   -0.014  0.878   -0.025  0.778   -0.013  0.881
18_62670367     18      62953374        C       G       0.034   0.748   0.006   0.951   0.014   0.879   0.006   0.952
MHC_154399      23      28449842        C       A       -0.056  0.562   -0.082  0.336   -0.080  0.344   -0.087  0.299
MHC_143922      23      28443060        T       A       -0.181  0.281   -0.170  0.252   -0.165  0.259   -0.149  0.310
MHC_59013       23      28388100        G       A       -0.188  0.188   -0.164  0.210   -0.175  0.176   -0.181  0.161
23_28651088     23      28586018        G       A       0.072   0.472   0.030   0.743   0.029   0.745   0.040   0.656
KIR_41480       18      62973689        G       T       -0.111  0.449   -0.034  0.793   -0.036  0.778   -0.074  0.559
5_99203402      5       99689138        C       T       -0.028  0.782   0.106   0.251   0.093   0.308   0.076   0.398
18_63122963     18      63369272        T       C       -0.044  0.596   -0.039  0.588   -0.036  0.619   -0.051  0.480
18_63082203     18      63331817        A       C       -0.102  0.412   -0.109  0.320   -0.107  0.325   -0.115  0.288
18_63084493     18      63333922        G       C       -0.233  0.262   -0.285  0.126   -0.274  0.137   -0.203  0.275
18_63417698     18      63640110        A       T       0.485   0.000   0.431   0.000   0.426   0.000   0.405   0.000
18_62559417     18      62851420        G       A       0.092   0.428   0.072   0.479   0.076   0.449   0.114   0.254
5_99412564      5       99896106        G       A       -0.290  0.095   -0.121  0.447   -0.130  0.411   -0.078  0.619
18_62994372     18      63251108        C       T       -0.029  0.804   -0.012  0.910   0.002   0.982   -0.014  0.890
18_62716825     18      62996065        T       G       -0.077  0.528   -0.054  0.619   -0.053  0.619   -0.051  0.637
LRC_98785       18      63202343        G       A       0.072   0.499   0.155   0.104   0.161   0.087   0.144   0.124
LRC_106729      18      63202745        C       G       0.139   0.156   0.137   0.113   0.128   0.137   0.105   0.217
LIB14427_MHC_57505      23      28512858        T       C       -0.096  0.608   -0.126  0.448   -0.131  0.420   -0.123  0.451
MHC_286137      23      28535115        T       C       0.205   0.263   0.186   0.246   0.215   0.174   0.153   0.325
LIB14427_MHC_3538       23      28445610        A       G       0.210   0.039   0.145   0.094   0.161   0.060   0.161   0.059
18_63131111     18      63376760        C       T       -0.052  0.701   -0.008  0.944   -0.015  0.897   -0.020  0.867
LIB14427_MHC_2500       23      28444317        A       T       -0.048  0.617   -0.061  0.468   -0.056  0.499   -0.051  0.536
18_63340289     18      63568977        C       T       -0.134  0.202   -0.050  0.583   -0.052  0.570   -0.054  0.546
18_63156196     18      63399811        A       C       0.041   0.708   -0.005  0.961   -0.006  0.952   0.002   0.987
5_99780983      5       100260661       A       C       0.058   0.722   0.112   0.439   0.111   0.441   0.134   0.353
23_28648192     23      28583505        G       T       -0.028  0.789   -0.051  0.588   -0.062  0.512   -0.036  0.699
LRC_81028       18      63201446        T       C       -0.047  0.582   -0.081  0.275   -0.079  0.285   -0.061  0.403
MHC_260913      23      28518788        G       A       -0.201  0.014   -0.183  0.012   -0.176  0.014   -0.174  0.014
5_99036250      5       99523739        C       T       0.208   0.089   0.205   0.052   0.209   0.045   0.183   0.079
18_62460291     18      62760331        C       T       -0.183  0.105   -0.183  0.066   -0.173  0.078   -0.147  0.133
LIB14427_MHC_45690      23      28498135        T       G       -0.046  0.639   -0.036  0.679   -0.028  0.742   -0.035  0.677
5_99714227      5       100194605       G       A       -0.039  0.885   0.189   0.438   0.168   0.487   0.197   0.409
MHC_208321      23      28484745        C       T       -0.001  0.995   0.013   0.913   0.022   0.850   0.015   0.896
18_63186702     18      63427843        C       T       -0.150  0.394   -0.122  0.435   -0.119  0.440   -0.127  0.405
LIB14427_MHC_36271      23      28486398        T       C       0.093   0.435   0.030   0.778   0.034   0.741   0.018   0.864
18_63399311     18      63623214        G       C       0.033   0.761   0.091   0.346   0.088   0.354   0.065   0.496
LRC_174904      18      63206188        A       G       0.005   0.970   0.053   0.640   0.068   0.545   0.032   0.774
MHC_43656       23      28378159        G       T       -0.058  0.615   -0.102  0.318   -0.114  0.265   -0.101  0.316
MHCclassI_MHC_392940    23      28604248        C       A       0.346   0.251   0.263   0.313   0.299   0.246   0.274   0.286
18_62774779     18      63049320        A       C       0.158   0.060   0.160   0.029   0.146   0.043   0.133   0.064
LIB14427_MHC_73766      23      28533120        A       G       0.103   0.402   0.047   0.666   0.056   0.603   0.040   0.712
5_99333595      5       99817965        A       C       0.067   0.575   0.094   0.368   0.076   0.462   0.066   0.522
18_63067935     18      63318706        A       C       0.083   0.442   0.148   0.126   0.156   0.103   0.137   0.149
18_63036451     18      63289775        A       G       0.057   0.570   0.118   0.184   0.126   0.153   0.116   0.187
MHCclassI_MHC_302664    23      28545813        A       C       0.045   0.674   0.038   0.690   0.042   0.652   0.016   0.864
18_63089154     18      63338205        C       T       -0.022  0.856   0.011   0.920   0.012   0.911   0.031   0.766
MHC_324231      23      28559773        C       T       -0.161  0.259   -0.160  0.209   -0.161  0.202   -0.121  0.342
5_99190989      5       99676855        T       A       0.594   0.001   0.635   0.000   0.661   0.000   0.618   0.000
MHC_115082      23      28424393        G       A       -0.120  0.246   -0.084  0.367   -0.092  0.324   -0.081  0.376
LIB14427_MHC_126886     23      28599312        A       G       0.153   0.262   0.135   0.262   0.159   0.185   0.168   0.157
5_99094858      5       99581732        A       G       -0.007  0.944   0.136   0.143   0.129   0.162   0.115   0.204
5_99437819      5       99921096        G       A       0.047   0.643   0.149   0.094   0.146   0.097   0.154   0.077
18_62644240     18      62929365        C       T       -0.069  0.622   -0.073  0.564   -0.051  0.680   -0.079  0.523
MHC_380177      23      28595986        T       C       -0.000  0.998   0.065   0.410   0.053   0.500   0.041   0.597
LRC_72198       18      63201000        G       C       0.108   0.292   0.087   0.340   0.071   0.435   0.102   0.260
MHCclassI_MHC_120784    23      28428083        T       A       -0.078  0.382   -0.085  0.284   -0.073  0.354   -0.081  0.303
18_62665698     18      62949083        G       A       0.102   0.395   0.111   0.297   0.114   0.279   0.110   0.294
MHC_400205      23      28608950        G       A       0.147   0.207   0.058   0.573   0.064   0.529   0.083   0.416
5_99706728      5       100187185       A       G       0.165   0.583   0.326   0.223   0.325   0.223   0.353   0.180
MHCclassI_MHC_395846    23      28606129        A       G       0.361   0.012   0.220   0.081   0.219   0.079   0.243   0.052
18_63385249     18      63610292        T       C       0.016   0.902   0.035   0.772   0.042   0.719   0.020   0.868
LRC_61352       18      63200453        G       A       0.080   0.643   0.026   0.861   0.014   0.926   0.108   0.471
18_63185710     18      63426932        A       G       -0.109  0.314   -0.058  0.550   -0.045  0.632   -0.059  0.535
LIB14427_MHC_122678     23      28594068        G       A       0.324   0.021   0.189   0.130   0.200   0.105   0.214   0.085
MHCclassI_MHC_181938    23      28467668        C       T       -0.082  0.636   -0.141  0.360   -0.126  0.405   -0.163  0.286
LIB14427_MHC_116810     23      28586756        G       A       -0.010  0.916   -0.038  0.647   -0.027  0.735   -0.036  0.651
5_99445508      5       99928705        T       A       0.064   0.598   0.029   0.780   0.025   0.814   0.055   0.598
5_99763326      5       100243190       C       G       0.196   0.414   0.079   0.700   0.084   0.680   0.096   0.635
18_62572950     18      62863855        C       T       -0.054  0.622   0.006   0.952   0.015   0.873   0.004   0.962
18_63114542     18      63361534        A       G       -0.011  0.957   -0.069  0.687   -0.091  0.590   -0.129  0.437
MHC_9213        23      28355864        A       G       -0.176  0.259   -0.147  0.299   -0.164  0.244   -0.153  0.271
```