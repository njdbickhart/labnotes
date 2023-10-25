# Imputation pipeline running


## Table of Contents

## Setting up the files

OK, so I need to perform the following data manipulations in order to prepare the pipeline for use:

* Convert the mix files into PED, subset the SNPs by the common core and convert back to mix
* Modify the gVCF file sample names
* Merge the gVCF files into a singular file
* Ensure that the SNP markers correspond to the Sscro11.1 assembly coordinates
* Merge the pedigrees and potentially add the sequence sample IDs to the pedigree
* Modify the SNP manifest files to fit within the workflow
* Replicate all input files in the input directory

> Anunna: /scratch/HG/bickha01/side_projects/duroc_imputation/snp_genotypes

### Determining the core set of SNPs for the snp chip samples

```bash
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl manifest_1612_1.txt.list manifest_1612_3.txt.list manifest_1612_4.txt.list
File Number 1: manifest_1612_1.txt.list
File Number 2: manifest_1612_3.txt.list
File Number 3: manifest_1612_4.txt.list
Set     Count
1       719
1;2     49
1;2;3   30859
1;3     17771
2       15894
2;3     8676
3       5322

# I found out that only 100 samples were genotyped on chip 3, so here is the overlap on 1  4:

perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl manifest_1612_1.txt.list manifest_1612_4.txt.list       File Number 1: manifest_1612_1.txt.list
File Number 2: manifest_1612_4.txt.list
Set     Count
1       768
1;2     48630
2       13998

# Creating the shared file list for the snps
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -l 1_2 manifest_1612_1.txt.list manifest_1612_4.txt.list > shared_1_4.list

# Converting to ped file format temporarily
perl -lane '@gsegs = split(//, $F[2]); print "$F[1]\t$F[0]\t0\t0\t1\t3\t" . join("\t", @gsegs);' < genotype_1612_1.dat > genotype_1612_1.ped
perl -lane 'print "$F[3]\t$F[1]\t0\t$F[4]";' < manifest_1612_1.txt > genotype_1612_1.map
../bin/plink2 --file genotype_1612_1 --allow-extra-chr --extract shared_1_4.list --out genotype_1612_1_shared --make-bed


perl -lane '@gsegs = split(//, $F[2]); print "$F[1]\t$F[0]\t0\t0\t1\t3\t" . join("\t", @gsegs);' < genotype_1612_4.dat > genotype_1612_4.ped
perl -lane 'print "$F[3]\t$F[1]\t0\t$F[4]";' < manifest_1612_4.txt > genotype_1612_4.map
../bin/plink2 --file genotype_1612_4 --allow-extra-chr --extract shared_1_4.list --out genotype_1612_4_shared --make-bed

# merging the files
../bin/plink2 --bfile genotype_1612_1_shared --bmerge genotype_1612_4_shared --recode 12 --out genotype_sharedv --allow-extra-chr

# Converting to the gtp format
perl -lane '$name = $F[1]; @genos = (); for($i = 6; $i < scalar(@F); $i += 2){$t = $F[$i] + $F[$i+1] - 2; if($t < 0){$t = 5;} push(@genos, $t);} print "$name " . join("", @genos);'  < genotype_sharedv.ped > ../input/gtp.mix

# Here is the script with the code to convert back to GTP format:
# /lustre/backup/HG/bickha01/multispecies_chip/mixblup/modified_mixblup.sh

# creating the manifest file
perl -e '$c = 1; while(<STDIN>){chomp; @s = split(/\t/); print "$c\t$s[1]\t$s[0]\t$s[3]\n"; $c++;}' < genotype_sharedv.map > ../input/manifest.txt

# Wait, I need to trim more SNPs that have faulty chromosome IDs
perl -lane 'if($F[2] eq "0" || $F[2] =~ /^NW/){print "$F[1]";}' < input/manifest.txt > extra_markers_to_remove.list
sbatch -N 1 -n 2 --mem=15000 -p main -q std --wrap="../bin/plink2 --file genotype_sharedv --exclude extra_markers_to_remove.list --recode 12 --allow-extra-chr --out genotype_sharedv_filtered"

perl -lane '$name = $F[0]; @genos = (); for($i = 6; $i < scalar(@F); $i += 2){$t = $F[$i] + $F[$i+1] - 2; if($t < 0){$t = 5;} push(@genos, $t);} print "$name " . join("", @genos);'  < genotype_sharedv_filtered.ped > ../input/gtp.mix
perl -e '$c = 1; while(<STDIN>){chomp; @s = split(/\t/); if($s[0] == 24){$s[0] = 23;} print "$c\t$s[1]\t$s[0]\t$s[3]\n"; $c++;}' < genotype_sharedv_filtered.map > ../input/manifest.txt

cat pedigree_1612_1.txt pedigree_1612_4.txt | perl -e '%seen; while(<STDIN>){chomp; @s = split(/\s+/); if(exists($seen{$s[0]})){next;} if(scalar(@s) < 6){print "$s[0]   $s[1] $s[3] $s[0] $s[2]\n";}else{print "$s[0] $s[4] $s[5] $s[1] $s[3] $s[0] $s[2]\n";} $seen{$s[0]} = 1;}' > ../input/ped.mix

# Create the animal ID conversion file. Since they are the same, we will just double print the same codes. This also establishes the reference panel for the imputation
perl -lane 'print "$F[0]\t$F[0]";' < ../calls/freebayes/target_list_pedigree_locs.tab > ../input/codes.txt

# We need to list all of the samples in the VCF for them to be taken into account 
head -n 1276 input/input.vcf | tail -n 1 | perl -ne '@s = split(/\t/); for($x = 9; $x < scalar(@s); $x++){print "$s[$x]\n";}' > input/sample_all.txt

# We want to calculate accuracy, so we need to remove samples from the reference panel to impute. Since our number of samples with genotypes and sequence are 66 in total, we will remove the first third and calculate accuracy based on these values.
head -n 33 input/codes.txt | perl -lane 'print $F[0];' > input/target_list.txt
```

### Merging the VCF files

> Anunna: /scratch/HG/bickha01/side_projects/duroc_imputation/calls/freebayes

```bash
conda activate /lustre/backup/HG/bickha01/conda_envs/imputation/

ls *.gz > vcf_file.list

for i in *.gz; do echo $i; sbatch -N 1 -n 3 --mem=8000 -p main -q std -t 1-0 --wrap="bcftools index --threads 3 -f $i"; done

sbatch -t 2-0 merger_bcftools.sh vcf_file.list combined_prefiltered.vcf

# That did not work very well -- the non reference sites polluted the output VCF file! I need to try a different tool.
# Going to try glnexus using the deep variant gvcf style
sbatch -t 2-0 merger_glnexus.sh vcf_file.list combined_prefiltered.bcf /scratch/HG/bickha01/side_projects/duroc_imputation/sscrofa11_1.fasta.bed

# NOPE! At least glnexus works, so I need to see if I can use something like DeepVariant on my future workflows

# Let's strip the gvcf info and combine the files. X needs to be converted to 23, and all sample names need to be stripped of - and _ characters
for i in *.gvcf.gz; do name=`echo $i | cut -d'.' -f1 | perl -ne 'chomp; $_ =~ s/[\-\_]//g; print $_;'`; echo $name; sbatch -N 1 -n 2 --mem=8000 -p main -q std --wrap="python ~/python_toolchain/sequenceData/freebayesGVCFtoVCF.py -f $i -c X,23 -s $name -o $name.vcf; bgzip $name.vcf; bcftools index $name.vcf.gz"; done

ls *.vcf.gz > real_vcf_file.list
sbatch -t 2-0 merger_bcftools.sh real_vcf_file.list combined_prefiltered.vcf
sbatch --dependency=afterok:47810181 -t 2-0 bcftools_filter.sh combined_prefiltered.vcf combined_anremoved.vcf.gz

# That took too long and was not able to finish. Let's separate it by chromosome
for i in `seq 1 18` 23; do echo $i; sbatch -t 2-0 merger_bcftools_bychr.sh real_vcf_file.list combined.$i.prefiltered.vcf $i; done
# filtering the vcfs
for i in combined.*.vcf; do name=`echo $i | cut -d'.' -f1,2`; echo $name; sbatch -t 2-0 bcftools_filter.sh $i $name.filtered.vcf; done

# Combining it all
sbatch -N 1 -n 25 --mem=15000 -p main -q std -t 1-0 --wrap="bcftools concat --threads 25 -O v -o combined.all.final.vcf combined.1.filtered.vcf combined.2.filtered.vcf combined.3.filtered.vcf combined.4.filtered.vcf combined.5.filtered.vcf combined.6.filtered.vcf combined.7.filtered.vcf combined.8.filtered.vcf combined.9.filtered.vcf combined.10.filtered.vcf combined.11.filtered.vcf combined.12.filtered.vcf combined.13.filtered.vcf combined.14.filtered.vcf combined.15.filtered.vcf combined.16.filtered.vcf combined.17.filtered.vcf combined.18.filtered.vcf combined.23.filtered.vcf"

sbatch -N 1 -n 25 --mem=15000 -p main -q std -t 1-0 --wrap='bcftools stats combined.all.final.vcf -s - > combined.all.final.vchk; plot-vcfstats -p combined_all_stats -v combined.all.final.vchk'

# Pulling the sample IDs that match genotypes for the accuracy estimation
for i in *.vcf.gz; do name=`echo $i | cut -d'.' -f1`; grep $name ../../snp_genotypes/pedigree*.txt.list; done > target_list.txt
```

## Running the pipeline

> Anunna: /scratch/HG/bickha01/side_projects/duroc_imputation

```bash
conda activate /lustre/backup/HG/bickha01/conda_envs/imputation/

# Addie has most of the workflow setup under a series of scripts. I will start with the initial script and see where the process takes us!
# The parameters follow the example setup by Addie in his workflow document
sbatch -o imputation_pipeline.log -t 3-0 ./bin/imputation_pipeline.sh 2 2 10 1 50 10 10 8
```