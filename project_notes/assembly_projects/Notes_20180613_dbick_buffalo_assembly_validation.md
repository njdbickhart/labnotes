# Buffalo assembly statistics and QV generation
---
**6/13/2018**

These are my notes on the generation of QV and quality metric statistics on the long-read Buffalo reference genome.


## Table of Contents
* [Preparing the files](#preparing)
	* [Caspur assembly alignment](#caspur)
* [Generating quality metrics](#quality)
	* [calculate_qv.sh](#calculatescript)

<a name="preparing"></a>
## Preparing files

I am going to start downloading and preparing the files for alignment and processing. 

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/sequence_data/b_bubalus

```bash
# Downloading the next seq files first...
perl -lane 'if($F[0] =~ /Run/){next;}else{system("sbatch --nodes=1 --ntasks-per-node=2 --mem=2000 -p short --wrap=\"module load sratoolkit/gcc/64/2.8.2-1; prefetch --max-size 55000000 $F[0]; fastq-dump.2 -I --gzip --split-files $F[0]\"");}' < buffalo_wgs_sra_dataset.tab

# I downloaded the reference genome fasta, and I need to transfer it to the cluster for indexing and alignment
# OK, now to index the assembly and run the alignment
unpigz water_buffalo_20180219_gapf_noMito_arrowRename4_pilon.fasta.gz

sbatch --nodes=1 --ntasks-per-node=2 --mem=15000 -p short --wrap="module load bwa; module load samtools; bwa index water_buffalo_20180219_gapf_noMito_arrowRename4_pilon.fasta; samtools faidx water_buffalo_20180219_gapf_noMito_arrowRename4_pilon.fasta"

ls *.fastq.gz | perl -e 'use Cwd; %h; while(<>){chomp; @b = split(/_/); push(@{$h{$b[0]}}, $_);} $cwd = cwd(); foreach my $k (keys(%h)){print "$cwd/" . $h{$k}->[0] . "\t$cwd/" . $h{$k}->[1] . "\t$k\tbuffalo\n";}' > raw_sra_buffalo_paired_end.tab

# Now to fix the SRA files so that BWA doesn't throw a fit
perl -lane 'system("sbatch --nodes=1 --mem=6000 --ntasks-per-node=2 -p short --wrap=\"python3 ~/python_toolchain/sequenceData/fixSRAFastqFiles.py -f $F[0] -r $F[1] -o $F[2]\_r -l $F[2]\_r.log\"");' < raw_sra_buffalo_paired_end.tab

# Adding the reformated files to the spreadsheet
perl -lane '$F[0] =~ s/_1.fastq.gz/_r.1.fq.gz/; $F[1] =~ s/_2.fastq.gz/_r.2.fq.gz/; if( -s $F[0] && -s $F[1]){print join("\t", @F);}' < raw_sra_buffalo_paired_end.tab > reformated_sra_buffal_paired_end.tab

# Now to run my alignment pipeline
python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b buffalo -t reformated_sra_buffal_paired_end.tab -f water_buffalo_20180219_gapf_noMito_arrowRename4_pilon.fasta -p short -m
```

<a name="caspur"></a>
#### Caspur 2.0 WB assembly

It would be a good point of comparison to identify the quality metrics for this assembly vs the older caspur buffalo assembly.

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Bubalus_bubalis/CHR_Un/bbu_ref_UMD_CASPUR_WB_2.0_chrUn.fa.gz

unpigz bbu_ref_UMD_CASPUR_WB_2.0_chrUn.fa.gz

# I need to get rid of the horrible, horrible fasta header formatting. Let's do this in an interactive node
srun -N 1 -p short --ntasks-per-node=1 --mem=10000 --pty bash
perl -ne 'if($_ =~ /^>/){@bsegs = split(/[\|\s+]/); print ">$bsegs[3]\n";}else{print $_;}' < bbu_ref_UMD_CASPUR_WB_2.0_chrUn.fa > bbu_ref_UMD_CASPUR_WB_2.0_chrUn.rfmt.fa

module load samtools
samtools faidx bbu_ref_UMD_CASPUR_WB_2.0_chrUn.rfmt.fa

module load bwa/gcc/64/0.7.12
sbatch --nodes=1 --ntasks-per-node=2 --mem=18000 -p short --wrap="bwa index bbu_ref_UMD_CASPUR_WB_2.0_chrUn.rfmt.fa"

# And now to run my alignment pipeline
python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b caspur -t reformated_sra_buffal_paired_end.tab -f bbu_ref_UMD_CASPUR_WB_2.0_chrUn.rfmt.fa -p short -m
```

<a name="quality"></a>
## Generating quality metrics

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/sequence_data/b_bubalus

```bash
# First, let's start the longest running task (slurm jobid: 203460)
module load freebayes/1.1.0
sbatch --nodes=1 --ntasks-per-node=3 --mem=35000 -p medium --wrap="freebayes -C 2 -0 -O -q 20 -z 0.10 -E 0 -X -u -p 2 -F 0.75 -f water_buffalo_20180219_gapf_noMito_arrowRename4_pilon.fasta -v water_buffalo_20180219.freebayes.vcf buffalo/buffalo/buffalo.sorted.merged.bam"

# Now to generate SV calls (slurm jobid: 203461)
module load lumpy/gcc/64/0.2.13
sbatch --nodes=1 --mem=20000 --ntasks-per-node=2 -p short --wrap="lumpyexpress -B buffalo/buffalo/buffalo.sorted.merged.bam -o water_buffalo_20180219.lumpy.vcf"

# Hmm, I need to ask the scinet team to install FRC_align or code it myself without the lib-gmp dependencies

# QV generation
sbatch calculate_qv.sh water_buffalo_20180219.freebayes.vcf buffalo/buffalo/buffalo.sorted.merged.bam water_buffalo_20180219.qv

#### for the Caspur assembly ####
sbatch --nodes=1 --ntasks-per-node=3 --mem=35000 -p medium --dependency=afterany:203505 --wrap="freebayes -C 2 -0 -O -q 20 -z 0.10 -E 0 -X -u -p 2 -F 0.75 -f bbu_ref_UMD_CASPUR_WB_2.0_chrUn.rfmt.fa -v bbu_ref_UMD_CASPUR_WB_2.0.freebayes.vcf caspur/buffalo/buffalo.sorted.merged.bam"

sbatch --nodes=1 --mem=20000 --ntasks-per-node=2 -p short --dependency=afterany:203505 --wrap="lumpyexpress -B caspur/buffalo/buffalo.sorted.merged.bam -o bbu_ref_UMD_CASPUR_WB_2.0.lumpy.vcf"

# QV generation
sbatch calculate_qv.sh bbu_ref_UMD_CASPUR_WB_2.0.freebayes.vcf caspur/buffalo/buffalo.sorted.merged.bam bbu_ref_UMD_CASPUR_WB_2.0.qv
```


And here is the code for the QV calculation script.

<a name="calculatescript"></a>
#### calculate_qv.sh

```bash
#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=10000
#SBATCH --ntasks-per-node=1
#SBATCH -p short

# $1 = vcf file
# $2 = input bam file
# $3 = output QV file

module load samtools

NUM_BP=`samtools depth $2 | perl -e '$c = 0; while(<>){chomp; @s = split(/\t/); if(scalar(@s) >= 3){$c++;}} print "$c\n";'`
echo "num bp: "$NUM_BP

NUM_SNP=`cat $1 |grep -v "#" | awk -F "\t" '{if (!match($NF, "0/1")) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$8}' | tr ';' ' ' | sed s/AB=//g | awk -v WEIGHT=0 '{if ($6 >= WEIGHT) print $0}' | awk -v SUM=0 '{if (length($4) == length($5)) { SUM+=length($4); } else if (length($4) < length($5)) { SUM+=length($5)-length($4); } else { SUM+=length($4)-length($5)}} END { print SUM}'`
echo "num snp: "$NUM_SNP

perl -e 'chomp(@ARGV); $ns = $ARGV[0]; $nb = $ARGV[1]; print (-10 * log($ns/$nb)/log(10)); print "\n";' $NUM_SNP $NUM_BP > $3
cat $3
```