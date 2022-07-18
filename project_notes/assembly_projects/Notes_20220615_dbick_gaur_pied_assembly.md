# Gaur Pied analysis 
---
*9/16/2020*

These are my notes on quality control and publication of the Gaur x Pied manuscript

## Table of Contents


## Data locations

* Original assembly: /project/gaur_genome_assembly/Gaur_x_Pied/Guppy/GxP_HiFi_17cells_binning/Pied_dam.contigs.pd2-all.salsa.edits.tidy.MT.fasta
* Illumina data folder: /project/gaur_genome_assembly/Gaur_x_Pied/Guppy/illumina_data/

## Themis-ASM run on Pied

> Ceres: /project/gaur_genome_assembly/Gaur_x_Pied/Guppy/themis_run

```bash
module load bedtools

# assessing repetitive content before running Themis
perl -e '<>; <>; <>; while(<>){chomp; $_ =~ s/^\s+//; @s = split(/\s+/, $_); print "$s[4]\t$s[5]\t$s[6]\t$s[9]\n";}' < ../GxP_HiFi_17cells_binning/RepeatMakser/Pied_dam.contigs.pd2-all.salsa.edits.tidy.fasta.out > pied_dam_repeats.bed

perl -lane 'print "$F[0]\t1\t$F[1]";' < ../GxP_HiFi_17cells_binning/Pied_dam.contigs.pd2-all.salsa.edits.tidy.MT.fasta.fai > pied_dam_chrlens.bed

bedtools merge -i pied_dam_repeats.bed | bedtools subtract -a pied_dam_chrlens.bed -b stdin > pied_dam_nonrepetitive.bed

perl -e '%h; while(<>){chomp; @s = split(/\t/); $h{$s[0]} += $s[2] - $s[1];} foreach $k (keys(%h)){print "$k\t$h{$k}\n";}' < pied_dam_nonrepetitive.bed > pied_dam_nonrepetitive.tab

perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %h; while(<IN>){chomp; @s = split(/\t/); $h{$s[0]}->[0] = $s[1];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $h{$s[0]}->[1] = $s[1];} close IN; foreach $k (keys(%h)){$t = $h{$k}->[0] / $h{$k}->[1]; print "$k\t" . $h{$k}->[1] . "\t$t\n";}' pied_dam_nonrepetitive.tab ../GxP_HiFi_17cells_binning/Pied_dam.contigs.pd2-all.salsa.edits.tidy.MT.fasta.fai | perl -lane 'if($F[2] <= 0.20){print $_;}' > highly_repetitive_contigs.tab
module load r/4.1.2
```

```R
library(ggplot2)
data <- read.delim("highly_repetitive_contigs.tab", header=FALSE, sep="\t")
colnames(data) <- c("Scaffold", "Length", "NonRep")

pdf(file="highly_repetitive_contigs.pdf", useDingbats=FALSE)
ggplot(data, aes(Length, NonRep, color=NonRep < 0.10)) + geom_point() + geom_hline(yintercept=0.10, linetype="dashed", color="red") + geom_vline(xintercept=100000, linetype="dashed", color="blue") + scale_x_log10()
dev.off()
```

Now to screen contigs that are less than 200 kb in length and have > 90% repetitive content. Then I'll run Snakemake and start the Themis run.

```bash
module load miniconda snakemake/6.12.1 java/1.8.0_121 samtools

perl -lane 'if ($F[1] <= 200000 && $F[2] < 0.1){print $F[0];}' < highly_repetitive_contigs.tab > contigs_to_filter.txt

python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f ../GxP_HiFi_17cells_binning/Pied_dam.contigs.pd2-all.salsa.edits.tidy.MT.fasta.fai -c 0 -d '\t' -l contigs_to_filter.txt -v | perl -lane 'print "$F[0]\t1\t$F[1]\t$F[0]\t1\t+";' > contigs_to_keep.bed

sbatch -N 1 -n 5 --mem=30000 -p memlimit -q short -t 1-0 --wrap="java -Xmx35g -jar ~/rumen_longread_metagenome_assembly/binaries/CombineFasta/store/CombineFasta.jar agp2fasta -f ../GxP_HiFi_17cells_binning/Pied_dam.contigs.pd2-all.salsa.edits.tidy.MT.fasta -b contigs_to_keep.bed -o Pied_dam.contigs.filtered.fasta"

illumina_folder=/project/gaur_genome_assembly/Gaur_x_Pied/Guppy/illumina_data


PATH=$PATH:/software/7/apps/merqury/1.1/:/software/7/apps/meryl/1.0/Linux-amd64/bin/
MERQURY=/software/7/apps/merqury/1.1
PERL5LIB=''
sbatch -N 1 -n 2 --mem=6000 -p priority -q msn --wrap='python3 ~/rumen_longread_metagenome_assembly/binaries/Themis-ASM/themisASM.py -a /project/gaur_genome_assembly/Gaur_x_Pied/Guppy/themis_run/Pied_dam.contigs.filtered.fasta -n piedmontese -a /project/gaur_genome_assembly/Gaur_x_Pied/Guppy/themis_run/ARS_UCD1.3_corrected.fasta -n arsucd13 -b mammalia_odb10 -f /project/gaur_genome_assembly/Gaur_x_Pied/Guppy/illumina_data/LIB105190_S1_L001_R1_001.fastq.gz,/project/gaur_genome_assembly/Gaur_x_Pied/Guppy/illumina_data/LIB105190_S1_L001_R2_001.fastq.gz -s lane1 -f /project/gaur_genome_assembly/Gaur_x_Pied/Guppy/illumina_data/LIB105190_S1_L002_R1_001.fastq.gz,/project/gaur_genome_assembly/Gaur_x_Pied/Guppy/illumina_data/LIB105190_S1_L002_R2_001.fastq.gz -s lane2 -f /project/gaur_genome_assembly/Gaur_x_Pied/Guppy/illumina_data/LIB105190_S1_L003_R1_001.fastq.gz,/project/gaur_genome_assembly/Gaur_x_Pied/Guppy/illumina_data/LIB105190_S1_L003_R2_001.fastq.gz -s lane3 -f /project/gaur_genome_assembly/Gaur_x_Pied/Guppy/illumina_data/LIB105190_S1_L004_R1_001.fastq.gz,/project/gaur_genome_assembly/Gaur_x_Pied/Guppy/illumina_data/LIB105190_S1_L004_R2_001.fastq.gz -s lane4 -c "sbatch --nodes={cluster.nodes} --ntasks-per-node={cluster.ntasks-per-node} --mem={cluster.mem} --partition={cluster.partition} -q {cluster.qos} -o {cluster.stdout}" -j 100 --resume'

# I used the F1 reads here. Let's try to collect the data on the DAM side instead
mv themis_summary.zip f1_reads_themis_summary.zip

ls  /project/gaur_genome_assembly/Gaur_x_Pied/dam/run?/*.gz | perl -e 'use File::Basename; while($f = <>){$r = <>; chomp($f, $r); $fn = basename($f); $rn = basename($r); @fsegs = split(/_/, $fn); @rsegs = split(/_/, $rn); @fpsegs = split(/\//, $f); @rpsegs = split(/\//, $r); $sname = "$fpsegs[5]\_$fsegs[2]"; print "-s $sname -f $f,$r ";}'

sbatch -N 1 -n 2 --mem=6000 -p priority -q msn --wrap='python3 ~/rumen_longread_metagenome_assembly/binaries/Themis-ASM/themisASM.py -a /project/gaur_genome_assembly/Gaur_x_Pied/Guppy/themis_run/Pied_dam.contigs.filtered.fasta -n piedmontese -a /project/gaur_genome_assembly/Gaur_x_Pied/Guppy/themis_run/ARS_UCD1.3_corrected.fasta -n arsucd13 -b mammalia_odb10 -s run1_L001 -f /project/gaur_genome_assembly/Gaur_x_Pied/dam/run1/LIB105188_S1_L001_R1_001.fastq.gz,/project/gaur_genome_assembly/Gaur_x_Pied/dam/run1/LIB105188_S1_L001_R2_001.fastq.gz -s run1_L002 -f /project/gaur_genome_assembly/Gaur_x_Pied/dam/run1/LIB105188_S1_L002_R1_001.fastq.gz,/project/gaur_genome_assembly/Gaur_x_Pied/dam/run1/LIB105188_S1_L002_R2_001.fastq.gz -s run1_L003 -f /project/gaur_genome_assembly/Gaur_x_Pied/dam/run1/LIB105188_S1_L003_R1_001.fastq.gz,/project/gaur_genome_assembly/Gaur_x_Pied/dam/run1/LIB105188_S1_L003_R2_001.fastq.gz -s run1_L004 -f /project/gaur_genome_assembly/Gaur_x_Pied/dam/run1/LIB105188_S1_L004_R1_001.fastq.gz,/project/gaur_genome_assembly/Gaur_x_Pied/dam/run1/LIB105188_S1_L004_R2_001.fastq.gz -s run2_L001 -f /project/gaur_genome_assembly/Gaur_x_Pied/dam/run2/LIB105188_S1_L001_R1_001.fastq.gz,/project/gaur_genome_assembly/Gaur_x_Pied/dam/run2/LIB105188_S1_L001_R2_001.fastq.gz -s run2_L002 -f /project/gaur_genome_assembly/Gaur_x_Pied/dam/run2/LIB105188_S1_L002_R1_001.fastq.gz,/project/gaur_genome_assembly/Gaur_x_Pied/dam/run2/LIB105188_S1_L002_R2_001.fastq.gz -s run2_L003 -f /project/gaur_genome_assembly/Gaur_x_Pied/dam/run2/LIB105188_S1_L003_R1_001.fastq.gz,/project/gaur_genome_assembly/Gaur_x_Pied/dam/run2/LIB105188_S1_L003_R2_001.fastq.gz -s run2_L004 -f /project/gaur_genome_assembly/Gaur_x_Pied/dam/run2/LIB105188_S1_L004_R1_001.fastq.gz,/project/gaur_genome_assembly/Gaur_x_Pied/dam/run2/LIB105188_S1_L004_R2_001.fastq.gz -c "sbatch --nodes={cluster.nodes} --ntasks-per-node={cluster.ntasks-per-node} --mem={cluster.mem} --partition={cluster.partition} -q {cluster.qos} -o {cluster.stdout}" -j 100 --resume'
```

> Ceres: /project/gaur_genome_assembly/Gaur_x_Pied/Guppy/dominette_themis

```bash
module load sratoolkit/2.10.9

fastq-dump --split-files SRR5753600

mkdir logs
module load miniconda snakemake/6.12.1 java/1.8.0_121 samtools
PATH=$PATH:/software/7/apps/merqury/1.1/:/software/7/apps/meryl/1.0/Linux-amd64/bin/
MERQURY=/software/7/apps/merqury/1.1
PERL5LIB=''

sbatch -N 1 -n 2 --mem=6000 -p priority -q msn --wrap='python3 ~/rumen_longread_metagenome_assembly/binaries/Themis-ASM/themisASM.py -a /project/gaur_genome_assembly/Gaur_x_Pied/Guppy/themis_run/Pied_dam.contigs.filtered.fasta -n piedmontese -a /project/gaur_genome_assembly/Gaur_x_Pied/Guppy/themis_run/ARS_UCD1.3_corrected.fasta -n arsucd13 -b mammalia_odb10 -f /project/gaur_genome_assembly/Gaur_x_Pied/Guppy/dominette_themis/SRR5753600_1.fastq,/project/gaur_genome_assembly/Gaur_x_Pied/Guppy/dominette_themis/SRR5753600_2.fastq -s dominette -c "sbatch --nodes={cluster.nodes} --ntasks-per-node={cluster.ntasks-per-node} --mem={cluster.mem} --partition={cluster.partition} -q {cluster.qos} -o {cluster.stdout}" -j 100 --resume'
```