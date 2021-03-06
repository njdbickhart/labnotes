# Protist genome assembly project
---
*10/4/2019*

These are my notes on how to sequence and assemble the rumen protist population.

## Table of Contents

## DNA isolation and extraction troubleshooting

So, we've done some extractions of the rumen protists and sequenced them on the GridION, but the read lengths are short. I want to see if we're just assembling the bacterial fragmentary proportion or if we're getting close to the real data. There are a few tests I'd like to run to check this. The first would be the use of GC data partitioning and the second would be kmer-based genome size estimates.

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover

```bash
module load miniconda jellyfish2/2.2.9
# First, let's concatenate the fastqs into separate pools 

for i in Proto_T3_Q2 Proto_T3_R2 Proto_T3_WL2 RCL_Proto1 RC_Proto1; do echo $i; cat $i/fastq_pass/pass/*.fastq > $i.concatenated.fastq; done

# Now to generate Jellyfish dbs
for i in Proto_T3_Q2 Proto_T3_R2 Proto_T3_WL2 RCL_Proto1 RC_Proto1; do echo $i; sbatch --nodes=1 --mem=100000 --ntasks-per-node=5 -p short --wrap="jellyfish count -m 21 -s 9G -t 5 -o ${i}.jf ${i}.concatenated.fastq"; done

source activate /KEEP/rumen_longread_metagenome_assembly/kat
for i in Proto_T3_Q2 Proto_T3_R2 Proto_T3_WL2 RCL_Proto1 RC_Proto1; do echo $i; sbatch --nodes=1 --mem=50000 --ntasks-per-node=15 -p short --wrap="source activate /KEEP/rumen_longread_metagenome_assembly/kat; kat hist -t 15 -o ${i}_kat_hist ${i}.jf"; done

# I rewrote my fastq python script to process fastqs individually
source activate /KEEP/rumen_longread_metagenome_assembly/seaborn

python3 ~/python_toolchain/sequenceData/fastqStats.py -f Proto_T3_Q2.concatenated.fastq -o Proto_T3_Q2
for i in Proto_T3_R2 Proto_T3_WL2 RCL_Proto1 RC_Proto1; do echo $i; python3 ~/python_toolchain/sequenceData/fastqStats.py -f ${i}.concatenated.fastq -o $i; done

## Adding new data from the new preps Jen performed
for i in Proto_QiaT1 B_Pip_Proto Proto_HMW1ug Proto_PCI_HMW_DNA Proto_PCI_Spun_DNA Proto_HMW2ug; do echo $i; ls $i/*/*pass/*.fastq | wc -l; done
Proto_QiaT1
472
B_Pip_Proto
1044
Proto_HMW1ug
2156
Proto_PCI_HMW_DNA
3388
Proto_PCI_Spun_DNA
1928
Proto_HMW2ug
2008

for i in Proto_QiaT1 B_Pip_Proto Proto_HMW1ug Proto_PCI_HMW_DNA Proto_PCI_Spun_DNA Proto_HMW2ug; do echo $i; cat $i/*/*pass/*.fastq > $i.concatenated.fastq; done
for i in Proto_QiaT1 B_Pip_Proto Proto_HMW1ug Proto_PCI_HMW_DNA Proto_PCI_Spun_DNA Proto_HMW2ug; do echo $i; sbatch --nodes=1 --mem=18000 --ntasks-per-node=2 -p msn -q msn --wrap="python3 ~/python_toolchain/sequenceData/fastqStats.py -f ${i}.concatenated.fastq -o $i"; done

python3 ~/python_toolchain/sequenceData/fastqStats.py -f New_cell_proto.concatenated.fastq -o New_cell_proto
```

## Preliminary assembly

We have enough data to start trying a "diagnostic" assembly of the reads. I want to try this with Canu and Flye, but it will be messy! Let's try running the flye assembly on one MSN node and canu on the rest. We'll see how that goes.

Even worse, there are signs that the blue pippen protist samples may have a contaminant because of the GC% profile! I will pretend like I didn't know this and try to assemble them separately. First, let's prepare the files for flye assembly and canu assembly.

> Ceres: /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover

```bash
# Very simple -- I'm going to concatenate only the high content protist runs
cat B_Pip_Proto.concatenated.fastq New_cell_proto.concatenated.fastq Proto_HMW1ug.concatenated.fastq Proto_HMW2ug.concatenated.fastq Proto_PCI_HMW_DNA.concatenated.fastq Proto_PCI_Spun_DNA.concatenated.fastq Proto_QiaT1.concatenated.fastq RCL_Proto1.concatenated.fastq RC_Proto1.concatenated.fastq > diagnostic_protist_combined.fastq
```

And now to run the two assemblies. Fingers crossed!

> Ceres: /project/rumen_longread_metagenome_assembly/assemblies/protists

```bash
module load miniconda/3.6 canu/1.8
source activate /KEEP/rumen_longread_metagenome_assembly/flye

sbatch --nodes=1 --mem=320000 --ntasks-per-node=70 -p msn -q msn flye -g 1.0g --nano-raw /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/diagnostic_protist_11_14_combined.fastq -t 70 -i 2 -m 4000 --asm-coverage 40 --meta -o flye_meta_diagprot1114

sbatch --nodes=1 --ntasks-per-node=2 --mem=10G -p short -q memlimit --wrap="canu -p canu_diagprotist -d canu_diagprotist genomeSize=1000m corOutCoverage=10000 corMhapSensitivity=high corMinCoverage=0 redMemory=32 oeaMemory=32 batMemory=200 'gridOptions=-p short -q memlimit' -nanopore-raw /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/diagnostic_protist_combined.fastq"

# OK, going to add more to the mix!
sbatch --nodes=1 --mem=1000 --ntasks-per-node=2 -p brief-low -q memlimit --wrap="cat diagnostic_protist_combined.fastq ULR_Bp_PCI.combined.fastq ULR_Zymo_3.combined.fastq ULR_Zymo_1.combined.fastq ULR_Bp_HMW.combined.fastq Zymo_proto_1.combined.fastq Zymo_proto_3.combined.fastq Blu_pip_proto_PCI.combined.fastq Blu_pip_proto_HMW.combined.fastq > diagnostic_protist_11_14_combined.fastq"


```

Now to align sequence reads to the contigs from our longread rumen metagenome dataset and then separate into mags using a generic tool.

> Ceres: /project/rumen_longread_metagenome_assembly/assemblies/protists

```bash
module load minimap2/2.6 samtools/1.9
sbatch --nodes=1 --mem=9000 --ntasks-per-node=1 -p short -q memlimit --wrap="module load bwa; bwa index flye_meta_diagprot1114/assembly.fasta"

sbatch --nodes=1 --mem=15000 --ntasks-per-node=4 -p short -q memlimit --wrap="minimap2 -ax map-ont flye_meta_diagprot1114/assembly.fasta /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/diagnostic_protist_11_14_combined.fastq | samtools sort -o flye_meta_diagprot1114.ontmaps.sort.bam -T temp_flye -; samtools index flye_meta_diagprot1114.ontmaps.sort.bam"

sbatch --nodes=1 --ntasks-per-node=1 -p brief-low -q memlimit --wrap="python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b flye_protist_aligns -t /project/rumen_longread_metagenome_assembly/sequence_data/ymprep_illumina_sequence_files.tab -f /project/rumen_longread_metagenome_assembly/assemblies/protists/flye_meta_diagprot1114/assembly.fasta -q memlimit -p short -m"
```

Now that I have the coverages from two different datasets, let's start the blobtools pipeline, metabat binning and MAGpy pipelines.

```bash
module load metabat/2.12.1 miniconda/3.6 usearch/11.0.667 hmmer3/3.2.1 blobtools/1.1.1 diamond/0.9.28

# Downloading blobtools diamond databases
sbatch -t 2-0 download_diamond_db.sh

# Generating coverage files
for i in flye_meta_diagprot1114.ontmaps.sort.bam flye_protist_aligns/YMpreprun3/YMpreprun3.sorted.merged.bam; do echo $i; sbatch --nodes=1 --mem=20000 --ntasks-per-node=2 -p brief-low -q memlimit --wrap="jgi_summarize_bam_contig_depths --outputDepth $i.cov $i"; done

# The minimap2 alignment of reads to the protist contigs didn't work very well! Curious... Going to just use the YMprep3 data

sbatch --nodes=1 -p brief-low -q memlimit --ntasks-per-node=15 --mem=50000 --wrap="metabat2 -i flye_meta_diagprot1114/assembly.fasta -a flye_protist_aligns/YMpreprun3/YMpreprun3.sorted.merged.bam.cov -o mags/fm14 -t 15 -v"

# Now let's try MAGpy. I need to copy over the config files and use them for setting up the session
cp /project/rumen_longread_metagenome_assembly/binaries/MAGpy/config.json ./
cp /project/rumen_longread_metagenome_assembly/binaries/MAGpy/MAGpy.json ./
cp -r /project/rumen_longread_metagenome_assembly/binaries/MAGpy/scripts ./
cp -r /project/rumen_longread_metagenome_assembly/binaries/MAGpy/envs ./
# For sourmash (because pip installed it in the worst spot)
export PATH=/home/derek.bickharhth/.local/bin:$PATH
source activate /KEEP/rumen_longread_metagenome_assembly/magpy/

# MAGpy run
# --conda --conda-prefix /KEEP/rumen_longread_metagenome_assembly/magpy
sbatch --nodes=1 --mem=5000 --ntasks-per-node=1 -p short -q memlimit snakemake  -s /project/rumen_longread_metagenome_assembly/binaries/MAGpy/MAGpy --cluster-config MAGpy.json --cluster 'sbatch --nodes=1 --ntasks-per-node={cluster.core} --mem={cluster.vmem} -p {cluster.proj} -q memlimit' --jobs 1000

# Diamond run before blobtools
sbatch -t 2-0 -p msn -q msn --nodes=1 --ntasks-per-node=30 --mem=100000 --wrap="diamond blastx --query flye_meta_diagprot1114/assembly.fasta --db uniprot_ref_proteomes.diamond.dmnd --threads 29 --outfmt 6 --sensitive --max-target-seqs 1 --evalue 1e-25 -o flye_meta_diagprot1114.diamondout.tsv"

# Blobtools taxify
sbatch -t 2-0 --nodes=1 --mem=20000 --ntasks-per-node=3 -p mem -q memlimit --wrap='blobtools taxify -f flye_meta_diagprot1114.diamondout.tsv -m uniprot_ref_proteomes.taxids -s 0 -t 2 -o flye_meta_diagprot1114_diamond_uniprot'

# Blobtools create
sbatch -t 2-0 --nodes=1 --mem=10000 --ntasks-per-node=1 -p mem -q memlimit --wrap='blobtools create -i flye_meta_diagprot1114/assembly.fasta -b flye_protist_aligns/YMpreprun3/YMpreprun3.sorted.merged.bam -t flye_meta_diagprot1114_diamond_uniprot.flye_meta_diagprot1114.diamondout.tsv.taxified.out -o flye_protist_diagprot_blobtools --db blob_ncbi.db'

sbatch -t 2-0 --nodes=1 --mem=10000 --ntasks-per-node=1 -p mem -q memlimit --wrap='blobtools plot -i flye_protist_diagprot_blobtools.blobDB.json --notitle -r superkingdom -o flye_pd_blobtools_supkingdom'

mkdir blob_plots
mv flye_pd_blobtools_supkingdom.flye_protist_diagprot_blobtools.blobDB* ./blob_plots/
mv flye_meta_diagprot1114_diamond_uniprot.flye_meta_diagprot1114.diamondout.tsv.taxified.out ./blob_plots/
mv flye_protist_diagprot_blobtools.YMpreprun3.sorted.merged.bam.cov ./blob_plots/

# I didn't get "blobs" but this is a known problem. the dingbats at CERES made the default python version 2.7 for blobtools!
# I need to run it with python3 from the get-go or else it silently fails.
```

## Actual assemblies

OK, we've generated a ton of data on these samples, let's start assembling!

#### NOTE: to clean up the workspace, I moved most of the preliminary protist analysis to a new folder! diag_protist

> Ceres: /project/rumen_longread_metagenome_assembly/assemblies/protists

```bash
# flye assembly
module load miniconda/3.6
source activate /KEEP/rumen_longread_metagenome_assembly/flye

### 7201 ###
sbatch --nodes=1 --mem=300000 --ntasks-per-node=70 -p msn -q msn flye -g 1.0g --nano-raw /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/cow7201_total_combined.2019.fastq -t 70 -m 5000 --meta --plasmids -o flye_meta_7201 --resume

# Read stats for future reference:
[2019-12-27 17:37:44] INFO: Total read length: 134283567513
[2019-12-27 17:37:44] INFO: Input genome size: 1000000000
[2019-12-27 17:37:45] INFO: Estimated coverage: 134
[2019-12-27 17:37:45] INFO: Reads N50/N90: 16874 / 2519

# Also for future reference at 90% overlap of disjointigs
grep 'Length:' /project/rumen_longread_metagenome_assembly/assemblies/protists/flye_meta_7201/flye.log | perl -lane 'print $F[1];' | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   10,020
Sum:    1,387,722,782
Minimum 11,203
Maximum 1,217,044
Average 138,495.287625
Median  115,649
Standard Deviation      101,898.306042
Mode(Highest Distributed Value) 39,410

# The plasmids command took 7 days and produced an 8.9 terabyte paf file! It's an all-vs-all alignment using minimap2! Let's drop that for now.
# I had to change the "stage_name" to "consensus" in the params.json file in the metaflye folder.
sbatch --nodes=1 --mem=300000 --ntasks-per-node=70 -p msn -q msn flye -g 1.0g --nano-raw /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/cow7201_total_combined.2019.fastq -t 70 -m 5000 --meta -o flye_meta_7201 --resume

### 8615 ###
sbatch --nodes=1 --mem=600000 --ntasks-per-node=70 -p mem -q memlimit flye -g 2.0g --nano-raw /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/Cow8615_10.raw.fastq.gz /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/Cow8615_11.raw.fastq.gz /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/Cow8615_12.raw.fastq.gz /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/Cow8615_13.raw.fastq.gz /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/Cow8615_14.raw.fastq.gz /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/Cow8615_1.raw.fastq.gz /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/Cow8615_2.raw.fastq.gz /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/Cow8615_3.raw.fastq.gz /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/Cow8615_4.raw.fastq.gz /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/Cow8615_5.raw.fastq.gz /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/Cow8615_6.raw.fastq.gz /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/Cow8615_7.raw.fastq.gz /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/Cow8615_8.raw.fastq.gz /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/Cow8615_9.raw.fastq.gz -t 70 -m 2000 --meta -o flye_meta_8615 --resume

### 8484 ###
sbatch --nodes=1 --mem=600000 --ntasks-per-node=80 -p priority-mem -q msn-mem flye -g 1.0g --nano-raw /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/Cow8484_10.raw.fastq.gz /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/Cow8484_11.raw.fastq.gz /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/Cow8484_12.raw.fastq.gz /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/Cow8484_13.raw.fastq.gz /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/Cow8484_14.raw.fastq.gz /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/Cow8484_15.raw.fastq.gz /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/Cow8484_16.raw.fastq.gz /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/Cow8484_1.raw.fastq.gz /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/Cow8484_2.raw.fastq.gz /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/Cow8484_3.raw.fastq.gz /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/Cow8484_4.raw.fastq.gz /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/Cow8484_5.raw.fastq.gz /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/Cow8484_6.raw.fastq.gz /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/Cow8484_7.raw.fastq.gz /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/Cow8484_8.raw.fastq.gz /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/Cow8484_9.raw.fastq.gz -t 70 -m 3000 --meta -o flye_meta_8484 --resume

```

## Assembly informatics

> Ceres: /project/rumen_longread_metagenome_assembly/assemblies/protists

```bash
module load diamond/0.9.28

## 7201 ##
sbatch -t 2-0 -p msn -q msn --nodes=1 --ntasks-per-node=30 --mem=100000 --wrap="diamond blastx --query flye_meta_7201/disjointig_metaflye_7201.fasta --db uniprot_ref_proteomes.diamond.dmnd --threads 29 --outfmt 6 --sensitive --max-target-seqs 1 --evalue 1e-25 -o flye_meta_7201.disjointig.diamondout.tsv"
sbatch -t 2-0 -p msn -q msn --nodes=1 --ntasks-per-node=30 --mem=100000 --wrap="diamond blastx --query flye_meta_7201/filtered_metaflye_7201.fasta --db uniprot_ref_proteomes.diamond.dmnd --threads 29 --outfmt 6 --sensitive --max-target-seqs 1 --evalue 1e-25 -o flye_meta_7201.filtered.diamondout.tsv"

sbatch --nodes=1 --mem=10000 --ntasks-per-node=5 -p short -q memlimit --wrap="python3 ~/python_toolchain/sequenceData/calcGCcontentFasta.py -f flye_meta_8484/assembly.fasta -o flye_meta_8484/assembly.fasta.gc"
sbatch --nodes=1 --mem=5000 --ntasks-per-node=2 -p short -q memlimit --wrap='cat /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/[Cc]ow8484_*/*/fastq_pass/*.fastq > cow8484_total_raw.fastq'

sbatch --nodes=1 --mem=45000 --ntasks-per-node=10 -p short -q memlimit --wrap="python3 ~/python_toolchain/sequenceData/calcGCcontentFasta.py -q cow8484_total_raw.fastq -o cow8484_total_raw.fastq.gc -t 10"

# comparing to older datasets
sbatch --nodes=1 --mem=45000 --ntasks-per-node=10 -p short -q memlimit --wrap="python3 ~/python_toolchain/sequenceData/calcGCcontentFasta.py -f /project/rumen_longread_metagenome_assembly/sequence_data/pilot_project/pacbio/rumen_pacbio_corrected.fasta -o rumen_pacbio_corrected.fasta.gc -t 10"
sbatch --nodes=1 --mem=45000 --ntasks-per-node=10 -p short -q memlimit --wrap="python3 ~/python_toolchain/sequenceData/calcGCcontentFasta.py -f /project/rumen_longread_metagenome_assembly/assemblies/pilot_project/pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa -o usda_pacbio_second_pilon_indelsonly.fa.gc -t 10"

# Formatting for R
perl -ne '@F = split(/\t/); $a = "protist_asm"; print "$a\t$F[1]\t$F[2]";' < flye_meta_8484.asm.gc > flye_meta_8484.asm.gc.format
perl -ne '@F = split(/\t/); $a = "protist_reads"; print "$a\t$F[1]\t$F[2]";' < cow8484_total_raw.fastq.gc >cow8484_total_raw.fastq.gc.format
perl -ne '@F = split(/\t/); $a = "pacbio_asm"; print "$a\t$F[1]\t$F[2]";' < usda_pacbio_second_pilon_indelsonly.fa.gc > usda_pacbio_second_pilon_indelsonly.fa.gc.format
perl -ne '@F = split(/\t/); $a = "pacbio_reads"; print "$a\t$F[1]\t$F[2]";' < rumen_pacbio_corrected.fasta.gc > rumen_pacbio_corrected.fasta.gc.format
cat cow8484_total_raw.fastq.gc.format flye_meta_8484.asm.gc.format rumen_pacbio_corrected.fasta.gc.format usda_pacbio_second_pilon_indelsonly.fa.gc.format > combined_dataset_gc.format
```

```R
library(ggplot2)
library(dplyr)
data <- read.delim("combined_dataset_gc.format", header=FALSE)
colnames(data) <- c("dataset", "length", "GC")

pdf(file="combined_dataset_gc_plot.pdf", useDingbats=FALSE)
p <- ggplot(data, aes(x=dataset, y=GC, fill=dataset)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, outlier.shape = NA, fill = "white")
p + labs(y="Avg GC ratio") + scale_fill_brewer(palette="Dark2") + theme_bw() + facet_wrap(~ dataset, nrow=2)
dev.off()

data <- mutate(data, Quant = cut(data$GC, quantile(data$GC, probs=0:5/5), include.lowest=TRUE))
data %>% group_by(Quant, dataset) %>% summarize(num = n(), sum=sum(as.numeric(length)))
# A tibble: 20 x 4
# Groups:   Quant [5]
   Quant         dataset           num         sum
   <fct>         <fct>           <int>       <dbl>
 1 [0,0.215]     pacbio_asm        142      435070
 2 [0,0.215]     pacbio_reads    78628   376325101
 3 [0,0.215]     protist_asm     16854   283287743
 4 [0,0.215]     protist_reads 8841159 46704737738
 5 (0.215,0.245] pacbio_asm        319     1427805
 6 (0.215,0.245] pacbio_reads   257300  1738547311
 7 (0.215,0.245] protist_asm      1456    29882964
 8 (0.215,0.245] protist_reads 8677347 51069841819
 9 (0.245,0.308] pacbio_asm       2821    28744257
10 (0.245,0.308] pacbio_reads   645409  4922436360
11 (0.245,0.308] protist_asm      1521    73756635
12 (0.245,0.308] protist_reads 8296054 20313470822
13 (0.308,0.463] pacbio_asm      25076   337785439
14 (0.308,0.463] pacbio_reads  1750843 14126558510
15 (0.308,0.463] protist_asm      5466   260584557
16 (0.308,0.463] protist_reads 7146962 12301776828
17 (0.463,0.98]  pacbio_asm      49312   708033673
18 (0.463,0.98]  pacbio_reads  3199502 26039567030
19 (0.463,0.98]  protist_asm     11608   642261921
20 (0.463,0.98]  protist_reads 5675230 15469355257
```

#### Hi-C library alignment

```bash
module load bwa samtools

for i in flye_meta_7201/*.fasta; do echo $i; sbatch --nodes=1 --mem=12000 --ntasks-per-node=1 -t 1-0 -p msn -q msn --wrap="bwa index $i"; done

sbatch --nodes=1 --mem=12000 --ntasks-per-node=2 -t 1-0 -p msn -q msn --dependency=afterok:1598335 --wrap="bwa mem flye_meta_7201/disjointig_metaflye_7201.fasta /project/rumen_longread_metagenome_assembly/sequence_data/protist_hic/rumen_7201_qc_S3HiC_R1.fastq.gz /project/rumen_longread_metagenome_assembly/sequence_data/protist_hic/rumen_7201_qc_S3HiC_R2.fastq.gz | samtools sort -T disjointig -o disjointig_metaflye_7201.hic.bam -"

sbatch --nodes=1 --mem=12000 --ntasks-per-node=2 -t 1-0 -p msn -q msn --dependency=afterok:1598336 --wrap="bwa mem flye_meta_7201/filtered_metaflye_7201.fasta /project/rumen_longread_metagenome_assembly/sequence_data/protist_hic/rumen_7201_qc_S3HiC_R1.fastq.gz /project/rumen_longread_metagenome_assembly/sequence_data/protist_hic/rumen_7201_qc_S3HiC_R2.fastq.gz | samtools sort -T filtered -o filtered_metaflye_7201.hic.bam -"

for i in disjointig_metaflye_7201.hic.bam filtered_metaflye_7201.hic.bam; do samtools index $i; done

conda activate /KEEP/rumen_longread_metagenome_assembly/hic_qc/
for i in *.bam; do echo $i; sbatch --nodes=1 --mem=10000 --ntasks-per-node=2 -p priority -q msn -t 1-0 --wrap="python ~/rumen_longread_metagenome_assembly/binaries/hic_qc/hic_qc.py -b $i -n -1 -o $i.hic -r --sample_type metagenome"; done
```

#### Downloading ENA dataset 

> Ceres: /project/rumen_longread_metagenome_assembly/assemblies/mick_rug

```bash
sbatch --nodes=1 --mem=10000 --ntasks-per-node=1 -p priority -q msn --wrap="/project/rumen_longread_metagenome_assembly/binaries/enaBrowserTools/python3/enaGroupGet -g analysis -f submitted -t PRJEB31266"
```