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

sbatch --nodes=1 --mem=450000 --ntasks-per-node=70 -p mem -q memlimit flye -g 1.0g --nano-raw /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/cow7201_total_combined.2019.fastq -t 70 -i 2 -m 4000 --meta --plasmids -o flye_meta_7201

# Read stats for future reference:
[2019-12-27 17:37:44] INFO: Total read length: 134283567513
[2019-12-27 17:37:44] INFO: Input genome size: 1000000000
[2019-12-27 17:37:45] INFO: Estimated coverage: 134
[2019-12-27 17:37:45] INFO: Reads N50/N90: 16874 / 2519

```