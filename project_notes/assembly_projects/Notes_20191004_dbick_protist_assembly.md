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

sbatch --nodes=1 --mem=320000 --ntasks-per-node=70 -p msn -q msn flye -g 1.0g --nano-raw /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/diagnostic_protist_combined.fastq -t 70 -i 2 -m 4000 --asm-coverage 40 --meta -o flye_meta_diagprotist

sbatch --nodes=1 --ntasks-per-node=2 --mem=10G -p short -q memlimit --wrap="canu -p canu_diagprotist -d canu_diagprotist genomeSize=1000m corOutCoverage=10000 corMhapSensitivity=high corMinCoverage=0 redMemory=32 oeaMemory=32 batMemory=200 'gridOptions=-p short -q memlimit' -nanopore-raw /project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/diagnostic_protist_combined.fastq"
```