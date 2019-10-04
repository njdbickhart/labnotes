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
```