# Protist HiFi assembly
---
*11/15/2021*

## Table of Contents

## Recap of prior results

Juliana was able to run metaFlye on the data, but obtained a terrible assembly despite the addition of new data. She also ran a prior version of hifiasm-meta and it crashed with a core dump. I need to find out if the reads are uniform in content before engaging in additional assemblies.

## Testing read dataset uniformity

Let's start out with a k-mer based test and then move to more advanced diagnostics. If one read set has completely different (> 5%) k-mers from the rest, that may impact assembly quality.

> Ceres: 

```bash
module load miniconda/3.6
conda activate /KEEP/rumen_longread_metagenome_assembly/meryl/

# creating the meryl dbs
mkdir meryl_7201; for i in /lustre/project/rumen_longread_metagenome_assembly/sequence_data/protist_ccs/*; do name=`basename $i | cut -d'.' -f1`; echo $name; sbatch -N 1 -n 35 --mem=150000 -p priority -q msn -t 2-0 --wrap="/lustre/project/rumen_longread_metagenome_assembly/binaries/meryl/build/bin/meryl k=21 count output meryl_7201/$name.meryl threads=35 memory=145 $i"; done

cd meryl_7201/
sbatch -N 1 -n 5 -p priority -q msn --mem=50000 -t 1-0 --wrap="python3 ~/python_toolchain/sequenceData/merylVennUpset.py -m /lustre/project/rumen_longread_metagenome_assembly/binaries/meryl/build/bin/meryl -d m54337U_210722_195630.meryl -d m54337U_210802_201210.meryl -d m54337U_210809_182434.meryl -d m54337U_211013_163649.meryl -d m54337U_211015_151518.meryl -d m54337U_211017_003233.meryl -o protist_7201_21mer_test"
cd ../

conda deactivate
# testing out one of the very odd fastq files that popped out of kmer testing

conda activate /KEEP/rumen_longread_metagenome_assembly/seaborn/
for i in m54337U_210809_182434.hifi_reads.fastq.gz m54337U_211013_163649.hifi_reads.fastq.gz m54337U_211015_151518.hifi_reads.fastq.gz; do name=`echo $i | cut -d'.' -f1`; echo $name; sbatch -N 1 -n 2 --mem=25000 -p priority -q msn --wrap="python3 ~/python_toolchain/sequenceData/ontFastqCompare.py -f /lustre/project/rumen_longread_metagenome_assembly/sequence_data/protist_ccs/m54337U_210802_201210.hifi_reads.fastq -s /lustre/project/rumen_longread_metagenome_assembly/sequence_data/protist_ccs/$i -o $name\_compare"; done
```

## Graph based assessment

I've checked Juliana's metaFlye assembly graph for some clues on how it assessed the protist genomes, and there is a sufficiently dense hairball that makes up 778 Mbp in the graph. I suspect that repeat structure really wrecked havoc here! I will try to use ONT reads to span the gaps and scaffold the contigs, using the graph alignment methodology that the T2T project uses.

> Ceres: /lustre/project/rumen_longread_metagenome_assembly/assemblies/protists

```bash
module load python_3/3.6.6 miniconda/3.6

conda activate /KEEP/rumen_longread_metagenome_assembly/graph_tools

#NOTE:  I had to remove allot of extra parameters that the T2T guys added in here. Not sure if they were using a different branch or modified version of GraphAligner 
sbatch -N 1 -n 70 --mem=300000 -p priority -q msn --wrap="GraphAligner -t 70 -g /lustre/project/rumen_longread_metagenome_assembly/analysis/buccal/test_juliana/Genome_Assembly/HIFI/metaflye/assembly_graph.gfa -f /lustre/project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/Cow7201_proto_data/Cow7201_total_combined.2019.fastq.gz -a cow7201_protist_hififlye_ontalign.gaf --seeds-mxm-length 30 --seeds-mem-count 1000 -b 15 --multimap-score-fraction 0.99 --precise-clipping 0.85 --min-alignment-score 5000"

python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f cow7201_protist_hififlye_ontalign.gaf -c 5 -d '\t' | perl -ne '$num = () = $_ =~ /edge_/g; if($num > 1){print $_;}' | head -n 100
```