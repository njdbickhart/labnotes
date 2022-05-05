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

awk '/^S/{print ">"$2"\n"$3}' /lustre/project/rumen_longread_metagenome_assembly/analysis/buccal/test_juliana/Genome_Assembly/HIFI/metaflye/assembly_graph.gfa | fold > protist_graph_edges.fa

# 70% of reads had a seed, but only 1.7% of reads had an alignment! I wonder if starting with the metaflye graph was the problem here. 
```

OK, back to the drawing board. Do we need to start from a different assumption, trim the graph some more or contact the developers. I think that I would be violating too many assumptions if I went straight into the MGB approach, so let's try to have a conversation with the assembly experts.

## Gathering assembly stats

> Ceres: /lustre/project/rumen_longread_metagenome_assembly/analysis/buccal/test_juliana/Genome_Assembly/HIFI

```bash
# preparing a set of the top contigs > 100 kb for analysis
rm metaflye_top_contigs.fasta; perl -lane 'if($F[0] =~ /\#/){next;} if($F[1] > 100000){system("samtools faidx metaflye/assembly.fasta $F[0] >> metaflye_top_contigs.fasta");}' < metaflye/assembly_info.txt

samtools faidx metaflye_top_contigs.fasta
wc -l metaflye_top_contigs.fasta.fai
	739 metaflye_top_contigs.fasta.fai


module load python_3/3.6.6 miniconda/3.6 samtools
#OK, let's coopt my downsample snakemake workflow that I used for the previous HIFI manuscript to start gathering as much information as possible here.
mkdir protist_class_wkflow

mv protist_class_wkflow/Downsample protist_class_wkflow/Snakefile
mkdir protist_class_wkflow/logs

sbatch -N 1 -n 2 --mem=10000 -p priority -q msn snakemake -s ./protist_class_wkflow/Snakefile --cluster-config ~/python_toolchain/snakeMake/hifiMAGManuscript/cluster.json --cluster "sbatch -N 1 --ntasks-per-node={cluster.ntasks-per-node} --mem={cluster.mem} -p priority -q msn -o {cluster.stdout}" -p --use-conda --jobs 250 --verbose --latency-wait 40
```

## Working notes

Attempted to characterize potential mitochondrial/hydrogenosome sequence. Best hit was to contig_22275, but hits were of dubious quality (27% divergence to other protist mitochondria). 

## Minimizer graph construction

OK, so I have a metaFlye assembly that was only "OK" and managed to assemble more bacteria in higher contigs. I suspect that a different k-value is needed to extend the assembly across highly repetitive sequence in the protists. I will pull only the HiFi reads that mapped to protist "nodes" in the assembly graph and then test out several different values of k in order to optimize the parameters for assembly.

Once I have the new graph, I could try taking the nodes and using them as "seeds" for additional metagenome assembly of the protists. 

Validation data could include some of Sharon Huws' [cDNA RNA-seq](https://www.ncbi.nlm.nih.gov/sra/?term=SAMN13506237). 

Let's pull reads from the nodes first and then characterize them.

> Ceres: /project/rumen_longread_metagenome_assembly/assembly/protist

```bash
module load python_3/3.6.6 miniconda/3.6 samtools minimap2 r/4.1.2
cat nodes.txt | perl -lane 'foreach $f (@F){$f =~ s/,//; $f =~ s/edge_//; print "$f";}' > protist_nodes.list

perl -e 'chomp(@ARGV); %h = {}; open(IN, "< $ARGV[0]"); while(<IN>){chomp; $h{$_} = 1;} close IN; open(IN, "< $ARGV[1]"); <IN>; while(<IN>){chomp; @s = split(/\t/); @j = split(/,/, $s[-1]); foreach $j (@j){if(exists($h{$j})){print "$s[0]\n"; last;}}} close IN;' protist_nodes.list /OLD/project/rumen_longread_metagenome_assembly/analysis/buccal/test_juliana/Genome_Assembly/HIFI/metaflye/assembly_info.txt > protist_contigs.list

samtools faidx -r protist_contigs.list /OLD/project/rumen_longread_metagenome_assembly/analysis/buccal/test_juliana/Genome_Assembly/HIFI/metaflye/assembly.fasta > protist_contigs.fasta

mkdir /90daydata/rumen_longread_metagenome_assembly/analysis/protist
mkdir /90daydata/rumen_longread_metagenome_assembly/analysis/protist/protist_contig_aligns
for i in /OLD/project/rumen_longread_metagenome_assembly/sequence_data/protist_ccs/*fastq*; do name=`basename $i | cut -d'.' -f1`; echo $name; sbatch -N 1 -n 5 -p priority -q msn --mem=22000 --wrap="minimap2 -x map-hifi -t 5 -o /90daydata/rumen_longread_metagenome_assembly/analysis/protist/protist_contig_aligns/$name.alignments.paf protist_contigs.fasta $i"; done;

mkdir protist_fastq_stats
conda activate /project/rumen_longread_metagenome_assembly/environments/seaborn/
for i in /OLD/project/rumen_longread_metagenome_assembly/sequence_data/protist_ccs/*fastq*; do name=`basename $i | cut -d'.' -f1`; echo $name; sbatch -N 1 -n 2 -p priority -q msn --mem=30000 --wrap="python3 ~/python_toolchain/sequenceData/fastqStats.py -f $i -o protist_fastq_stats/$name"; done

mkdir /90daydata/rumen_longread_metagenome_assembly/analysis/protist/protist_contig_rlists
for i in /OLD/project/rumen_longread_metagenome_assembly/sequence_data/protist_ccs/*fastq*; do name=`basename $i | cut -d'.' -f1`; echo $name; sbatch -N 1 -n 5 -p priority -q msn --mem=22000 --wrap="python3 filterHifiRDNAByQPafsF.py /90daydata/rumen_longread_metagenome_assembly/analysis/protist/protist_contig_aligns/$name.alignments.paf /90daydata/rumen_longread_metagenome_assembly/analysis/protist/protist_contig_rlists/$name.reads.list"; done;

conda activate /project/rumen_longread_metagenome_assembly/environments/seaborn/
mkdir /90daydata/rumen_longread_metagenome_assembly/analysis/protist/protist_contig_reads
for i in /OLD/project/rumen_longread_metagenome_assembly/sequence_data/protist_ccs/*fastq*; do name=`basename $i | cut -d'.' -f1`; echo $name; sbatch -N 1 -n 2 -p priority -q msn --mem=25000 --wrap="python3 ~/python_toolchain/sequenceData/filterFastaqFromList.py -q $i -l /90daydata/rumen_longread_metagenome_assembly/analysis/protist/protist_contig_rlists/$name.reads.list -o /90daydata/rumen_longread_metagenome_assembly/analysis/protist/protist_contig_reads/$name.rdna.fastq; gzip /90daydata/rumen_longread_metagenome_assembly/analysis/protist/protist_contig_reads/$name.rdna.fastq"; done
```

OK, I've kept only 20-35% of the total number of HiFi reads this time. I don't know how lengthy the repeat segments will be, so let's try a couple kmer histogram plots to try to suss things out.

> Ceres: /90daydata/rumen_longread_metagenome_assembly/analysis/protist/protist_contig_reads/

```bash
module load jellyfish2/2.2.9

# Ugh, need to unzip these first
for i in *.fastq.gz; do sbatch -N 1 -n 2 --mem=5000 -p priority -q msn --wrap="gunzip $i"; done

for i in *.fastq; do echo -n "$i "; done; echo

for i in 14 21 28; do echo $i; sbatch --nodes=1 --mem=100000 --ntasks-per-node=15 --wrap="jellyfish count -m $i -s 90G -t 15 -o protist_reads_${i}.jf -C m54337U_210722_195630.rdna.fastq m54337U_210802_201210.rdna.fastq m54337U_210809_182434.rdna.fastq m54337U_211013_163649.rdna.fastq m54337U_211015_151518.rdna.fastq m54337U_211017_003233.rdna.fastq"; done
```

Let's queue up some preliminary graphs.

> Ceres: /90daydata/rumen_longread_metagenome_assembly/analysis/protist

```bash
module load python_3/3.6.6 miniconda/3.6 samtools minimap2 r/4.1.2
conda activate /project/rumen_longread_metagenome_assembly/environments/verkko
mkdir protist_mbg

for i in protist_contig_reads/*.fastq; do echo -n "-i $i "; done; echo
sbatch -N 1 -n 35 --mem=300000 -p priority -q msn -o mbg_2001_1500.slurm.out --wrap="MBG -t 35 -k 2001 -w 1500 -r 15000 --output-sequence-paths protist_mbg/protist_2001_1500_paths.gaf  --out protist_mbg/protist_2001_1500_graph.gfa -i protist_contig_reads/m54337U_210722_195630.rdna.fastq -i protist_contig_reads/m54337U_210802_201210.rdna.fastq -i protist_contig_reads/m54337U_210809_182434.rdna.fastq -i protist_contig_reads/m54337U_211013_163649.rdna.fastq -i protist_contig_reads/m54337U_211015_151518.rdna.fastq -i protist_contig_reads/m54337U_211017_003233.rdna.fastq"

sbatch -N 1 -n 35 --mem=300000 -p priority -q msn -o mbg_3501_3000.slurm.out --wrap="MBG -t 35 -k 3501 -w 3000 -r 15000 --output-sequence-paths protist_mbg/protist_3501_3000_paths.gaf  --out protist_mbg/protist_3501_3000_graph.gfa -i protist_contig_reads/m54337U_210722_195630.rdna.fastq -i protist_contig_reads/m54337U_210802_201210.rdna.fastq -i protist_contig_reads/m54337U_210809_182434.rdna.fastq -i protist_contig_reads/m54337U_211013_163649.rdna.fastq -i protist_contig_reads/m54337U_211015_151518.rdna.fastq -i protist_contig_reads/m54337U_211017_003233.rdna.fastq"

# Let's try some different of k as the 3501 graph was resolved, but a bit too separate in nodes!
for i in 501 1001 2501 4001 4501 5001; do j=`echo $i | perl -lane '$F[0] = int($F[0] / 2); print $F[0];'`; echo "$i $j"; sbatch -N 1 -n 35 --mem=300000 -p priority -q msn -o mbg_${i}_${j}.slurm.out --wrap="MBG -t 35 -k $i -w $j -r 15000 --output-sequence-paths protist_mbg/protist_${i}_${j}_paths.gaf  --out protist_mbg/protist_${i}_${j}_graph.gfa -i protist_contig_reads/m54337U_210722_195630.rdna.fastq -i protist_contig_reads/m54337U_210802_201210.rdna.fastq -i protist_contig_reads/m54337U_210809_182434.rdna.fastq -i protist_contig_reads/m54337U_211013_163649.rdna.fastq -i protist_contig_reads/m54337U_211015_151518.rdna.fastq -i protist_contig_reads/m54337U_211017_003233.rdna.fastq"; done 
```

Testing a run of HiFiasm-Meta on the reads just for giggles.

```bash
module load python_3/3.6.6 miniconda/3.6 samtools minimap2 r/4.1.2
sbatch -N 1 -n 70 --mem=1000000 -p priority-mem -q msn-mem --wrap="/project/rumen_longread_metagenome_assembly/binaries/hifiasm-meta/hifiasm_meta -o protist_hifiasm -t 70 protist_contig_reads/m54337U_210722_195630.rdna.fastq protist_contig_reads/m54337U_210802_201210.rdna.fastq protist_contig_reads/m54337U_210809_182434.rdna.fastq protist_contig_reads/m54337U_211013_163649.rdna.fastq protist_contig_reads/m54337U_211015_151518.rdna.fastq protist_contig_reads/m54337U_211017_003233.rdna.fastq"
```