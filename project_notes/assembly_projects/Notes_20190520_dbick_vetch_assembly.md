# Vetch Genome assembly
---
*5/20/2019*

## Table of Contents

## Testing the variance of the different Vetch libraries

OK, so Lisa tried out several different methods for diluting the DNA preps prior to library creation, and we need to see if there are any obvious biases in these methods that could cause us to choose to avoid them in the future!

We're just going to do the Mash screen to compare them (and to compare them against a red clover sample as an outgroup!).

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/sequence_data/red_clover_nano

```bash
# Creating a merged fastq file
for i in VetchBooger VetchMN VetchVortex VetchZymo; do echo $i; cat $i/*/fastq_pass/*.fastq > $i/$i.combined.pass.fastq; done

for i in VetchBooger VetchMN VetchVortex VetchZymo; do echo $i; sbatch --nodes=1 --mem=10000 --ntasks-per-node=4 -p msn --wrap="~/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.0/mash sketch -o $i/$i.combined.msh -p 4 -s 100000 -r -m 4 -g 420M $i/$i.combined.pass.fastq"; done
~/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.0/mash paste vetch_combined VetchBooger/VetchBooger.combined.msh VetchMN/VetchMN.combined.msh VetchZymo/VetchZymo.combined.msh VetchVortex/VetchVortex.combined.msh clover14.msh
~/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.0/mash dist -t vetch_combined.msh VetchBooger/VetchBooger.combined.msh VetchMN/VetchMN.combined.msh VetchZymo/VetchZymo.combined.msh VetchVortex/VetchVortex.combined.msh clover14.msh > vetch_combined.dist

perl -ne '$_ =~ s/(Vetch.+)\/.+fastq/$1/g; print $_;' < vetch_combined.dist > vetch_combined.rfmt.dist

```

```R
library(dplyr)
library(ggplot2)
library(reshape)


vdist <- read.delim("vetch_combined.rfmt.dist", header=TRUE)
rownames(vdist) <- vdist$query
vdist <- vdist[,c(2:6)]
vdist.m <- melt(as.matrix(vdist))

pdf(file="vetch_distance_heatmap.pdf", useDingbats=FALSE)
ggplot(vdist.m, aes(X1, X2)) + geom_tile(aes(fill=value), colour = "white") + scale_fill_gradient(low = "white", high = "steelblue") + labs(title = "Vetch Prep Dataset Mash Distances")
dev.off()
```

## Starting the assembly process

I'm going to unpack what I have and run an assembly over the long weekend.

> Ceres: /home/derek.bickharhth/forage_assemblies/sequence_data

```bash
for i in Vetch*.tar.gz; do echo $i; sbatch --nodes=1 --mem=1000 --ntasks-per-node=1 -p msn --wrap="tar -xvf $i"; done

# Now to concatenate the fastqs
sbatch --nodes=1 --mem=5000 --ntasks-per-node=1 -p msn --wrap='for i in Vetch*/*/fastq_pass/*.fastq; do cat $i; done > vetch_combined_reads.fastq'
```

Now to run Flye on Vetch to test it out.

> Ceres: /project/forage_assemblies/assemblies/vetch

```bash
module load miniconda/3.6
conda activate /KEEP/rumen_longread_metagenome_assembly/flye

sbatch --nodes=1 --mem=300000 --ntasks-per-node=70 -p msn flye --nano-raw /home/derek.bickharhth/forage_assemblies/sequence_data/vetch_combined_reads.fastq -g 2000m -t 70 -i 2 -o vetch_flye

## Ran out of memory! Going to try to increase a slight bit, but if it fails I'll have to try an alternate strategy
sbatch --nodes=1 --mem=370000 --ntasks-per-node=70 -p msn flye --nano-raw /home/derek.bickharhth/forage_assemblies/sequence_data/vetch_combined_reads.fastq -g 2000m -t 70 -i 2 -o vetch_flye --resume

## That also ran out of memory. Going to queue up on the MEM nodes for the long haul
sbatch --nodes=1 --mem=900000 --ntasks-per-node=70 -p mem flye --nano-raw /home/derek.bickharhth/forage_assemblies/sequence_data/vetch_combined_reads.fastq -g 2000m -t 70 -i 2 -o vetch_flye --resume

## I am going to try a co-assembly approach with the Noble vetch reads as well
sbatch --nodes=1 --mem=900000 --ntasks-per-node=70 -p priority-mem -q msn-mem flye --nano-raw /project/forage_assemblies/sequence_data/combined_both_vetch.fastq.gz -g 2000m -t 70 -i 2 -m 8000 --asm-coverage 60 -o vetch_limit_flye --resume
```

Now trying Canu

```bash
module load canu/1.8 java/1.8.0_121

# De novo correction
sbatch --nodes=1 --ntasks-per-node=10 --mem=10G -p short --wrap="canu -correct -p vetch_pre -d vetch_pre genomeSize=2000m corMhapSensitivity=high corMinCoverage=0 corOutCoverage=100 saveReadCorrections=true 'gridOptions=-p short' -nanopore-raw /project/forage_assemblies/sequence_data/vetch_combined_reads.fastq"

### NOT RUN YET replace fastq###
sbatch --nodes=1 --ntasks-per-node=30 --mem=10G -p short --wrap="canu -p vetch_canu -d vetch_canu genomeSize=2000m correctedErrorRate=0.120 'corMhapOptions=--threshold 0.8 --num-hashes 512 --ordered-sketch-size 1000 --ordered-kmer-size 14' 'gridOptions=-p short' -nanopore-raw /project/forage_assemblies/sequence_data/vetch_combined_reads.fastq"

# No correction; assume "polyploid"
sbatch --nodes=1 --ntasks-per-node=10 --mem=10G -p short --wrap="canu -p vetch_nocorrect -d vetch_nocorrect genomeSize=2000m corMhapSensitivity=normal corOutCoverage=200 'batOptions=-dg 3 -db 3 -dr 1 -ca 500 -cp 50' saveReadCorrections=true 'gridOptions=-p short' -nanopore-raw /project/forage_assemblies/sequence_data/vetch_combined_reads.fastq"

# Combined; new version of canu; 
sbatch --nodes=1 --ntasks-per-node=10 --mem=10G -p short -q memlimit --wrap="canu -p vetch_canu_combined -d vetch_canu_combined genomeSize=2000m corMhapSensitivity=normal corOutCoverage=200 'batOptions=-dg 3 -db 3 -dr 1 -ca 500 -cp 50' saveReadCorrections=true 'gridOptions=-p priority -q msn' -nanopore-raw /project/forage_assemblies/sequence_data/combined_both_vetch.fastq.gz"
```

## Assembly statistics

Obviously, I should have tried this first, but I'm a very impulsive person! Let's try to generate some new statistics in order to compare the quality of the assemblies vs the quality of the reads.

> Ceres: /project/forage_assemblies/analysis/vetch/assembly_stats

### Read set Kmer analysis

```bash
module load jellyfish2/2.2.9

# Making generators so that Jellyfish can use the reads
echo "gunzip -c /project/forage_assemblies/sequence_data/combined_both_vetch.fastq.gz" > combined_vetch_generator.g
echo "gunzip -c /project/forage_assemblies/sequence_data/nobel_vetch_combined_reads.fastq.gz" > nobel_vetch_generator.g
echo "gunzip -c /project/forage_assemblies/sequence_data/vetch_combined_reads.fastq.gz" > dfrc_vetch_generator.g

# Testing permutations of kmers and input reads
for j in 21 32 64; do echo $j; for i in combined_vetch nobel_vetch dfrc_vetch; do echo $i; sbatch --nodes=1 --mem=320000 --ntasks-per-node=20 -p priority -q msn --wrap="jellyfish count -m $j -s 100M -t 20 --bf-size 10G -C -g ${i}_generator.g -o ${i}_reads_${j}.jf"; done; done

# Generating histo files for assemblytics
for i in *.jf; do name=`echo $i | cut -d'.' -f1`; echo $name; sbatch --nodes=1 --mem=15000 --ntasks-per-node=10 -p priority -q msn --wrap="jellyfish histo -t 10 $i > $name.histo"; done

sbatch --nodes=1 --mem=320000 --ntasks-per-node=20 -p priority -q msn --wrap="jellyfish count -m 17 -s 100M -t 20 --bf-size 10G -C -g dfrc_vetch_generator.g -o dfrc_vetch_reads_17.jf"
```

### Assembly-level Kmer analysis

```bash
module load jellyfish2/2.2.9

# Making a generator for Jellyfish
echo "gunzip -c /project/forage_assemblies/assemblies/vetch/vetch_flye/assembly.fasta.gz" > vetch_flye_assembly.g

# Note: flye selected a kmer size of 17
sbatch --nodes=1 --mem=320000 --ntasks-per-node=20 -p priority -q msn --wrap="jellyfish count -m 17 -s 100M -t 20 --bf-size 10G -C -g vetch_flye_assembly.g -o vetch_flye_assembly_17.jf"
```

### Merqury

```bash
module load canu/2.0 r/3.6.1 java_8_sdk/1.8.0_121 bedtools/2.25.0 samtools/1.9

export MERQURY=/project/forage_assemblies/binaries/merqury

sh $MERQURY/best_k.sh 2000000000
	genome: 2000000000
	tolerable collision rate: 0.001
	20.4308

ls /project/forage_assemblies/sequence_data/combined_both_vetch.fastq.gz > input.fofn

# Creating Meryl db for combined both vetch dataset
sbatch --nodes=1 --mem=10000 --ntasks-per-node=1 -p priority -q msn $MERQURY/_submit_build.sh 21 input.fofn combined_both_vetch
```

### Assemblytics

```bash
conda activate /KEEP/rumen_longread_metagenome_assembly/mummer

sbatch --nodes=1 --mem=150000 --ntasks-per-node=2 -p priority -q msn --wrap="nucmer -maxmatch -l 100 -c 500 /project/forage_assemblies/assemblies/vetch/vetch_flye/assembly.fasta /project/forage_assemblies/assemblies/vetch/vetch_limit_flye/00-assembly/draft_assembly.fasta -prefix vetch_vs_vetch"
```

### Purge dups

```bash
module load miniconda/3.6
conda activate /KEEP/rumen_longread_metagenome_assembly/purge_dups

sbatch --nodes=1 --mem=35000 --ntasks-per-node=12 -p priority -q msn --wrap='minimap2 -x map-ont -t 12 /project/forage_assemblies/assemblies/vetch/vetch_flye/assembly.fasta /project/forage_assemblies/sequence_data/vetch_combined_reads.fastq.gz | gzip -c - > dfrc_vetch_readcov.paf.gz; purge_dups/bin/pbcstat dfrc_vetch_readcov.paf.gz; purge_dups/bin/calcuts PB.stat > cutoffs 2>calcults.log'

sbatch --nodes=1 --mem=35000 --ntasks-per-node=12 -p priority -q msn --wrap='purge_dups/bin/split_fa /project/forage_assemblies/assemblies/vetch/vetch_flye/assembly.fasta > dfrc_vetch_flye.split; minimap2 -x asm5 -t 12 -DP dfrc_vetch_flye.split dfrc_vetch_flye.split | gzip -c - > dfrc_vetch_flye.split.self.paf.gz'

sbatch --nodes=1 --mem=35000 --ntasks-per-node=2 --dependency=afterok:2744137:2744138 -p priority -q msn --wrap='purge_dups/bin/purge_dups -2 -T cutoffs -c PB.base.cov dfrc_vetch_flye.split.self.paf.gz > dups.bed 2> purge_dups.log'

sbatch --nodes=1 --mem=35000 --ntasks-per-node=2 --dependency=afterok:2744142 -p priority -q msn --wrap='purge_dups/bin/get_seqs dups.bed /project/forage_assemblies/assemblies/vetch/vetch_flye/assembly.fasta'

perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/assembly_scripts/calculateContigN50.pl -i purged.fa.fai
Total length:   2772992890
Num contigs:    12805
N50 length:     1386686605
N50 value:      447489
L50 value:      1864
Max:    2895845
Min:    16

perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/assembly_scripts/calculateContigN50.pl -i hap.fa.fai
Total length:   846472118
Num contigs:    20106
N50 length:     423302256
N50 value:      72980
L50 value:      3701
Max:    745525
Min:    66


# I am going to run a KAT analysis against both the purged assembly and the 
```

### Busco analysis

```bash
module load busco4

# Necessary configuration
busco_configurator.py /software/7/apps/busco4/4.0.2/config/config.ini /project/forage_assemblies/analysis/vetch/assembly_stats/busco_config.ini
export BUSCO_CONFIG_FILE=/project/forage_assemblies/analysis/vetch/assembly_stats/busco_config.ini

cp -Rp /software/7/apps/augustus/3.3.2/config/ /project/forage_assemblies/analysis/vetch/assembly_stats/AUGUSTUS_CONFIG
export AUGUSTUS_CONFIG_PATH=/project/forage_assemblies/analysis/vetch/assembly_stats/AUGUSTUS_CONFIG

sbatch --nodes=1 --mem=40000 --ntasks-per-node=10 -p priority -q msn --wrap="busco -m genome -i purged.fa -f --offline -l /reference/data/BUSCO/v4/lineages/eudicots_odb10 -c 10 -o purged_vetch"
sbatch --nodes=1 --mem=40000 --ntasks-per-node=10 -p priority -q msn --wrap="busco -m genome -i /project/forage_assemblies/assemblies/vetch/vetch_flye/assembly.fasta -f --offline -l /reference/data/BUSCO/v4/lineages/eudicots_odb10 -c 10 -o dfrc_vetch"
```

## Recalling reads

It turns out that I basecalled all of the vetch reads using the 3.0 version of guppy and they're now up to version 3.6. I just redownloaded guppy and want to give this another go. Let's try it out.

> Ceres: /lustre/project/forage_assemblies/sequence_data/

```bash
# Calling all with the new config High accuracy config file
for i in `ls -d /lustre/project/forage_assemblies/sequence_data/Vetch*/*/fast5_p*`; do echo $i; sbatch --nodes=1 --mem=100000 --ntasks-per-node=25 -p priority -q msn --wrap="/lustre/project/rumen_longread_metagenome_assembly/binaries/ont-guppy-cpu/bin/guppy_basecaller -i $i -s $i/../../g3_6 -c /lustre/project/rumen_longread_metagenome_assembly/binaries/ont-guppy-cpu/data/dna_r9.4.1_450bps_hac.cfg --cpu_threads_per_caller 25"; done

# It was taking a LOONG time, so I requeued a few with specs derived from biostars
for i in B2 B3 B4 B5 B6 B7 B8 B9 Booger MN Vortex Zymo; do name=Vetch${i}; echo $name; sbatch --nodes=1 --mem=100000 --ntasks-per-node=25 -p priority -q msn --wrap="/lustre/project/rumen_longread_metagenome_assembly/binaries/ont-guppy-cpu/bin/guppy_basecaller -i $name/*/fast5_pass -s $name/g3_6 -c /lustre/project/rumen_longread_metagenome_assembly/binaries/ont-guppy-cpu/data/dna_r9.4.1_450bps_hac.cfg --cpu_threads_per_caller 1 --num_callers 25"; done

# testing a new config for VetchB2 (which failed due to memory):
sbatch --nodes=1 --mem=300000 --ntasks-per-node=60 -p priority -q msn --wrap="/lustre/project/rumen_longread_metagenome_assembly/binaries/ont-guppy-cpu/bin/guppy_basecaller -i VetchB2/*/fast5_pass -s VetchB2/g3_6 -c /lustre/project/rumen_longread_metagenome_assembly/binaries/ont-guppy-cpu/data/dna_r9.4.1_450bps_hac.cfg --cpu_threads_per_caller 3 --num_callers 20"

# VetchB7 also failed. Not sure why it's taking so long to call these reads!
# Another new config for that one
sbatch --nodes=1 --mem=300000 --ntasks-per-node=60 -p priority -q msn --wrap="/lustre/project/rumen_longread_metagenome_assembly/binaries/ont-guppy-cpu/bin/guppy_basecaller -i VetchB7/*/fast5_pass -s VetchB7/g3_6 -c /lustre/project/rumen_longread_metagenome_assembly/binaries/ont-guppy-cpu/data/dna_r9.4.1_450bps_hac.cfg --cpu_threads_per_caller 6 --num_callers 10"

# It looks like it takes a long cpu time to call each read and each caller loads up until it runs out of memory. Let's try to redo this
for i in VetchBooger VetchB8 VetchB9 VetchB5 VetchB6 VetchB3; do echo $i; rm -r $i/g3_6; sbatch --nodes=1 --mem=300000 --ntasks-per-node=70 -p priority -q msn --wrap="/lustre/project/rumen_longread_metagenome_assembly/binaries/ont-guppy-cpu/bin/guppy_basecaller -i $i/*/fast5_pass -s $i/g3_6 -c /lustre/project/rumen_longread_metagenome_assembly/binaries/ont-guppy-cpu/data/dna_r9.4.1_450bps_hac.cfg --cpu_threads_per_caller 14 --num_callers 5"; done


# This is still taking a LOOONG time! Let's try the fast calling for the time being?
sbatch --nodes=1 --mem=300000 --ntasks-per-node=70 -p priority -q msn --wrap="/lustre/project/rumen_longread_metagenome_assembly/binaries/ont-guppy-cpu/bin/guppy_basecaller -i VetchB10/*/fast5_pass -s VetchB10/g3_6_f -c /lustre/project/rumen_longread_metagenome_assembly/binaries/ont-guppy-cpu/data/dna_r9.4.1_450bps_fast.cfg --cpu_threads_per_caller 14 --num_callers 5"

for i in VetchB2 VetchB3 VetchB4 VetchB5 VetchB6 VetchB7 VetchB8 VetchB9 VetchBooger VetchMN VetchVortex VetchZymo; do echo $i; sbatch --nodes=1 --mem=300000 --ntasks-per-node=70 -p priority -q msn --wrap="/lustre/project/rumen_longread_metagenome_assembly/binaries/ont-guppy-cpu/bin/guppy_basecaller -i $i/*/fast5_pass -s $i/g3_6_f -c /lustre/project/rumen_longread_metagenome_assembly/binaries/ont-guppy-cpu/data/dna_r9.4.1_450bps_fast.cfg --cpu_threads_per_caller 14 --num_callers 5"; done

# Now, let's run the assembly again with Flye
conda activate /KEEP/rumen_longread_metagenome_assembly/flye

sbatch --nodes=1 --mem=900000 --ntasks-per-node=70 -p priority-mem -q msn-mem flye --nano-raw /lustre/project/forage_assemblies/sequence_data/vetch_36_recalled_reads.fastq -g 2000m -t 70 -o vetch_recalled_flye
```

#### Vetch HIC

> Ceres: /lustre/project/forage_assemblies/sequence_data/vetch_hic

```bash
# queueing up the alignments
module load bwa miniconda/3.6

ls /lustre/project/forage_assemblies/sequence_data/vetch_hic/*.fastq.gz > template.tab

# Starting the alignments
for i in vetch_flye vetch_recalled_flye; do p="/lustre/project/forage_assemblies/assemblies/vetch/"; echo "$p/$i"; sbatch --nodes=1 --ntasks-per-node=2 -p priority -q msn -t 1-0 --wrap="bwa index $p/$i/assembly.fasta; python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b $i -t template.tab -f $p/$i/assembly.fasta -q msn -p priority -e '2-0' -m"; done

# Damn, the script needs readname sorted bams
sbatch --nodes=1 --mem=25000 --ntasks-per-node=10 -p priority -q msn --wrap="samtools sort -n -T vetch_flye/S3HIC/vtemp -o vetch_flye/S3HIC/vetch_flye_readname.bam -@ 10 vetch_flye/S3HIC/vetch_flye.merged.bam"
sbatch --nodes=1 --mem=25000 --ntasks-per-node=10 -p priority -q msn --wrap="samtools sort -n -T vetch_recalled_flye/S3HIC/vtemp -o vetch_recalled_flye/S3HIC/vetch_recalled_flye_readname.bam -@ 10 vetch_recalled_flye/S3HIC/vetch_recalled_flye.merged.bam"

conda activate /lustre/project/forage_assemblies/sequence_data/vetch_hic/hic_qc/hic_qc
python hic_qc/hic_qc.py -b vetch_recalled_flye/S3HIC/vetch_recalled_flye_readname.bam -n 1000000 -r
mv Read_mate_dist* ./vetch_recalled_flye/

python hic_qc/hic_qc.py -b vetch_flye/S3HIC/vetch_flye_readname.bam -n 1000000 -r
```

And now to collate results from the dist.tsv files.

```bash
perl -e '%v = ("perc_pairs_on_same_strand_hq" => ">1.5%", "perc_informative_read_pairs" => ">5.0%", "perc_pairs_intra_hq_gt10kbp" => ">3.0%", "perc_intercontig_pairs_hq_gt10kbp" => ">2.5%", "proximo_usable_rp_hq_per_ctg_gt_5k" => ">600.0", "perc_noninformative_read_pairs" => "<=50%", "perc_zero_dist_pairs" => "<=20%", "perc_mapq0_reads" => "<=20%", "perc_unmapped_reads" => "<=10%", "perc_split_reads" => "<10%", "perc_intercontig_pairs" => "10-60%", "perc_intercontig_pairs_hq" => "10-60%"); chomp(@ARGV); %data; foreach $i (@ARGV){open(IN, "< $i"); while(<IN>){chomp; @s = split(/\t/); if(exists($v{$s[0]})){push(@{$data{$s[0]}}, $s[1]);} } close IN;} print "Category\t" . join("\t", @ARGV) . "\tExpected\n"; foreach $k (keys(%data)){print "$k\t" . join("\t", @{$data{$k}}) . "\t$v{$k}\n";}' vetch_flye_dist.tsv vetch_recalled_dist.tsv
Category        vetch_flye_dist.tsv     vetch_recalled_dist.tsv Expected

proximo_usable_rp_hq_per_ctg_gt_5k      0.94    0.74    >600.0
perc_pairs_on_same_strand_hq    0.67%   0.98%   >1.5%
perc_informative_read_pairs     16.27%  12.90%  >5.0%
perc_noninformative_read_pairs  69.06%  72.89%  <=50%

perc_pairs_intra_hq_gt10kbp     0.38%   0.56%   >3.0%
perc_intercontig_pairs  54.57%  55.84%  10-60%
perc_intercontig_pairs_hq       12.29%  12.77%  10-60%
perc_intercontig_pairs_hq_gt10kbp       12.00%  11.85%  >2.5%

perc_unmapped_reads     1.36%   1.06%   <=10%
perc_zero_dist_pairs    23.15%  21.48%  <=20%
perc_split_reads        38.79%  41.35%  <10%
perc_mapq0_reads        39.51%  46.35%  <=20%
```

#### Illumina QC runs.

Now that the illumina WGS data is here, let's run it through my qc pipeline.

> Ceres: /lustre/project/forage_assemblies/assemblies/vetch/vetch_flye

```bash
module load miniconda/3.6

sbatch --nodes=1 --mem=5000 --ntasks-per-node=2 -p priority -q msn snakemake --cluster-config ~/python_toolchain/snakeMake/assemblyValidation/cluster.json --cluster "sbatch --nodes={cluster.nodes} --ntasks-per-node={cluster.ntasks-per-node} --mem={cluster.mem} --partition={cluster.partition} -q {cluster.qos} -o {cluster.stdout}" -p --jobs 250 -s ~/python_toolchain/snakeMake/assemblyValidation/assemblyValidation --use-conda
```
> Ceres: /lustre/project/forage_assemblies/assemblies/vetch/vetch_recalled_flye

```bash
module load miniconda/3.6

sbatch --nodes=1 --mem=5000 --ntasks-per-node=2 -p priority -q msn snakemake --cluster-config ~/python_toolchain/snakeMake/assemblyValidation/cluster.json --cluster "sbatch --nodes={cluster.nodes} --ntasks-per-node={cluster.ntasks-per-node} --mem={cluster.mem} --partition={cluster.partition} -q {cluster.qos} -o {cluster.stdout}" -p --jobs 250 -s ~/python_toolchain/snakeMake/assemblyValidation/assemblyValidation --use-conda
```