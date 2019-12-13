# MetaFlye analysis of sheep dataset
---
*12/5/2019*

These are my notes and plans for validating the MetaFlye data for use in the reviewer rebuttals.

## Table of Contents


## Analysis plan

I will need to divy up tasks and run them in order to generate the results needed.

* Run CheckM on both assemblies
* Run BlobTools on both assemblies
	* Taxify contigs and derive contig taxonomic profiles
	* Derive a consensus taxonomic affiliation per bin and plot CheckM results
* Run Minimap2 using CCS reads on both assemblies
* Align Hi-C reads to both assemblies
	* Run the viral association script
	* Count number of split-read or hard-clipped CCS reads
	* Adapt DESMAN to pick strains from the alignments
* Run prodigal on the contigs


## Running the analysis

#### CCS read alignments

> Ceres: /project/forage_assemblies/sheep_project

```bash
module load minimap2/2.6
for i in *.tar.gz; do echo $i; sbatch --nodes=1 --mem=1000 --ntasks-per-node=1 -p msn --wrap="tar -xvf $i"; done

# OK, so the contigs are not properly binned, but given their sizes, this is probably OK.
# Let me start the long-term alignment first so that I can analyze this next week
# I will use this for blobtools "coverage" estimates as well
for i in canu.contigs.fasta flye_contigs.fasta; do echo $i; sbatch --nodes=1 --mem=30000 --ntasks-per-node=3 -p msn --wrap="minimap2 -x map-pb $i /project/rumen_longread_metagenome_assembly/sheep_poop/sheep_poop_CCS.fastq.gz > $i.ccs.paf"; done

# The blobtools coverage estimates didn't work very well with the CCS reads in sam format. Maybe I can resolve this with scripting?
```

#### WGS read alignments

```bash
# BWA indexing
for i in canu.contigs.fasta flye_contigs.fasta ; do echo $i; sbatch --nodes=1 --mem=15000 --ntasks-per-node=1 -p msn -q msn --wrap="module load bwa; bwa index $i"; done

# queueing alignment jobs
python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -f /project/forage_assemblies/sheep_project/flye_contigs.fasta -t /project/rumen_longread_metagenome_assembly/sheep_poop/hic_links.tab -b flye_hic -p short -q memlimit -m
python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -f /project/forage_assemblies/sheep_project/flye_contigs.fasta -t /project/rumen_longread_metagenome_assembly/sheep_poop/wgs_reads.tab -b flye_wgs -p short -q memlimit -m

python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -f /project/forage_assemblies/sheep_project/canu.contigs.fasta -t /project/rumen_longread_metagenome_assembly/sheep_poop/hic_links.tab -b canu_hic -p short -q memlimit -m
python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -f /project/forage_assemblies/sheep_project/canu.contigs.fasta -t /project/rumen_longread_metagenome_assembly/sheep_poop/wgs_reads.tab -b canu_wgs -p short -q memlimit -m
```

#### ORF annotation

I am going to run both assemblies through Prodigal. The ORF data will be extremely useful downstream.

```bash
module load prodigalorffinder/2.6.3

for i in canu.contigs.fasta flye_contigs.fasta; do echo $i; sbatch --nodes=1 --mem=100000 --ntasks-per-node=2 -p msn --wrap="prodigal -a $i.prod.prottrans -c -d $i.prod.genenuc -f gff -i $i -o $i.prod.out -p meta"; done
```

#### Blobtools

```bash
module load miniconda
conda activate /KEEP/rumen_longread_metagenome_assembly/blobtools

# Diamond run before blobtools
for i in canu.contigs.fasta flye_contigs.fasta; do echo $i; sbatch -t 2-0 -p msn -q msn --nodes=1 --ntasks-per-node=30 --mem=100000 --wrap="diamond blastx --query $i --db /project/rumen_longread_metagenome_assembly/assemblies/protists/uniprot_ref_proteomes.diamond.dmnd --threads 29 --outfmt 6 --sensitive --max-target-seqs 1 --evalue 1e-25 -o $i.diamondout.tsv"; done

# NOTE: I had to clone the blobtools repo because Conda screwed up the libraries for the program
# Blobtools taxify
for i in canu.contigs.fasta flye_contigs.fasta; do echo $i; sbatch -t 2-0 --nodes=1 --mem=20000 --ntasks-per-node=3 -p msn -q msn --wrap="./blobtools/blobtools taxify -f $i.diamondout.tsv -m /project/rumen_longread_metagenome_assembly/assemblies/protists/uniprot_ref_proteomes.taxids -s 0 -t 2 -o ${i}_unip"; done

# Copying the database file
cp /project/rumen_longread_metagenome_assembly/assemblies/protists/blob_ncbi.db ./

# Blobtools create
for i in canu flye; do sbatch -t 2-0 --nodes=1 --mem=10000 --ntasks-per-node=1 -p msn -q msn --wrap="./blobtools/blobtools create -i $i.contigs.fasta -b ${i}_wgs/63/63.sorted.merged.bam -t $i.contigs.fasta_unip.$i.contigs.fasta.diamondout.tsv.taxified.out -o $i.blobtools --db blob_ncbi.db"; done

# Blobtools view and plot
for i in canu flye; do sbatch -t 2-0 --nodes=1 --mem=10000 --ntasks-per-node=1 -p msn -q msn --wrap="./blobtools/blobtools plot -i $i.blobtools.blobDB.json --notitle -r superkingdom -o ${i}_supkingdom"; sbatch -t 2-0 --nodes=1 --mem=10000 --ntasks-per-node=1 -p msn -q msn --wrap="./blobtools/blobtools plot -i $i.blobtools.blobDB.json --notitle -r phylum -o ${i}_phylum"; sbatch -t 2-0 --nodes=1 --mem=10000 --ntasks-per-node=1 -p msn -q msn --wrap="./blobtools/blobtools view -i $i.blobtools.blobDB.json -o ${i}_table -r all"; done

mkdir blob_plots
mv *.stats.txt ./blob_plots/
mv *.png ./blob_plots/

# Taking stock of superkingdoms
python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f canu_table.canu.blobtools.blobDB.table.txt -c 5 -i '#' -d '\t' -m
|Entry     | Value|
|:---------|-----:|
|Bacteria  | 18821|
|Archaea   |   881|
|Eukaryota |   296|
|Viruses   |   208|
|no-hit    |    71|

python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f flye_table.flye.blobtools.blobDB.table.txt -c 5 -i '#' -d '\t' -m
|Entry     | Value|
|:---------|-----:|
|Bacteria  | 15751|
|Archaea   |   584|
|no-hit    |   249|
|Eukaryota |   205|
|Viruses   |   123|
```

#### Network plot and read alignment

I added a rudimentary network plot to my pipeline script. Let's see if it works!

```bash
source activate /KEEP/rumen_longread_metagenome_assembly/seaborn/
module load samtools/1.9 minimap2/2.6

perl -lane 'if($F[0] =~ /^#/){next;} if($F[5] eq "Viruses"){print "$F[0]\t$F[1]";}' < canu_table.canu.blobtools.blobDB.table.txt > canu.viral.contigs.list
perl -lane 'if($F[0] =~ /^#/){next;} if($F[5] eq "Viruses"){print "$F[0]\t$F[1]";}' < flye_table.flye.blobtools.blobDB.table.txt > flye.viral.contigs.list

# Preparing list of viral contig fastas
perl -lane 'print $F[0];' < canu.viral.contigs.list | xargs -I {} samtools faidx canu.contigs.fasta {} >> canu.viral.contigs.fa
perl -lane 'print $F[0];' < flye.viral.contigs.list | xargs -I {} samtools faidx flye.contigs.fasta {} >> flye.viral.contigs.fa

sbatch --nodes=1 --mem=30000 --ntasks-per-node=4 -p msn -q msn -J canuflye --wrap="python3 ~/rumen_longread_metagenome_assembly/binaries/RumenLongReadASM/viralAssociationPipeline.py -a canu.contigs.fasta -g canu.viral.contigs.fa -b canu_table.canu.blobtools.blobDB.table.txt -i canu_hic/63/63.sorted.merged.bam -v canu.viral.contigs.list -l /project/rumen_longread_metagenome_assembly/sheep_poop/sheep_poop_CCS.fastq -m /software/7/apps/minimap2/2.6/minimap2 -o canu.contigs.vassoc"

sbatch --nodes=1 --mem=30000 --ntasks-per-node=4 -p msn -q msn -J flyeflye --wrap="python3 ~/rumen_longread_metagenome_assembly/binaries/RumenLongReadASM/viralAssociationPipeline.py -a flye.contigs.fasta -g flye.viral.contigs.fa -b flye_table.flye.blobtools.blobDB.table.txt -i flye_hic/63/63.sorted.merged.bam -v flye.viral.contigs.list -l /project/rumen_longread_metagenome_assembly/sheep_poop/sheep_poop_CCS.fastq -m /software/7/apps/minimap2/2.6/minimap2 -o flye.contigs.vassoc"
```

### TODO:
### Run DESMAN
### Run a script to calculate Hi-C intercontig mapping rates
### Run a script to calculate number of chimeric CCS read mappings