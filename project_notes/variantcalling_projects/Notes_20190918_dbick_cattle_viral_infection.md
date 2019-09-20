# Novel Pathogen detection
---
*9/18/2019*

These are my notes on screening WGS reads for novel viral infections in cattle.

## Table of Contents

## Preparing for the workflow

We need to start generating files before we can make the donuts! Let's prep everything that we need:

> Ceres: /project/rumen_longread_metagenome_assembly/Aspen_Workman

```bash
module load bwa

# Making the alignment spreadsheet
ls *.fastq.gz > fastq_spreadsheet.tab
vim fastq_spreadsheet.tab

# Indexing the reference
sbatch --nodes=1 --ntasks-per-node=2 --mem=15000 -t 1-0 -p msn --wrap="bwa index /project/cattle_genome_assemblies/dominette/symposium_comparison/ARS-UCD1.2_Btau5.0.1Y.fa"

python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b aspen -t fastq_spreadsheet.tab -f /project/cattle_genome_assemblies/dominette/symposium_comparison/ARS-UCD1.2_Btau5.0.1Y.fa -p msn -e '2-0' -m

rm aspen/Aspen/*.sorted.bam*
```

Now, I know that this is a small number of samples (n=1) but let's queue it up with the Snakemake pipeline I created previously.

```bash
mkdir logs
python3 ~/python_toolchain/snakeMake/readScrape/scripts/generateJSONForPipeline.py aspen/Aspen/ > default.json

# I added in my previous centrifuge db and binary location
# I need to modify the cluster.json data to make this work with the time limits on CERES:
cp ~/python_toolchain/snakeMake/readScrape/cluster.json ./
vim cluster.json

# I basically gave everything 2 days max

sbatch --nodes=1 --mem=5000 --ntasks-per-node=2 -p msn -t 2-0 snakemake --cluster-config cluster.json --cluster "sbatch --nodes={cluster.nodes} --ntasks-per-node={cluster.ntasks-per-node} --mem={cluster.mem} --partition={cluster.partition} -o {cluster.stdout} -t {cluster.time}" --jobs 999 -s ~/python_toolchain/snakeMake/readScrape/readScrape

# OK, after debugging a bit, let's count the number of reads that mapped to each scaffold:
samtools sort -T temp -o filtered/Aspen/Aspen.readcontig.sorted.bam filtered/Aspen/Aspen.readcontig.bam
samtools index filtered/Aspen/Aspen.readcontig.sorted.bam
samtools idxstats filtered/Aspen/Aspen.readcontig.sorted.bam
Aspen.scaf.1    5463    8798    0
Aspen.scaf.2    2886    388159  0
Aspen.scaf.3    2189    9271    0
Aspen.scaf.4    1863    515382  0
Aspen.scaf.5    1660    102     0
Aspen.scaf.6    1657    16130   0
Aspen.scaf.7    1610    3530    0
Aspen.scaf.8    1608    4083    0
Aspen.scaf.9    1566    148     0
Aspen.scaf.10   1445    1044195 0
Aspen.scaf.11   1406    60      0
Aspen.scaf.12   1305    741     0
Aspen.scaf.13   1202    790     0
Aspen.scaf.14   1073    407     0
Aspen.scaf.15   1039    220223  0
Aspen.scaf.16   1004    483     0
*       0       0       666068

# Interesting. Let's collate this all together as a single table. 
samtools idxstats filtered/Aspen/Aspen.readcontig.sorted.bam > filtered/Aspen/Aspen.readcontig.sorted.bam.stats

python3 ~/python_toolchain/utils/tabFileLeftJoinTable.py -f filtered/Aspen/Aspen.readcontig.sorted.bam.stats -f classification/Aspen/Aspen.centrifuge.out -c 0 -o aspen_filtered_combine.tab
perl -e 'print "\| Contig \| Size \| Mapped Reads \| seqID \| taxID \| Score \| refLen \| qLen \|\n"; $t = ":-- \|"; print "\| " . $t x 8 . "\n"; while(<>){chomp; @s = split(/\t/); if($s[0] eq "*"){next;} print "\| " . join(" \| ", @s[(0..5, 7, 8)]) . "\|\n";}' < aspen_filtered_combine.tab

perl -e 'print "Contig\tSize\tMapped Reads\tseqID\ttaxID\tTopHit\tScore\trefLen\tqLen\n"; chomp(@ARGV); open(IN, "< $ARGV[0]"); %data; <IN>; while(<IN>){chomp; @s = split(/\t/); $data{$s[1]} = $s[0];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); if($s[0] eq "*"){next;} $n = $data{$s[4]}; print join("\t", @s[(0..4)]) . "\t$n\t$s[5]\t$s[7]\t$s[8]\n";} close IN;' classification/Aspen/Aspen.centrifuge.report aspen_filtered_combine.tab > aspen_filtered_combine.plusTax.tab
```

| Contig | Size | Mapped Reads | seqID | taxID | Name | Score | refLen | qLen |
| :-- |:-- |:-- |:-- |:-- |:-- | :--- |:-- |:-- |
| Aspen.scaf.1 | 5463 | 8798 | NZ_LT703009.1 | 1922226 | Anderseniella sp. Alg231-50 | 20008589 | 5447 | 5463|
| Aspen.scaf.10 | 1445 | 1044195 | genus | 497 | Psychrobacter | 756128 | 950 | 1445|
| Aspen.scaf.11 | 1406 | 60 | NC_000068.7 | 10090 | Mus musculus | 86815 | 427 | 1406|
| Aspen.scaf.12 | 1305 | 741 | species | 1747 | Cutibacterium acnes | 1563076 | 1304 | 1305|
| Aspen.scaf.13 | 1202 | 790 | genus | 1372 | Planococcus | 127312 | 530 | 1202|
| Aspen.scaf.14 | 1073 | 407 | species | 1747 | Cutibacterium acnes | 557505 | 837 | 1073|
| Aspen.scaf.15 | 1039 | 220223 | NC_018708.1 | 358220 | Acidovorax sp. KKS102 | 580644 | 777 | 1039|
| Aspen.scaf.16 | 1004 | 483 | species | 1747 | Cutibacterium acnes | 976144 | 1003 | 1004|
| Aspen.scaf.2 | 2886 | 388159 | NZ_CP015367.1 | 418223 | Methylobacterium phyllosphaerae | 725023 | 2128 | 2886|
| Aspen.scaf.3 | 2189 | 9271 | NZ_CP033200.1 | 1270 | Micrococcus luteus | 541696 | 751 | 2189|
| Aspen.scaf.4 | 1863 | 515382 | species | 47884 | Pseudomonas taetrolens | 647505 | 867 | 1863|
| Aspen.scaf.5 | 1660 | 102 | NC_000067.6 | 10090 | Mus musculus | 50516 | 403 | 1660|
| Aspen.scaf.6 | 1657 | 16130 | NZ_CP015367.1 | 418223 | Methylobacterium phyllosphaerae | 81225 | 300 | 1657|
| Aspen.scaf.7 | 1610 | 3530 | NC_018708.1 | 358220 | Acidovorax sp. KKS102 | 208849 | 472 | 1610|
| Aspen.scaf.8 | 1608 | 4083 | NZ_CP018863.1 | 37928 | Arthrobacter crystallopoietes  | 186624 | 447 | 1608|
| Aspen.scaf.9 | 1566 | 148 | NC_000068.7 | 10090 | Mus musculus | 33033 | 361 | 1566|

```bash
# Now to give a graphical representation of the total data, let's create a Krona plot
perl -lane 'if($F[0] eq "Contig"){next;} print "$F[0]\t$F[4]\t$F[2]";' < aspen_filtered_combine.plusTax.tab > aspen_filtered_combine.prekrona

/project/rumen_longread_metagenome_assembly/binaries/Krona/bin/ktImportTaxonomy -o aspen_filtered_contigs.krona.html aspen_filtered_combine.prekrona
```