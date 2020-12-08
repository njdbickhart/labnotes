# Nelore ReadScrape pipeline
---
*11/2/2020*

## Table of Contents


## Preparing files and running pipeline

OK, Yuri couldn't send me the bams but he could send me the readscrape partial files. I will be adjusting the pipeline to work from there.

> Ceres: /lustre/project/rumen_longread_metagenome_assembly/analysis/yuri_readscrape

```bash
module load miniconda/3.6
for i in *.gz; do echo $i; gunzip $i; done

mkdir rawreads; mkdir fakes;
for i in NEL* Sao*; do name=`echo $i | cut -d'_' -f1`; echo $name; done | sort | uniq| perl -ne 'chomp; print "$_ ";'; echo

for i in NEL1 NEL10 NEL11 NEL12 NEL13 NEL14 NEL15 NEL16 NEL17 NEL18 NEL19 NEL2 NEL20 NEL21 NEL22 NEL23 NEL24 NEL3 NEL4 NEL5 NEL6 NEL7 NEL8 NEL9 Sao8131 Sao8243 Sao8254 Sao8344; do mkdir rawreads/$i; touch fakes/$i.bam; echo "  \"$i\" : \"/lustre/project/rumen_longread_metagenome_assembly/analysis/yuri_readscrape/fakes/${i}.bam\","; mv ${i}_* rawreads/$i/; done

## This was a redo because snakemake deleted all of my files!
for i in NEL1 NEL10 NEL11 NEL12 NEL13 NEL14 NEL15 NEL16 NEL17 NEL18 NEL19 NEL2 NEL20 NEL21 NEL22 NEL23 NEL24 NEL3 NEL4 NEL5 NEL6 NEL7 NEL8 NEL9 Sao8131 Sao8243 Sao8254 Sao8344; do echo "  \"$i\" : \"/lustre/project/rumen_longread_metagenome_assembly/analysis/yuri_readscrape/fakes/${i}.bam\","; cp ${i}_* rawreads/$i/; gunzip rawreads/$i/*; done

module load bedtools/2.25.0 bwa samtools centrifuge cd-hit
conda activate /KEEP/rumen_longread_metagenome_assembly/snakemake
sbatch --nodes=1 --mem=5000 --ntasks-per-node=2 -p priority -q msn snakemake --cluster-config readScrape/cluster.json --cluster "sbatch --nodes={cluster.nodes} --ntasks-per-node={cluster.ntasks-per-node} --mem={cluster.mem} --partition=priority -q msn -o {cluster.stdout}" --jobs 999 -s readScrape/readScrape
```

## Screening assembled regions

OK, it finished! Let's see if we have any clusters that align with the ASIP region. ASIP will be defined as 35kb upstream and downstream of the gene itself. These are the approximate coordinates: chr13:63,608,465-63,692,383

```bash
# I had to extend the region to find something!
for i in genotypes/*.genotype.tab; do perl -lane 'if($F[0] eq "sample"){next;}elsif($F[5] eq "13" && $F[6] < 69692383 && $F[7] > 60608465){print "$F[1]";}' < $i; done | python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f stdin -c 0 -m
|Entry       | Value|
|:-----------|-----:|
|Cluster8857 |     2|
|Cluster6188 |     1|
...

grep Cluster8857 genotypes/genotype_summary.tab
Cluster8857     NEL23.scaf.740  1054    1035.578947368421       12.394038420371745      19      99.89105263157896       99.23   NEL1;NEL10;NEL11;NEL15;NEL2;NEL20;NEL21;NEL22;NEL23;NEL24;NEL3;NEL5;NEL6;NEL7;NEL8;NEL9;Sao8243;Sao8254;Sao8344

grep Cluster8857 genotypes/association_counts.tab
Cluster8857     9       X;5     3;4     13;3    2;2     25;1    1;1     4;1     20;1    11;1


grep Cluster8857 genotypes/*.genotype.tab
genotypes/NEL1.genotype.tab:NEL1        Cluster8857     NEL1.scaf.877   1045    100.0   25      19195115        19195115        4
genotypes/NEL10.genotype.tab:NEL10      Cluster8857     NEL10.scaf.816  1028    99.9    13      62082682        62082682        2
genotypes/NEL11.genotype.tab:NEL11      Cluster8857     NEL11.scaf.917  1025    99.9    2       8250944 8250945 4
genotypes/NEL15.genotype.tab:NEL15      Cluster8857     NEL15.scaf.772  1021    99.9    X       2238872 2238918 5
genotypes/NEL2.genotype.tab:NEL2        Cluster8857     NEL2.scaf.956   1031    100.0   X       2238872 2238920 3
genotypes/NEL20.genotype.tab:NEL20      Cluster8857     NEL20.scaf.4824 1045    99.23   X       2238907 2238922 3
genotypes/NEL21.genotype.tab:NEL21      Cluster8857     NEL21.scaf.797  1036    99.9    1       69528537        69528537        1
genotypes/NEL22.genotype.tab:NEL22      Cluster8857     NEL22.scaf.696  1048    100.0   3       65802718        65802718        2
genotypes/NEL23.genotype.tab:NEL23      Cluster8857     NEL23.scaf.740  1054    100.0   2       8250945 8250945 1
genotypes/NEL24.genotype.tab:NEL24      Cluster8857     NEL24.scaf.894  1028    99.9    3       65802718        65802718        3
genotypes/NEL3.genotype.tab:NEL3        Cluster8857     NEL3.scaf.968   1039    99.9    4       25561998        25561998        1
genotypes/NEL5.genotype.tab:NEL5        Cluster8857     NEL5.scaf.956   1043    100.0   20      29317941        29317950        2
genotypes/NEL6.genotype.tab:NEL6        Cluster8857     NEL6.scaf.905   1049    100.0   X       2238872 2238872 3
genotypes/NEL7.genotype.tab:NEL7        Cluster8857     NEL7.scaf.967   1010    99.9    13      49456693        49456693        1
genotypes/NEL8.genotype.tab:NEL8        Cluster8857     NEL8.scaf.857   1022    99.8    X       2238879 2238922 5
genotypes/NEL9.genotype.tab:NEL9        Cluster8857     NEL9.scaf.856   1035    99.9    13      62082682        62082683        3
genotypes/Sao8243.genotype.tab:Sao8243  Cluster8857     Sao8243.scaf.827        1018    99.9    11      89852525        89852525        1
genotypes/Sao8254.genotype.tab:Sao8254  Cluster8857     Sao8254.scaf.814        1051    99.9    3       65802718        65802718        2
genotypes/Sao8344.genotype.tab:Sao8344  Cluster8857     Sao8344.scaf.735        1048    99.9    3       65802718        65802718        4
```

## Grabbing the paragon for blasting

Now all that remains to do is to pull the sequence for the best cluster (**NEL23.scaf.740**) and blast it. Maybe there will be an exon in there?

```bash
samtools faidx assembly/NEL23/NEL23.renamed.fasta NEL23.scaf.740

```

Blast hits were found to Yak and sheep, but no annotation to compare with. 