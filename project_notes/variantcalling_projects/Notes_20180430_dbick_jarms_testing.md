# Testing JaRMS
---
*4/30/2018*

These are my notes on the testing and simulation of JaRMs for future deployment.

## Table of Contents


## Generating simulation output

I have some scripts that I can use to make high throughput simulations of sequence data for testing with JaRMs. To keep the test files small, I will generate sequence from the first 10 megabases of chr29 and pepper that with some example regions for the algorithm.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/jarms_testing

```bash
samtools faidx /mnt/nfs/nfs2/Genomes/umd3_kary_unmask_ngap.fa chr28:1-10000000 > chr28_first_tenmb.fa
samtools faidx /mnt/nfs/nfs2/Genomes/umd3_kary_unmask_ngap.fa chr29:1-10000000 > chr29_first_tenmb.fa
cat chr28_first_tenmb.fa chr29_first_tenmb.fa > test_chrs_ref.fa

# Let's make a list of the non-gap regions first.
java -jar /mnt/nfs/nfs2/bickhart-users/binaries/GetMaskBedFasta/store/GetMaskBedFasta.jar -f test_chrs_ref.fa -o test_chrs_ref.gaps.bed -s test_chrs_ref.gaps.stats

bedtools subtract -a chrlens.bed -b test_chrs_ref.gaps.bed > test_chrs_ref.nongap.bed

# Now I will generate a template input file with simulated regions for assessment
perl generateTANDAndDelBedRegions.pl test_chrs_ref.nongap.bed test_chrs_ref.template.bed 10

# I want several regions that are comparable across samples to see if the INI model accurrately detects them
mkdir sim_bams
for i in `seq 1 6`; do echo $i; cp test_chrs_ref.template.bed ./sim_bams/sample_${i}_sites.bed; done

# Let's make a list of sites from the existing template
```

chr | start | end | mut | samples
:---| ---: | ---: | :---| :----
chr28|524975|567748|DEL| all
chr28|793775|930573|DEL| 1, 3
chr28|3409967|3437802|TAND| 2
chr28|3460971|3557628|TAND| 5,6
chr28|4264645|4458163|TAND|all
chr28|5754604|5755850|TAND| 1
chr28|8690624|8808241|DEL| 6
chr28|9234676|9252751|TAND| 3
chr29|2631263|2951482|DEL|	2
chr29|4283939|4472931|DEL| 4
chr29|5544899|5558117|DEL| 5
chr29|6084195|6086996|DEL| all
chr29|9210317|9211109|TAND| all


```bash
# Now to generate all of the fasta files and then create the reads
for i in `seq 1 6`; do samp="sample_"$i; echo $samp; sbatch generate_test_bams.sh $samp sim_bams/${samp}_sites.bed test_chrs_ref.fa; done

```