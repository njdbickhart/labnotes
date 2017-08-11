# Rumen microbiome technology bias check
---
*8/9/2017*

These are my notes and commands for testing the bias of different sequencing platforms in the collection and determining if there are major limitations to each.

## Mash profile comparisons

I am going to test for kmer cardinality differences in the methods using MASH. My concern is that certain kmers may be underrepresented in the different methods and this is the best method for finding out!

Some oddities of MASH: it doesn't read gzipped files and each separate file is considered an "island" in each sketch comparison. I need to merge files to keep them in the same profile for comparison. The great news is that the sketches I create can be transferred and used on other systems, making them great portable compressed profiles.

#### Fastq merger and preparation

```bash
cat *.fastq > nanopore_yu_and_morrison_3.combined.fq
for i in illuminaR*/*.fastq.gz; do echo $i; gunzip $i; done

# I wrote a python script to filter the reads
python3 ../../../programs_source/python_toolchain/sequenceData/filterNextseqFastqFiles.py -f YMPrepCannula_S1_L001_R1_001.fastq -r YMPrepCannula_S1_L001_R2_001.fastq -o YMRewriteFilterGbase
Processed 29316700 reads and identified 884938 single G-base and 3990662 both G-base artifacts

# Automating the rest
for i in illuminaR1/YMPrepCannula_S1_L002_R1_001.fastq,illuminaR1/YMPrepCannula_S1_L002_R2_001.fastq illuminaR1/YMPrepCannula_S1_L003_R1_001.fastq,illuminaR1/YMPrepCannula_S1_L003_R2_001.fastq illuminaR1/YMPrepCannula_S1_L004_R1_001.fastq,illuminaR1/YMPrepCannula_S1_L004_R2_001.fastq illuminaR2/YMPrepCannula_S1_L001_R1_001.fastq,illuminaR2/YMPrepCannula_S1_L001_R2_001.fastq illuminaR2/YMPrepCannula_S1_L002_R1_001.fastq,illuminaR2/YMPrepCannula_S1_L002_R2_001.fastq illuminaR2/YMPrepCannula_S1_L003_R1_001.fastq,illuminaR2/YMPrepCannula_S1_L003_R2_001.fastq illuminaR2/YMPrepCannula_S1_L004_R1_001.fastq,illuminaR2/YMPrepCannula_S1_L004_R2_001.fastq; do file1=`echo $i | cut -d',' -f1`; file2=`echo $i | cut -d',' -f2`; folder=`echo $i | cut -d'/' -f1`; echo $folder; lane=`basename $file1 | cut -d'_' -f3`; python3 ../../programs_source/python_toolchain/sequenceData/filterNextseqFastqFiles.py -f $file1 -r $file2 -o ${folder}"/YMRewriteFilterGbase."${lane} ; done
illuminaR1
Processed 52887277 reads and identified 1863215 single G-base and 26864810 both G-base artifacts
illuminaR1
Processed 26327191 reads and identified 572981 single G-base and 1135473 both G-base artifacts
illuminaR1
Processed 28095512 reads and identified 793048 single G-base and 3164439 both G-base artifacts
illuminaR2
Processed 63812154 reads and identified 808176 single G-base and 34846 both G-base artifacts
illuminaR2
Processed 63221220 reads and identified 892873 single G-base and 85939 both G-base artifacts
illuminaR2
Processed 63346742 reads and identified 777543 single G-base and 17110 both G-base artifacts
illuminaR2
Processed 62625319 reads and identified 818123 single G-base and 34208 both G-base artifacts
```

#### Mash sketch generation

I want to create Mash sketches that conform with Serge's metagenomics estimates and are comparable to the sketches I've generated [previously](https://github.com/njdbickhart/labnotes/blob/master/project_notes/metagenomics/Notes_20161219_dbick_metagenomics_software_test_notes.md#mash). So I will use the same settings as before.

```bash
for i in *.fastq; do name=`echo $i | cut -d'.' -f1`; echo $name; mash sketch -p 3 -k 21 -s 10000 -r -m 2 -o $name $i; done

for i in cheryl/*/*.subreads.bam; do name=`basename $i | cut -d'.' -f1`; echo $name; samtools fastq $i > $name.fq; mash sketch -p 3 -k 21 -s 10000 -r -m 2 -o $name $name.fq; rm $name.fq; done

# Sketching all Illumina reads from stdin
for i in illuminaR*/YMRewriteFilterGbase*fq; do cat $i; done | mash sketch -p 3 -k 21 -s 10000 -r -m 2 -o illumina_filtered_non-interleaved -

for i in *.bam; do name=`echo $i | cut -d'.' -f1`; echo $name; samtools fastq $i | mash sketch -p 3 -k 21 -s 10000 -r -m 2 -o $name - ; done

for i in ./*.bam; do name=`basename $i | cut -d'.' -f1`; echo $name; samtools fastq $i | mash sketch -p 3 -k 21 -s 10000 -r -m 2 -o $name - ; done

for i in `find ./nanopore -name "*.fastq"`; do cat $i ; done | mash sketch -p 3 -k 21 -s 10000 -r -m 2 -o nanopore_yu_morrison_fastq -