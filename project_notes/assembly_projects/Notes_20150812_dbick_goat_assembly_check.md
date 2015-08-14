# Identifying problems with PacBio reference assemblies
---
*8/12/2015*

My goal is to identify any putative INDELs that missed quiver consensus.

> Blade14: /mnt/nfs/nfs2/GoatData/Ilmn

```bash
# Calling SNPs and INDELs with samtools by default
samtools mpileup -ugf ../Goat-Genome-Assembly/Papadum-v3/papadum-v3s.ctg.fa Goat250-PBv3-ctg.bam | bcftools call -vmO z -o Goat250-PBv3-ctg.vcf.gz

# Just checking initial stats
gunzip -c Goat250-PBv3-ctg.vcf.gz | grep 'INDEL' | wc -l
	3,357,559

# A test case
samtools faidx ../Goat-Genome-Assembly/Papadum-v3/papadum-v3s.ctg.fa utg42009:2337272-2337278
	>utg42009:2337272-2337278
	AGGGGGC

# The illumina reads show that the region should have 6 G's, not 5
utg42009        2337272 .       AGGGGG  AGGGGGG 89      .       INDEL;IDV=8;IMF=0.8;DP=10;VDB=0.0849961;SGB=-0.651104;MQSB=1;MQ0F=0;AC=2;AN=2;DP4=0,0,2,6;MQ=60   GT:PL   1/1:116,24,0
```


--
*8/13/2015*

OK, let's filter this out to get the best INDEL calls so that they can be corrected in the assembly. Since there are over 3 million calls, let's profile what types of INDELs there are in a script. Here are the types that could be present:

* simple INS homopolymeric
* simple DEL homopolymeric
* simple INS 
* simple DEL
* complex INS
* complex DEL


In my above terms, "simple" means a single base extension. I do not expect many "complex" homopolymeric INDELs, but if there is a large amount of "complex" variants, I might just pull them out as a separate functional category. Samtools HTSLIB also has better VCF info fields for INDELs. Here are ones that I will attempt to use to filter away putative bad calls:

|Field | Type | Description|
|:-- | :-- | :--|
|IDV | Integer | Maximum number of reads supporting an indel |
| IMF | Float | Maximum fraction of reads supporting an indel |
| DP | Integer | Raw read depth |

I will also filter away any "heterozygous" INDELs as these might just be genuine assembly fork sites.

I created a script that automates data collection and filters the VCF simultaneously

> Blade14: /mnt/nfs/nfs2/GoatData/Ilmn

```bash
# stats are written to stdout
perl ~/perl_toolchain/assembly_scripts/pacbioINDELCompFilter.pl -v Goat250-PBv3-ctg.vcf.gz -o Goat250-PBv3-ctg.indel.filtered.vcf > Goat250-PBv3-ctg.indel.filtered.vcf.stats

```

Oops! The file got corrupted! I had to recall using my snp pipeline

```bash
perl ~/perl_toolchain/sequence_data_scripts/samtoolsSNPFork.pl -r ../Goat-Genome-Assembly/Papadum-v3/papadum-v3s.ctg.fa -i Goat250-PBv3-ctg.bam -o Goat250-PBv3-ctg -n 20 -t 1

bcftools view -o Goat250-PBv3-ctg.samtools.merged.derek.vcf.gz -O z Goat250-PBv3-ctg.samtools.merged.bcf

perl ~/perl_toolchain/assembly_scripts/pacbioINDELCompFilter.pl -v Goat250-PBv3-ctg.samtools.merged.derek.vcf.gz -o Goat250-PBv3-ctg.indel.filtered.vcf > Goat250-PBv3-ctg.indel.filtered.vcf.stats

```

*8/14/2015*

--

The script finished and now the stats look like the following:

|Category|Count|
|:---|---:|
Simple_INS   |   179,251
Simple_DEL    |  1939
Simple_HOM_INS|  2,779,135
Simple_HOM_DEL|  7653
Complex_INS   |  48,844
Complex_DEL   |  7274
**Total_INDELS**  |  **3,024,096**

Correcting the simple homopolymeric insertions would reduce the INDEL count down to 244,961, with 73% of those remaining INDELs being simple, one-base insertions compared to the Illumina data.