# Brangus haploid assembly analysis
---
*11/3/2017*

These are my notes on the analysis of the haploid Brangus reference assembly. Paternal haplotypes are derived from an Angus and the maternal is from a Brahman.

## Table of Contents


## Assembly download, prep and marker analysis

I will download the assembly, map the recmap markers to it and then try to estimate the divergence of the assembly from the recmap.

> Assembler 2: /mnt/nfs/nfs2/bickhart-users/brangus_asm 

```bash

# Queueing analysis on the paternal haplotype assembly
sbatch --mem=20000 --ntasks-per-node=1 --nodes=1 --wrap="module load samtools bwa java; wget https://gembox.cbcb.umd.edu/seqdata/bos_taurus/asm/paternal.arrow.fasta; samtools faidx paternal.arrow.fasta; bwa index paternal.arrow.fasta; perl /mnt/nfs/nfs2/dbickhart/dominette_asm/recombination/alignAndOrderSnpProbes.pl -a paternal.arrow.fasta -p /mnt/nfs/nfs2/dbickhart/dominette_asm/recombination/rcmap_manifest_correct.sorted.fa -o paternal.arrow.recmap; bwa mem paternal.arrow.fasta /mnt/nfs/nfs2/dbickhart/dominette_asm/recombination/rcmap_manifest_correct.sorted.fa paternal.arrow.recmap.sam; touch paternal.arrow.recmap.dummy.jf; java -jar /mnt/nfs/nfs2/bickhart-users/cattle_asms/assembly_revision/CombineFasta.jar misassembly -s paternal.arrow.recmap.sam -j paternal.arrow.recmap.dummy.jf -f paternal.arrow.fasta -o paternal.arrow.misassembly"

# Queueing the analysis on the maternal haplotype assembly
sbatch --mem=20000 --ntasks-per-node=1 --nodes=1 --wrap="module load samtools bwa java; wget https://gembox.cbcb.umd.edu/seqdata/bos_taurus/asm/maternal.arrow.fasta; samtools faidx maternal.arrow.fasta; bwa index maternal.arrow.fasta; perl /mnt/nfs/nfs2/dbickhart/dominette_asm/recombination/alignAndOrderSnpProbes.pl -a maternal.arrow.fasta -p /mnt/nfs/nfs2/dbickhart/dominette_asm/recombination/rcmap_manifest_correct.sorted.fa -o maternal.arrow.recmap; bwa mem maternal.arrow.fasta /mnt/nfs/nfs2/dbickhart/dominette_asm/recombination/rcmap_manifest_correct.sorted.fa maternal.arrow.recmap.sam; touch maternal.arrow.recmap.dummy.jf; java -jar /mnt/nfs/nfs2/bickhart-users/cattle_asms/assembly_revision/CombineFasta.jar misassembly -s maternal.arrow.recmap.sam -j maternal.arrow.recmap.dummy.jf -f maternal.arrow.fasta -o maternal.arrow.misassembly"

# And I also have the Angus Brahman falcon unzip
bwa index angusBrahmanF1_FALCONUnzip_arrow.fa

# Queueing up the HD probes
module load bwa; for i in angusBrahmanF1_FALCONUnzip_arrow.fa maternal.arrow.fasta paternal.arrow.fasta; do echo $i; name=`echo $i | cut -d'.' -f1`; perl /mnt/nfs/nfs2/dbickhart/dominette_asm/recombination/alignAndOrderSnpProbes.pl -a $i -p /mnt/nfs/nfs2/dbickhart/dominette_asm/recombination/BovineHD_B1.probseq.rev1coords.fa -o ${name}.HDprobes ; done

# I also screwed up the CombineFasta execution. Rerunning
bwa mem paternal.arrow.fasta /mnt/nfs/nfs2/dbickhart/dominette_asm/recombination/rcmap_manifest_correct.sorted.fa > paternal.arrow.recmap.sam
bwa mem maternal.arrow.fasta /mnt/nfs/nfs2/dbickhart/dominette_asm/recombination/rcmap_manifest_correct.sorted.fa > maternal.arrow.recmap.sam

java -jar /mnt/nfs/nfs2/bickhart-users/cattle_asms/assembly_revision/CombineFasta.jar missassembly -s paternal.arrow.recmap.sam -j paternal.arrow.recmap.dummy.jf -f paternal.arrow.fasta -o paternal.arrow.misassembly
java -jar /mnt/nfs/nfs2/bickhart-users/cattle_asms/assembly_revision/CombineFasta.jar missassembly -s maternal.arrow.recmap.sam -j maternal.arrow.recmap.dummy.jf -f maternal.arrow.fasta -o maternal.arrow.misassembly

# Grabbing the unmapped reads for comparisons
perl -lane 'if($F[1] eq "*"){print "$F[0]";}' < maternal.HDprobes.tab > maternal.HDprobes.missing.list
perl -lane 'if($F[1] eq "*"){print "$F[0]";}' < paternal.HDprobes.tab > paternal.HDprobes.missing.list
perl -lane 'if($F[1] eq "*"){print "$F[0]";}' < angusBrahmanF1_FALCONUnzip_arrow.HDprobes.tab > angusBrahmanF1_FALCONUnzip_arrow.HDprobes.missing.list

perl ~/sperl/bed_cnv_fig_table_pipeline/nameListVennCount.pl angusBrahmanF1_FALCONUnzip_arrow.HDprobes.missing.list maternal.HDprobes.missing.list paternal.HDprobes.missing.list
File Number 1: angusBrahmanF1_FALCONUnzip_arrow.HDprobes.missing.list
File Number 2: maternal.HDprobes.missing.list
File Number 3: paternal.HDprobes.missing.list
Set     Count
1       67
1;2     35
1;2;3   354
1;3     4
2       3846
2;3     16
3       2793

# Interesting! So, the majority of the probes are unique to each
perl -lane 'if($F[1] eq "*"){print "$F[0]";}' < ../cattle_asms/assembly_revision/ARS-UCD1.0.18.HDprobes.tab > ARS-UCD1.0.18.HDprobes.missing.list
perl ~/sperl/bed_cnv_fig_table_pipeline/nameListVennCount.pl angusBrahmanF1_FALCONUnzip_arrow.HDprobes.missing.list maternal.HDprobes.missing.list paternal.HDprobes.missing.list ARS-UCD1.0.18.HDprobes.missing.list
File Number 1: angusBrahmanF1_FALCONUnzip_arrow.HDprobes.missing.list
File Number 2: maternal.HDprobes.missing.list
File Number 3: paternal.HDprobes.missing.list
File Number 4: ARS-UCD1.0.18.HDprobes.missing.list
Set     Count
1       56
1;2     9
1;2;3   266
1;2;3;4 88
1;2;4   26
1;3     1
1;3;4   3
1;4     11
2       3421
2;3     14
2;3;4   2
2;4     425
3       2714
3;4     79
4       291
```

Now to generate a venn for this.

```R
library(VennDiagram)
pdf(file="brangus_missing_HD_probes.pdf", useDingbats=FALSE)
venn <- draw.triple.venn(area1 = 460, area2 = 4251, area3 = 3167, n12 = 389, n23 = 370, n13 = 358, n123 = 354, category = c("angusBrahmanF1_FALCONUnzip", "maternal", "paternal"), fill = c("red", "green", "blue"), cex = 2, cat.cex = 2, cat.col = c("red", "green", "blue"))
dev.off()

# And for the quad comparison
pdf(file="brangus_missing_HD_probes_dominettecomp.pdf", useDingbats=FALSE)
venn <- draw.quad.venn(area1 = 460, area2 = 4251, area3 = 3167, area4 = 925, n12= 389, n13 = 358, n14 = 128, n23 = 370, n24 = 541, n34 = 172, n123 = 354, n124 = 114, n134 = 91, n234 = 90, n1234 = 88, category = c("F1_FALCONUnzip", "maternal", "paternal", "Dominettev18"), fill = c("red", "green", "blue", "orange"), cex = 2, cat.cex = 2, cat.col = c("red", "green", "blue", "orange"))
dev.off()
```
Checking the intersection of PAR snps

```bash
dos2unix par_HD_snps.list
