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

# 