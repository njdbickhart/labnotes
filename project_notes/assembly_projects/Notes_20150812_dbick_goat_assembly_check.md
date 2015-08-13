# Identifying problems with PacBio reference assemblies
---
*8/12/2015*

My goal is to identify any putative INDELs that missed quiver consensus.

> Blade14: /mnt/nfs/nfs2/GoatData/Ilmn

```bash
# Calling SNPs and INDELs with samtools by default
samtools mpileup -ugf ../Goat-Genome-Assembly/Papadum-v3/papadum-v3s.ctg.fa Goat250-PBv3-ctg.bam | bcftools call -vmO z -o Goat250-PBv3-ctg.vcf.gz

```
