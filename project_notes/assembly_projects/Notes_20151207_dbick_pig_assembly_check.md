# Pig assembly validation work
---
*12/07/2015*

These are my notes on identifying misassemblies on the pig assembly for Tim Smith and Alan Archibald.

Let's download the data from the SRA first -- I'm going to have to use an unrelated pig in order to get this to work.

OK, I tried the hard way (using NCBI's SRA tools), but the "virtual database" error was bugging me. Instead, I downloaded from their ftp site using wget:

> /mnt/iscsi/vnx_gliu_7/pig

```bash
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/ERR/ERR977/ERR977386/ERR977386.sra
```

All SRA accessions follow the same scheme on the NCBI ftp site:

> ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/{SRR|ERR|DRR}/{first 6 characters of accession}/{accession}/{accession}.sra

```bash
bwa mem -t 10 -M /mnt/nfs/nfs2/dbickhart/pig/pig_03Aug2015_1Wx8C.fasta ERR977386.fastq.gz > pig_03Aug2015_ERR977386.sam
```

## Generating a repeatmasked reference

```bash
# I had to rename the scaffold due to its incredible length!
/mnt/nfs/nfs2/bickhart-users/binaries/RepeatMasker/RepeatMasker -pa 60 -species pig -no_is -dir scaffold_rmask SScrofa.1703USDA.meta.all.qpjp_scf13_edited_at_LARGE1_to_incorporate_scf33.fasta

gunzip GCA_002844635.1_USMARCv1.0_genomic.fna.gz; gunzip GCA_002844635.1_USMARCv1.0_rm.out.gz

# Removing the long chr names and masking lowercase sequence
perl -ne 'chomp; if($_ =~ /^>/){@s = split(/\s+/); print "$s[0]\n";}else{$_ =~ tr/actg/NNNN/; print "$_\n";}' < GCA_002844635.1_USMARCv1.0_genomic.fna > GCA_002844635.1_USMARCv1.0_nmasked.fna

samtools faidx GCA_002844635.1_USMARCv1.0_nmasked.fna
java -jar /mnt/nfs/nfs2/bickhart-users/binaries/GetMaskBedFasta/store/GetMaskBedFasta.jar -f GCA_002844635.1_USMARCv1.0_nmasked.fna -o GCA_002844635.1_USMARCv1.0_nmasked.gaps.bed -s GCA_002844635.1_USMARCv1.0_nmasked.gaps.stast

samtools faidx scaffold_rmask/SScrofa.1703USDA.meta.all.qpjp_scf13_edited_at_LARGE1_to_incorporate_scf33.fasta.masked
cp scaffold_rmask/SScrofa.1703USDA.meta.all.qpjp_scf13_edited_at_LARGE1_to_incorporate_scf33.fasta.masked ./
echo -e "SScrofa.1703USDA.meta.all.qpjp_scf13_edited_at_LARGE1_to_incorporate_scf33.fasta.masked\tCM009090.1" > correction.list

perl reorder_fasta.pl USDA_MARC_Sscrofa_repmasked_chr5_replaced.fasta GCA_002844635.1_USMARCv1.0_nmasked.fna correction.list
```