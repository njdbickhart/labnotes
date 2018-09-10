# Creating a custom SNPEff database from a gff3 file
---
*9/10/2018*

These are my notes on how to generate the binned SNPeff db for the ARS-UCDv1.2 reference.

## Table of Contents
* [Generating the custom db information](#generation)
* [How to install the custom database](#installation)

<a name="generation"></a>
## Generating the custom db information

These are my notes on how to 

> Assembler2: /mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi

```bash
for i in `seq 1 29` X MT; do name=bt_ref_ARS-UCD1.2_chr${i}.fa; echo $name; perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); open(OUT, ">> $ARGV[1]"); while(<IN>){chomp; @F = split(/\s+/); if($F[0] =~ /^>/){@bsegs = split(/\s+/, $F[0]); $bsegs[0] =~ s/ref\|//; $bsegs[0] =~ s/\|//g; print ">$ARGV[2]\n"; $bsegs[0] =~ s/>//; print OUT "$ARGV[2]\t$bsegs[0]\n";}else{print "$F[0]\n";}}' $name ARS-UCD_ncbi_to_chrnum.tab $i >> ARSUCD1.2.current_ref.fa; done

perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %h; while(<IN>){chomp; @s = split(/\t/); $h{$s[1]} = $s[0];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); if(exists($h{$s[0]})){$t = $h{$s[0]}; shift(@s); unshift(@s, $t); print join("\t", @s) . "\n";}elsif($s[0] =~ /^#/){print join("\t", @s) . "\n";}} close IN;' ARS-UCD_ncbi_to_chrnum.tab ../ars_ucd_125/ref_ARS-UCD1.2_top_level.gff3 > ../ars_ucd_125/ARSUCD1.2_annotation_regnum.gff3
```

> Assembler2: /mnt/nfs/nfs2/bickhart-users/binaries/snpEff/data/ARSUCD1.2/

```bash
cp /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_125/ARSUCD1.2_annotation_regnum.gff3 ./genes.gff
echo "##FASTA" >> ./genes.gff
cat /mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARSUCD1.2.current_ref.fa >> ./genes.gff
```

> Assembler2: /mnt/nfs/nfs2/bickhart-users/binaries/snpEff/

```bash
echo "ARSUCD1.2" >> snpEff.config
java -jar snpEff.jar build -gff3 ARSUCD1.2
tar -czvf ARSUCD1.2.tar.gz ARSUCD1.2
```

<a name="installation"></a>
## How to install the custom database

First, enter into your SNPeff installation directory (the base directory. You should see the "SnpEFF.jar" file if you do an "ls" in the current directory). If it doesn't already contain a "data" folder, then make one:

```bash
mkdir data
```

Next, you will copy the ARSUCD1.2.tar.gz tarball into the data directory and unpack.

```bash
cp ARSUCD1.2.tar.gz data/
cd data
tar -xvf ARSUCD1.2.tar.gz
```

If you are using the current version of SnpEff, all you will need to do to install the database will be to modify your SNPEff config file as follows:

```bash
# In the SNPeff installation directory:
echo "ARSUCD1.2" >> snpEff.config

# It won't hurt to try to build the database again:
java -jar snpEff.jar build -gff3 ARSUCD1.2
```

In order to annotate VCF files with the new database, run the following command:

```bash
java -jar (path)/(to)/snpEff.jar ARSUCD1.2 (yourvcf).vcf > (youroutput).vcf
```