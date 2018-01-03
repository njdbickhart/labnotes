# IGC kmer analysis in ruminants
---
*1/2/2018*

These are my notes on a kmer-based analysis of missing regions in ruminant genomes -- most immediately in goats.

## Table of Contents



## Goat Analysis

I selected a group of 100 goat WGS datasets from SRA using the default "goat" and "DNA" filters. The SRA entries were as follows:


In order to limit the amount of space the analysis files will take on the server, I am going to pipe the SRA content into jellyfish/Mash.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/goat_projects

```bash
# Generating the mash profiles
perl -ne '$_ =~ s/\"//g; @s = split(/,/, $_); print join("\t", @s);' < goat_wgs_sra.csv > goat_wgs_sra.tab

perl -lane '@s = split(/,/, $F[0]); if($s[0] =~ /^Run/){next;} print $s[0]' < goat_wgs_sra_runifo.csv | xargs -I {} sbatch --nodes=1 --ntasks-per-node=2 --mem=5000 --wrap="fastq-dump -Z {} | /mnt/nfs/nfs2/bickhart-users/binaries/mash-Linux64-v1.1.1/mash sketch -b 1G -s 25000 -o {} - "

# Generating the jf databases
perl -lane '@s = split(/,/, $F[0]); if($s[0] =~ /^Run/){next;} print $s[0]' < goat_wgs_sra_runifo.csv | xargs -I {} sbatch --nodes=1 --ntasks-per-node=2 --mem=110000 -p assemble3 --wrap="module load sratoolkit/2.8.1-3; module load jellyfish/2.2.3; jellyfish count -s 3G -m 21 --bf-size 100G -t 10 -o {}.jf <(fastq-dump -Z {} )";

perl -lane '@s = split(/,/, $F[0]); if($s[0] =~ /^Run/){next;} print $s[0]' < goat_wgs_sra_runifo.csv | xargs -I {} sbatch --nodes=1 --ntasks-per-node=2 --mem=160000 -p assemble3 --wrap="module load sratoolkit/2.8.1-3; module load jellyfish/2.2.3; echo 'fastq-dump -Z {}' > {}.g; jellyfish count -s 3G -m 21 --bf-size 100G -t 10 -o {}.jf -g {}.g";

# Mashing some cattle profiles as an outgroup
sbatch dominette_hiseq_mash.sh
```

And here is the code in that simple script:

```bash
#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=15000
#SBATCH --ntasks-per-node=10
#SBATCH -p assemble3

for i in /mnt/nfs/nfs2/brosen/projects/ARS-UCD/Dominette_HiSeq_data/*.fastq.gz
do
        gunzip -c $i
done | /mnt/nfs/nfs2/bickhart-users/binaries/mash-Linux64-v1.1.1/mash sketch -b 1G -s 25000 -p 9 -o dominette_hiseq_mash -
```

Now to generate the distance matrix for mash and plot it as a PCA plot.

```bash
ls *.msh > combined_mash_sketches.list

/mnt/nfs/nfs2/bickhart-users/binaries/mash-Linux64-v1.1.1/mash paste -l combined_mash_sketches combined_mash_sketches.list
for i in `cat combined_mash_sketches.list`; do echo $i; /mnt/nfs/nfs2/bickhart-users/binaries/mash-Linux64-v1.1.1/mash dist -p 10 -t combined_mash_sketches.msh $i > ${i}.dist; done

# Now to generate a group assignment file
perl -ne '@b = split(/,/); print "$b[0]\t$b[10]\t$b[14]\t$b[19]\t$b[20]\t$b[24]\t$b[26]\t$b[29]\t$b[31]\t$b[33]\t$b[34]\t$b[40]\t$b[41]\n";' < goat_wgs_sra_runifo.csv > goat_wgs_sra_runifo.abbrev.tab

perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); @names; while(<IN>){chomp; @d = split(/\./); push(@names, $d[0]);} close IN; open(IN, "< $ARGV[1]"); <IN>; %h; while(<IN>){chomp; @s = split(/\t/); $h{$s[0]} = $s[7];} close IN; foreach my $b (@names){$d = ($h{$b} =~ /SAN/)? "SANEN" : "CASHMERE"; print "$b\t$h{$b}\t$d\n";}' combined_mash_sketches.list goat_wgs_sra_runifo.abbrev.tab > combined_mash_groups.tab

# I had to edit a few names to make it work. Now to make the distance matrix
perl -lane 'open(IN, "$F[0].dist"); $h = <IN>; $d = <IN>; chomp $d; @dsegs = split(/\t/, $d); @nm = split(/\./, $F[0]); $dsegs[0] = $nm[0]; print join("\t", @dsegs);' < combined_mash_sketches.list > combined_mash_distances.matrix
```

Loading things into R so that I can generate the plots and perform the analysis.

```R
library(vegan)
library(dplyr)
mash_data <- read.delim("combined_mash_distances.matrix", header=FALSE)
entries <- mash_data[,1]
mash_data.format <- select(mash_data, -V1)

rownames(mash_data.format) <- entries
colnames(mash_data.format) <- entries

groups <- read.delim("combined_mash_groups.tab", header=FALSE)
colnames(groups) <- c("run", "sample", "group")
groups$pch[groups$group == "COW"] = 17
groups$pch[groups$group == "CASHMERE"] = 0
groups$pch[groups$group == "SANEN"] = 16
groups$pch[groups$group == "BOER"] = 3

groups$col[groups$group == "COW"] <- "black"
groups$col[groups$group == "CASHMERE"] <- "dodgerblue"
groups$col[groups$group == "SANEN"] <- "darkorange2"
groups$col[groups$group == "BOER"] <- "chartreuse3"

pcoa.result <- wcmdscale(mash_data.format, eig=TRUE, add="lingoes", k=75)

pdf(file="goat_pcoa_plot.pdf", useDingbats=FALSE)
plot(pcoa.result$points, type = "p", col=groups$col, pch = groups$pch, asp =1/1)
ordihull(pcoa.result$points, groups = groups$group, label = FALSE, col = c("black", "darkorange2", "dodgerblue", "chartreuse3"))
dev.off()
```

The results of the ordination showed a clear division between the different breeds, and some variability within them. I need to spruce up the quality of the plot, but it's a good start. Now to 

### KAT analysis of unique kmers

Rather than roll my own, I'm actually happy to see that someone has generated a tool for comparison of kmers in assemblies.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/goat_projects

```bash
# Just a trial run here...
/mnt/nfs/nfs2/bickhart-users/binaries/KAT/kat-2.3.4/src/kat comp -t 15 -m 21 -n -p 'png' SRR3144624.jf SRR5557710.jf

# I want to see if any of the sequence datasets have unique kmers not present in the assembly
# first, unzip the assembly to a file
gunzip -c /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v13/papadum-v13.full.fa.gz > ./ARS1.goat.reference.fasta
# generate the jellyfish database for goat
jellyfish count -s 3G -m 21 -t 20 -o ARS1.goat.reference.jf ./ARS1.goat.reference.fasta

# generating kmer count plots for each dataset
for i in *.jf; do name=`echo $i | cut -d'.' -f1`; echo $name; sbatch --nodes=1 --mem=45000 --ntasks-per-node=15 --wrap="/mnt/nfs/nfs2/bickhart-users/binaries/KAT/kat-2.3.4/src/kat comp -t 15 -o $name.kat $i ARS1.goat.reference.jf"; done
```
