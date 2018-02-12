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

# The snps aren't present on my HD probe list, so I have to do this the hard way
samtools faidx ARS-UCDv1.0.18.base.fasta X:1-5700000 > ARS-UCD1.0.18.par.region.fasta

sh ../binaries/run_nucmer_plot_automation_script.sh maternal.arrow.fasta ARS-UCD1.0.18.par.region.fasta
```


And generation of a circos diagram with the following data lines:

1. Small box Ideogram of HD chromosome locations
2. Staggered histogram contig locations (+ = top, - = bottom)
3. Links of conflict contigs (> 30% in two locations) with text of contig names

```bash
head -n 29 ../cattle_asms/assembly_revision/ARS-UCD1.0.18.base.fasta.fai | perl -lane 'print "chr - $F[0] $F[0] 0 $F[1] chr$F[0]";' > ARS-UCD1.0.18.circos.karyotype

perl ~/sperl/assembly_scripts/ConvertSNPProbeAlignsToCircos.pl -s maternal.HDprobes.segs -l maternal.arrow.fasta.fai -c maternal.HDprobes.conflicts -k ARS-UCD1.0.18.circos.karyotype -o maternal.HDmap
uniq maternal.HDmap.ctxt > tmp
mv tmp maternal.HDmap.ctxt
perl ../binaries/circos-0.69-6/bin/circos -conf maternal.HDmap.conf

# Now to do this for the paternal haplotype
perl ~/sperl/assembly_scripts/alignAndOrderSnpProbes.pl -a paternal.arrow.fasta -p /mnt/nfs/nfs2/dbickhart/dominette_asm/recombination/BovineHD_B1.probseq.rev1coords.fa -o paternal.HDprobes

perl ~/sperl/assembly_scripts/ConvertSNPProbeAlignsToCircos.pl -s paternal.HDprobes.segs -l paternal.arrow.fasta.fai -c paternal.HDprobes.conflicts -k ARS-UCD1.0.18.circos.karyotype -o paternal.HDmap
uniq paternal.HDmap.ctxt > tmp
mv tmp paternal.HDmap.ctxt

# Had a weird bug where unmapped chromosomes were being placed as incomplete entries
perl ../binaries/circos-0.69-6/bin/circos -conf paternal.HDmap.conf
```

#### Identifying consecutive blocks of missing HD markers

```bash
# Determining marker blocks
perl -e '@blocks = (0); $int = 0; while(<>){chomp; @s = split(/\t/); if(!$int && $s[1] eq "*"){$int = 1; $blocks[-1] += 1;}elsif($int && $s[1] eq "*"){$blocks[-1] += 1;}elsif($int && $s[1] ne "*"){$int = 0; push(@blocks, 0);}} foreach $v (@blocks){print "$v\n";}' < paternal.HDprobes.tab | perl ~/sperl/bed_cnv_fig_table_pipeline/statStd.pl
total   904
Minimum 0
Maximum 344
Average 3.503319
Median  1
Standard Deviation      20.751869
Mode(Highest Distributed Value) 1

# That's actually a big block!
perl -e '@blocks = (0); @chrs = (""); @pos = (""); $int = 0; while(<>){chomp; @s = split(/\t/); if(!$int && $s[1] eq "*"){$int = 1; $blocks[-1] += 1; $chrs[-1] .= $s[4]; $pos[-1] .= $s[5];}elsif($int && $s[1] eq "*"){$blocks[-1] += 1; $pos[-1] .= ";$s[5]";}elsif($int && $s[1] ne "*"){$int = 0; push(@blocks, 0); push(@chrs, ""); push(@pos, "");}} print "chr\tstart\tend\tsize\thdprobenum\n"; for($x = 0; $x < scalar(@blocks); $x++){if($blocks[$x] < 100){next;}@segs = split(";", $pos[$x]); $size = $segs[-1] - $segs[0]; print "$chrs[$x]\t$segs[0]\t$segs[-1]\t$size\t$blocks[$x]\n";}' < paternal.HDprobes.tab > paternal.HDprobes.missing.segs.tab

# Total missing is about 4 megabases

perl -e '@blocks = (0); @chrs = (""); @pos = (""); $int = 0; while(<>){chomp; @s = split(/\t/); if(!$int && $s[1] eq "*"){$int = 1; $blocks[-1] += 1; $chrs[-1] .= $s[4]; $pos[-1] .= $s[5];}elsif($int && $s[1] eq "*"){$blocks[-1] += 1; $pos[-1] .= ";$s[5]";}elsif($int && $s[1] ne "*"){$int = 0; push(@blocks, 0); push(@chrs, ""); push(@pos, "");}} print "chr\tstart\tend\tsize\thdprobenum\n"; for($x = 0; $x < scalar(@blocks); $x++){if($blocks[$x] < 100){next;}@segs = split(";", $pos[$x]); $size = $segs[-1] - $segs[0]; print "$chrs[$x]\t$segs[0]\t$segs[-1]\t$size\t$blocks[$x]\n";}' < maternal.HDprobes.tab > maternal.HDprobes.missing.segs.tab

# again, about 4 megabases

# Combined tab file:
perl -e 'chomp(@ARGV); my %h; foreach my $af (@ARGV){open(IN, "< $af"); while(<IN>){chomp; @s = split(/\t/); push(@{$h{$s[0]}}, $s[1], $s[2]);} close IN;} foreach my $p (sort {$a cmp $b}keys(%h)){@vals = @{$h{$p}}; print "$p\t" . join("\t", @vals); print "\n";}' angusBrahmanF1_FALCONUnzip_arrow.HDprobes.tab paternal.HDprobes.tab maternal.HDprobes.tab > joined_list_HD_marker_coords.tab
```

#### aligning the Brahman x Angus assembly to v23

```bash
for i in *.fasta; do samtools faidx $i; done

# I am going to use a script that automates the comparison of similar chromosomes


```


## Hi-C scaffolding check, assignment and validation

I am going to take the scaffold assignments from the Hi-C data, order them and score them based on expected chromosome coverage, errors and (eventually) use a greedy algorithm to assign them into the best haplotype scaffolds.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/cattle_asms/angus_x_brahman/hic_testing

```bash
# Downloading the necessary files
wget https://gembox.cbcb.umd.edu/triobinning/scaffolding/f1_dam_3ddna.out
wget https://gembox.cbcb.umd.edu/triobinning/scaffolding/f1_dam_phase.out
wget https://gembox.cbcb.umd.edu/triobinning/scaffolding/f1_dam_phasevssalsa.out
wget https://gembox.cbcb.umd.edu/triobinning/scaffolding/f1_dam_salsa.out
wget https://gembox.cbcb.umd.edu/triobinning/scaffolding/f1_sire_3ddna.out
wget https://gembox.cbcb.umd.edu/triobinning/scaffolding/f1_sire_phase.out
wget https://gembox.cbcb.umd.edu/triobinning/scaffolding/f1_sire_phasevssalsa.out
wget https://gembox.cbcb.umd.edu/triobinning/scaffolding/f1_sire_salsa.out
I am redoing the HD probe alignment to generate a stats file that allows the assignment of scaffolds to chromosome units.

sbatch --nodes=1 --mem=15000 --ntasks-per-node=1 --partition=assemble3 --wrap="perl /mnt/nfs/nfs2/dbickhart/dominette_asm/recombination/alignAndOrderSnpProbes.pl -a bostaurus_angus.fasta -p /mnt/nfs/nfs2/dbickhart/dominette_asm/recombination/BovineHD_B1.probseq.rev1coords.fa -o bostaurus_angus.HDProbes"

sbatch --nodes=1 --mem=15000 --ntasks-per-node=1 --partition=assemble3 --wrap="perl /mnt/nfs/nfs2/dbickhart/dominette_asm/recombination/alignAndOrderSnpProbes.pl -a bostaurus_brahma.fasta -p /mnt/nfs/nfs2/dbickhart/dominette_asm/recombination/BovineHD_B1.probseq.rev1coords.fa -o bostaurus_brahma.HDProbes"

# First, I made modifications to my previous segmentation script to use the mashmap alignments that Serge generated
# Removing the alignments to the unplaced contigs
grep -v 'Scb' f1_dam_3ddna.out > f1_dam_3ddna.out.filt
perl ~/sperl/assembly_scripts/alignAndOrderSnpProbes.pl -g f1_dam_3ddna.out.filt -o f1_dam_3ddna.filt

# Going to try to automate the process
for i in *.out; do grep -v 'Scb' $i > $i.filt; perl ~/sperl/assembly_scripts/alignAndOrderSnpProbes.pl -g $i.filt -o $i.filt; done


```

## Repeatmasker analysis and putative active retroviral element

Serge found a poly-A region in the assembly that NCBI was trying to tell us was a Solid adaptor sequence. Here is the region in the Angus (sire) assembly:

> tig00020254    84863977    15197319..15197349    adaptor:NGB00970.1


I need to repeatmask and pull relevant repetitive data from the region directly flanking that candidate. The Poly-A tail should be on the 3' end of a relatively "intact" ERV for this to be convincing.

> /mnt/nfs/nfs2/bickhart-users/cattle_asms/angus_x_brahman

```bash
# Generating repeatmasker files
# First, let's get rid of the annoying arrow tags
perl -ne '$_ =~ s/\|arrow\|arrow//; print $_;' < bostaurus_angus.fasta > bostaurus_angus.reformat.fasta

mkdir angus_rm
sbatch --nodes=1 --ntasks-per-node=60 --mem=100000 -p assemble1 --wrap="/mnt/nfs/nfs2/bickhart-users/binaries/RepeatMasker/RepeatMasker -pa 60 -species cow -no_is -dir angus_rm bostaurus_angus.reformat.fasta"

