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

# I am using bedops to preserve as much information as possible from the masking
/mnt/nfs/nfs2/bickhart-users/binaries/bin/rmsk2bed < angus_rm/bostaurus_angus.reformat.fasta.out > angus_rm/bostaurus_angus.reformat.fasta.rm.bed

# Bedops has an "optional" 16th column that gives bedtools a fit
perl -lane 'print "$F[0]\t$F[1]\t$F[2]\t$F[3]";' < angus_rm/bostaurus_angus.reformat.fasta.rm.bed > angus_rm/bostaurus_angus.reformat.fasta.rm.simp.bed
echo -e "tig00020254\t15197319\t15197349" | bedtools intersect -a stdin -b angus_rm/bostaurus_angus.reformat.fasta.rm.simp.bed -wa -wb
tig00020254     15197319        15197349        tig00020254     15197050        15197340        L1_BT

# It's a line1 nearby, let's expand the search a bit to see if there are complete elements nearby
echo -e "tig00020254\t15190319\t15207349" | bedtools intersect -a stdin -b angus_rm/bostaurus_angus.reformat.fasta.rm.simp.bed -wa -wb
tig00020254     15190319        15207349        tig00020254     15193973        15194183        L1-2_BT
tig00020254     15190319        15207349        tig00020254     15194206        15194346        MIRb
tig00020254     15190319        15207349        tig00020254     15194577        15194795        Charlie10a
tig00020254     15190319        15207349        tig00020254     15195078        15195267        Bov-tA3
tig00020254     15190319        15207349        tig00020254     15195280        15195394        L1MC3
tig00020254     15190319        15207349        tig00020254     15195464        15195704        MIR
tig00020254     15190319        15207349        tig00020254     15197050        15197340        L1_BT  <- direct intersection
tig00020254     15190319        15207349        tig00020254     15197937        15198085        L1MA9
tig00020254     15190319        15207349        tig00020254     15198091        15198214        Bov-tA2
tig00020254     15190319        15207349        tig00020254     15198319        15198568        CHR-2B
tig00020254     15190319        15207349        tig00020254     15198770        15198887        L2a
```

Bed ops column headers can be found at [this webpage](http://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/rmsk2bed.html). I think that the identification of poly-As or poly-As with some hamming distance cutoff may be a good criterion to identify recent transposition events in the genome. 

## Creating full assemblies based on scaffold winners

This will get a bit more complex, however, I need to generate separate full references for each haplotype-resolved assembly. Most of the "winner" chromosomes are easy as they are single scaffolds.

Let's start with those and then work on the remaining chromosomes by hand.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/cattle_asms/angus_x_brahman/hic_testing

```bash
# How to get single scaffold winners
perl -e '%row = ("salsa" => 3, "phase" => 8, "3ddna2" => 13, "3ddna" => 18); <>; while(<>){chomp; @s = split(/\t/); if($s[$row{$s[1]}] == 1){$scaff = $s[$row{$s[1]} - 1]; print "$s[0]\t$s[1]\t$scaff\n";}}' < dam_scaffs_consensus.retry2.tab
perl -e '%row = ("salsa" => 3, "phase" => 8, "3ddna2" => 13, "3ddna" => 18); <>; while(<>){chomp; @s = split(/\t/); if($s[$row{$s[1]}] == 1){$scaff = $s[$row{$s[1]} - 1]; print "$s[0]\t$s[1]\t$scaff\n";}}' < sire_scaffs_consensus.retry2.tab

# Let's work on the above formula to try to get the segments going that we need for an agp format
# Well, I don't have all of the information to generate a true AGP, so lets go with a transitional format first
# The last columns contain the technology name of the scaffolds
perl -e '%row = ("salsa" => 3, "phase" => 8, "3ddna2" => 13, "3ddna" => 18); <>; while(<>){chomp; @s = split(/\t/); $scaff = $s[$row{$s[1]} - 1]; @ss = split(/;/, $scaff); $p = 0; for(my $x = 0; $x < scalar(@ss); $x++){$j = $ss[$x]; $p++; $sname = $j =~ s/,.+$//; print "$s[0]\t0\t0\t$p\tA\t$j\t0\t0\t?\t$s[1]\n"; unless($x + 1 >= scalar(@ss)){$p++; print "$s[0]\t0\t0\t$p\tU\t100\tscaffold\tyes\tmap\n";}}}' < dam_scaffs_consensus.retry2.tab > dam_scaffs_consensus.retry2.pagp

perl -e '%row = ("salsa" => 3, "phase" => 8, "3ddna2" => 13, "3ddna" => 18); <>; while(<>){chomp; @s = split(/\t/); $scaff = $s[$row{$s[1]} - 1]; @ss = split(/;/, $scaff); $p = 0; for(my $x = 0; $x < scalar(@ss); $x++){$j = $ss[$x]; $p++; $sname = $j =~ s/,.+$//; print "$s[0]\t0\t0\t$p\tA\t$j\t0\t0\t?\t$s[1]\n"; unless($x + 1 >= scalar(@ss)){$p++; print "$s[0]\t0\t0\t$p\tU\t100\tscaffold\tyes\tmap\n";}}}' < sire_scaffs_consensus.retry2.tab > sire_scaffs_consensus.retry2.pagp

# Now to download the individual fastas and then create a script to process them.
# I plan to fasta index them, then grep out the individual components for subsequent entry into SIFF
# First, I need to determine actual alignments and breakpoints of scaffolds in the winners categories of my tables and use that info to update the agp files.
module load samtools; for i in *.fasta; do echo $i; samtools faidx $i; done

# Testing assumption that scaffolds are present only once in the pseudoagp files I created
perl ~/sperl/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f dam_scaffs_consensus.retry2.pagp -c 5

# Notes: the salsa dam scaffold (scaffold_608) has a portion of the telomere of chr1 at 119355000 - 119500365, and should be broken on chr2 at 1-119354684 and 119517501
# Notes: the 3ddna_v2 sire scaffold (HiC_scaffold_18) has a portion of the near end of chr15 in the beginning of the scaffold. 

perl prepare_master_fasta.pl sire_scaffs_consensus.retry2.pagp sire_scaffold_fastas.tab sire_scaffold_remake

# the script worked, but it needs allot of manual edits to make the final fasta. At least I have the components ready!
# Note, I noticed for the Sire that Phase  had the shorted condensed segments, around the same regions as the PAR in the X. Since this is the sire haplotype, the only remaining content would be the Y or PAR. I relabelled and used Phase instead of 3dDNA here in the agp

# Now to do the same with the dam agp file
perl prepare_master_fasta.pl dam_scaffs_consensus.retry2.pagp dam_scaffold_fastas.tab dam_scaffold_remake
```

#### Changes list for sire
	* PAR changes listed above
	* Removed 3ddna2 scaffold 3 from chr1 (small alignment
	* Note: chr2 scaffold phase is circularized
	* Note: chr3 salsa scaffold 475 needs to be in the middle of scaffold 824
	* Note: chr7 phase scaffold is all sorts of messed up
	* Note: same phase issue with chr10
	* note: Phase chr14 has same issues
	* note: same with phase chr18
	* Removed Salsa chr21 scaffold_1306 (small alignment)
	

#### Changes list for dam
	* Removed salsa scaffolds 556 and 608 from chr1 (small alignment)
	* Note: a portion of salsa scaffold_608 contains the chr1 telomere -- it was not modified or removed in chr2
	* Note: phase chr5 scaffold has lots of orientation errors
	* Note: salsa chr10 scaffold is circularized
	* Note: salsa chr12 scaffold_671 is circularized
	* Note: phase chr22 is circularized
	* Note: phase chr23 scaffold_513 is in the middle of scaffold_921 and at the beginning of the chr
	* Note: 3ddna chr28 asm_hic_scaffold_21 is circularized

```bash
# Adding final components to the fastas, fixing the agps and generating the renamed, concatenated fasta
# For Sire
perl -ne '$_ =~ s/\?/+/g; print $_;' < sire_scaffold_remake.agp > temp.agp
mv temp.agp sire_scaffold_remake.agp

# adding the PAR/Y scaffolds
samtools faidx f1_sire_phase.fasta PGA_scaffold28__120_contigs__length_19643512 PGA_scaffold27__152_contigs__length_25761177 > sire_scaffold_remake.adds.fasta
cat sire_scaffold_remake.fasta sire_scaffold_remake.adds.fasta > sire_scaffold_remake.full.fasta
samtools faidx /mnt/nfs/nfs2/bickhart-users/cattle_asms/angus_x_brahman/hic_testing/sire_scaffold_remake.full.fasta

# OK, here goes nothing!
java -Xmx100g -jar CombineFasta.jar agp2fasta -f sire_scaffold_remake.full.fasta -a sire_scaffold_remake.agp -o sire_best_scaffold_reference.fa
# It worked!

# For Dam
perl -ne '$_ =~ s/\?/+/g; print $_;' < dam_scaffold_remake.agp > temp.agp
mv temp.agp dam_scaffold_remake.agp
samtools faidx dam_scaffold_remake.fasta

# OK, this is the dam fasta generation
java -Xmx100g -jar CombineFasta.jar agp2fasta -f dam_scaffold_remake.fasta -a dam_scaffold_remake.agp -o dam_best_scaffold_reference.fa
```

#### Tracking contig alignments and finding missing sequence

I am going to align the original sire and dam scaffolds to the "best scaffold" fasta and track their locations. I will use minimap to try to keep track of everything.

| Haplotype assembly| contig length | scaffold length | difference|
|:---| ---:|---:|---:|
|Sire | 2573806573 | 2541684135 | 32,122,438 |
|Dam | 2678769087 | 2678438895 | 330,192 |


```bash
for i in f1_dam_3ddna.fasta f1_dam_3ddna_v2.fasta f1_dam_phase.fasta f1_dam_salsa.fasta; do echo $i; sbatch --nodes=1 --ntasks-per-node=1 --mem=25000 -p assemble1 --wrap="/mnt/nfs/nfs2/bickhart-users/binaries/minimap2/minimap2 -x asm5 ../bostaurus_brahma.fasta $i > $i.brahma.paf"; done

for i in f1_sire_3ddna.fasta f1_sire_3ddna_v2.fasta f1_sire_phase.fasta f1_sire_salsa.fasta; do echo $i; sbatch --nodes=1 --ntasks-per-node=1 --mem=25000 -p assemble1 --wrap="/mnt/nfs/nfs2/bickhart-users/binaries/minimap2/minimap2 -x asm5 ../bostaurus_angus.fasta $i > $i.angus.paf"; done

# Now I need to develop a script that converts the scaffolds into bed coordinates for comparisons
# First, let's get rid of the annoying, annoying pipes in the damn contig names!
for i in *.paf; do perl -ne '$_=~ s/\|arrow\|arrow//; print $_;' < $i > temp; mv temp $i; done

# Due to alignment ambiguity, I want to "squash" the first few bases and last few bases in the paf alignments to make sure that we're not missing thousands of bases from the starts and ends of the contigs in my comparison
for i in *.paf; do perl -lane 'if($F[7] < 1000){$F[7] = 1;} if(abs($F[8] - $F[6]) < 1000){ $F[8] = $F[6];} print join("\t", @F);' < $i > $i.tpaf; done

# Now to create a translation script
ls f1_dam*.tpaf > dam_scaffolding_tech_pafs.tab
perl convert_agp_to_bed.pl dam_scaffolding_tech_pafs.tab dam_scaffold_remake.agp dam_scaffold_remake.contigmap.bed
perl -lane '$F[0] =~ s/\|arrow\|arrow//; print "$F[0]\t1\t$F[1]";' < ../bostaurus_brahma.fasta.fai > dam_assembly_contig_lengths.bed

# OK, so now it's just a bedtools exercise
bedtools subtract -a dam_assembly_contig_lengths.bed -b dam_scaffold_remake.contigmap.bed | perl -lane '$l = $F[2] - $F[1]; print "$_\t$l";' > dam_assembly_contig_len_missing.bed
cat dam_assembly_contig_len_missing.bed | perl ~/sperl/bed_cnv_fig_table_pipeline/bed_length_sum.pl
        Interval Numbers:       1764
        Total Length:           67078728	<- OK, this is double the amount predicted to be in the other animal!
        Length Average:         38026.4897959184
        Length Median:          26274.5
        Length Stdev:           161369.420124495
        Smallest Length:        73
        Largest Length:         5525021

bedtools merge -i dam_scaffold_remake.contigmap.bed -c 4 -o distinct > dam_scaffold_remake.contigmap.merged.bed

bedtools merge -i dam_scaffold_remake.contigmap.bed -c 4 -o distinct | grep ',' | perl ~/sperl/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 3
Entry   Count
3ddna,3ddna2    17
3ddna,3ddna2,salsa      3
3ddna,phase     2
3ddna,salsa     22
3ddna2,phase    4
3ddna2,salsa    5

bedtools merge -i dam_scaffold_remake.contigmap.bed -c 4 -o distinct | grep ',' | perl -lane 'system("grep $F[0] dam_scaffold_remake.contigmap.bed");' > dam_scaffold_remake.duptable.bed

perl -e '%c; while(<>){chomp; @s = split(/\t/); $c{$s[0]}->{$s[3]}->{$s[5]} += $s[2] - $s[1];} print "Contig\tTech\tChr\tLength\n"; foreach my $t (keys(%c)){foreach my $s (keys(%{$c{$t}})){foreach my $c (keys(%{$c{$t}->{$s}})){print "$t\t$s\t$c\t" . $c{$t}->{$s}->{$c} . "\n";}}}' < dam_scaffold_remake.duptable.bed > dam_scaffold_remake.duptable.tab

perl -e '<>; $sum = 0; while($f =<>){$s = <>; chomp $f, $s; @j = split(/\t/, $f); @h = split(/\t/, $s); $min = ($j[-1] < $h[-1])? $j[-1] : $h[-1]; $sum += $min;} print "$sum\n";' < dam_scaffold_remake.duptable.tab
82,231,839

# The main issue was 3ddna. Let's see how the scaffolds perform in terms of duplication when one of the 3ddna scaffolds is removed.

```

#### Checking duplication after selective scaffold reduction

I am going to remove a scaffolding technology to see how much duplication that removes.

```bash
# The pipeline for the dam
cp dam_scaffs.list dam_scaffs.no3d.list
perl compare_algns_per_chrassignment.pl -l dam_scaffs.no3d.list -f ../../ars_ucd_123/ARS-UCDv1.0.23.fasta.fai -t 50000 -o dam_scaffs_consensus.no3ddna.tab

perl -e '%row = ("salsa" => 3, "phase" => 13, "3ddna2" => 8); <>; while(<>){chomp; @s = split(/\t/); $scaff = $s[$row{$s[1]} - 1]; @ss = split(/;/, $scaff); $p = 0; for(my $x = 0; $x < scalar(@ss); $x++){$j = $ss[$x]; $p++; $sname = $j =~ s/,.+$//; print "$s[0]\t0\t0\t$p\tA\t$j\t0\t0\t?\t$s[1]\n"; unless($x + 1 >= scalar(@ss)){$p++; print "$s[0]\t0\t0\t$p\tU\t100\tscaffold\tyes\tmap\n";}}}' < dam_scaffs_consensus.no3ddna.tab > dam_scaffs_consensus.no3ddna.pagp

perl prepare_master_fasta.pl dam_scaffs_consensus.no3ddna.pagp dam_scaffold_fastas.tab dam_scaffold_no3ddna
# I still need to check scaffold orientation from the alignments, but that can be changed later on. Let's check the duplication/deletion counts
# I needed to delete scaffold_608 from the end of chr1 from the salsa scaffolder. Also 556. Basically, it looks like I will be following the same changes list as above!

perl convert_agp_to_bed.pl dam_scaffolding_tech_pafs.tab dam_scaffold_no3ddna.agp dam_scaffold_no3ddna.contigmap.bed

bedtools subtract -a dam_assembly_contig_lengths.bed -b dam_scaffold_no3ddna.contigmap.bed | perl ~/sperl/bed_cnv_fig_table_pipeline/bed_length_sum.pl
        Interval Numbers:       1409
        Total Length:           54437778
        Length Average:         38635.75443577
        Length Median:          27992
        Length Stdev:           125096.609047426
        Smallest Length:        73
        Largest Length:         3582770

# That's a little better. Let's see where the duplicates exist again
bedtools merge -i dam_scaffold_no3ddna.contigmap.bed -c 4 -o distinct | grep ',' | perl -lane 'system("grep $F[0] dam_scaffold_no3ddna.contigmap.bed");' > dam_scaffold_no3ddna.duptable.bed
perl -e '%c; while(<>){chomp; @s = split(/\t/); $c{$s[0]}->{$s[3]}->{$s[5]} += $s[2] - $s[1];} print "Contig\tTech\tChr\tLength\n"; foreach my $t (keys(%c)){foreach my $s (keys(%{$c{$t}})){foreach my $c (keys(%{$c{$t}->{$s}})){print "$t\t$s\t$c\t" . $c{$t}->{$s}->{$c} . "\n";}}}' < dam_scaffold_no3ddna.duptable.bed > dam_scaffold_no3ddna.duptable.tab

# Obviously, some of these are really small alignments that only need to be counted once
perl -e '<>; $sum = 0; while($f =<>){$s = <>; chomp $f, $s; @j = split(/\t/, $f); @h = split(/\t/, $s); $min = ($j[-1] < $h[-1])? $j[-1] : $h[-1]; $sum += $min;} print "$sum\n";' < dam_scaffold_no3ddna.duptable.tab
9,903,170

# Checking the number of contigs missing:
bedtools intersect -a dam_assembly_contig_lengths.bed -b dam_scaffold_no3ddna.contigmap.bed -v | wc -l
966
bedtools intersect -a dam_assembly_contig_lengths.bed -b dam_scaffold_no3ddna.contigmap.bed -v | perl ~/sperl/bed_cnv_fig_table_pipeline/bed_length_sum.pl
        Interval Numbers:       966
        Total Length:           43,025,319
        Length Average:         44539.6677018634
        Length Median:          33831.5
        Length Stdev:           96728.8496640135
        Smallest Length:        1039
        Largest Length:         2591174

# Yup, it's the smaller ones!
bedtools intersect -a dam_assembly_contig_lengths.bed -b dam_scaffold_no3ddna.contigmap.bed -v | perl -lane 'print "$F[2]";' > dam_no3ddna_totallymissing_contiglens.list
perl -lane 'print "$F[1]";' < ../bostaurus_brahma.fasta.fai > dam_total_contiglens.list

############################
#### Quick mashmap test ####
############################
# I just wanted to see how mashmap stacks up against my minimap2 alignments
# I think that minimap2 may have aligned incredibly short fragments of the contigs over repetitive regions
# Mashmap has a "one-to-one" mapping that would eliminate these spurious alignments
for i in f1_dam_3ddna_v2.fasta f1_dam_phase.fasta f1_dam_salsa.fasta; do echo $i; echo ../bostaurus_brahma.fasta; sbatch --nodes=1 --ntasks-per-node=5 --mem=25000 -p assemble1 --wrap="../../../binaries/mashmap-Linux64-v2.0/mashmap -r ../bostaurus_brahma.fasta -q $i -t 5 -f one-to-one -o $i.mashmap"; done

for i in *.mashmap; do echo $i; perl -ne '$_ =~ s/\|arrow\|arrow//; print $_;' < $i > temp; mv temp $i; done
ls *.mashmap > dam_scaffolding_tech_mashmap.tab
perl convert_agp_to_bed.pl dam_scaffolding_tech_mashmap.tab dam_scaffold_no3ddna.agp dam_scaffold_no3d_mash.contigmap.bed

bedtools subtract -a dam_assembly_contig_lengths.bed -b dam_scaffold_no3d_mash.contigmap.bed| perl ~/sperl/bed_cnv_fig_table_pipeline/bed_length_sum.pl
        Interval Numbers:       1942
        Total Length:           55459811	<- about 1 mbp larger than lasttime
        Length Average:         28558.0901132853
        Length Median:          19922.5
        Length Stdev:           108011.018957934
        Smallest Length:        1
        Largest Length:         3582717

bedtools merge -i dam_scaffold_no3d_mash.contigmap.bed -c 4 -o distinct | grep ',' | perl -lane 'system("grep $F[0] dam_scaffold_no3ddna.contigmap.bed");' > dam_scaffold_no3d_mash.duptable.bed
perl -e '%c; while(<>){chomp; @s = split(/\t/); $c{$s[0]}->{$s[3]}->{$s[5]} += $s[2] - $s[1];} print "Contig\tTech\tChr\tLength\n"; foreach my $t (keys(%c)){foreach my $s (keys(%{$c{$t}})){foreach my $c (keys(%{$c{$t}->{$s}})){print "$t\t$s\t$c\t" . $c{$t}->{$s}->{$c} . "\n";}}}' < dam_scaffold_no3d_mash.duptable.bed > dam_scaffold_no3d_mash.duptable.tab

```

In summary, the winner scaffold strategy runs into problems when I use an error-correcting scaffolder with another non-error corrector. The number of putative errors is similar to what we found in Goat (30 - 100 ish), but defies easy correction. I'm going to try a quick comparison with just phase and salsa (no 3ddna at all). 

```bash
cp dam_scaffs.no3d.list dam_scaffs.sans3d.list
perl compare_algns_per_chrassignment.pl -l dam_scaffs.sans3d.list -f ../../ars_ucd_123/ARS-UCDv1.0.23.fasta.fai -t 50000 -o dam_scaffs_consensus.sans3ddna.tab

perl -e '%row = ("salsa" => 3, "phase" => 8); <>; while(<>){chomp; @s = split(/\t/); $scaff = $s[$row{$s[1]} - 1]; @ss = split(/;/, $scaff); $p = 0; for(my $x = 0; $x < scalar(@ss); $x++){$j = $ss[$x]; $p++; $sname = $j =~ s/,.+$//; print "$s[0]\t0\t0\t$p\tA\t$j\t0\t0\t?\t$s[1]\n"; unless($x + 1 >= scalar(@ss)){$p++; print "$s[0]\t0\t0\t$p\tU\t100\tscaffold\tyes\tmap\n";}}}' < dam_scaffs_consensus.sans3ddna.tab > dam_scaffs_consensus.sans3ddna.pagp

perl prepare_master_fasta.pl dam_scaffs_consensus.sans3ddna.pagp dam_scaffold_fastas.tab dam_scaffold_sans3ddna

# I made the chr1 edits listed above, but it looks like chr7 and 15 have problems too
# Chr7 removed salsa scaffold_815 from chr7  (short alignment)
# chr15 the first 50 megs of salsa scaffold_539 map to chr13  <- corrected in the agp
# chr15 the last 230 kb of salsa scaffold_513 maps to chr23		<- removed in the agp

perl convert_agp_to_bed.pl dam_scaffolding_tech_mashmap.tab dam_scaffold_sans3ddna.agp dam_scaffold_sans3ddna.contigmap.bed

bedtools subtract -a dam_assembly_contig_lengths.bed -b dam_scaffold_sans3ddna.contigmap.bed | perl ~/sperl/bed_cnv_fig_table_pipeline/bed_length_sum.pl
        Interval Numbers:       1748
        Total Length:           47861769		<- the best one yet!
        Length Average:         27380.8747139588
        Length Median:          17700
        Length Stdev:           80070.0348587128
        Smallest Length:        1
        Largest Length:         2591174

bedtools merge -i dam_scaffold_sans3ddna.contigmap.bed -c 4 -o distinct | grep ',' | perl -lane 'system("grep $F[0] dam_scaffold_no3ddna.contigmap.bed");' > dam_scaffold_sans3d_mash.duptable.bed
```

## Generating salsa scaffolds and assessing scaffold contiguity

I am going to hijack my pipeline to create chromosome scaffolds from the Salsa scaffold list with minimal modifications.

```bash
# Starting with the dam
cp dam_scaffs.sans3d.list dam_scaffs.salsa.list

perl compare_algns_per_chrassignment.pl -l dam_scaffs.salsa.list -f ../../ars_ucd_123/ARS-UCDv1.0.23.fasta.fai -t 50000 -o dam_scaffs_consensus.salsa.tab
# 11 chromosomes were a single salsa scaffold

perl -e '%row = ("salsa" => 3, "phase" => 8); <>; while(<>){chomp; @s = split(/\t/); $scaff = $s[$row{$s[1]} - 1]; @ss = split(/;/, $scaff); $p = 0; for(my $x = 0; $x < scalar(@ss); $x++){$j = $ss[$x]; $p++; $sname = $j =~ s/,.+$//; print "$s[0]\t0\t0\t$p\tA\t$j\t0\t0\t?\t$s[1]\n"; unless($x + 1 >= scalar(@ss)){$p++; print "$s[0]\t0\t0\t$p\tU\t100\tscaffold\tyes\tmap\n";}}}' < dam_scaffs_consensus.salsa.tab > dam_scaffs_consensus.salsa.pagp

perl prepare_master_fasta.pl dam_scaffs_consensus.salsa.pagp dam_scaffold_fastas.tab dam_scaffold_salsa

perl -ne '$_ =~ s/\?/+/g; print $_;' < dam_scaffold_salsa.agp > temp.agp
mv temp.agp dam_scaffold_salsa.agp

samtools faidx dam_scaffold_salsa.fasta
java -Xmx100g -jar CombineFasta.jar agp2fasta -f dam_scaffold_salsa.fasta -a dam_scaffold_salsa.agp -o dam_scaffold_salsa_only.ref.fasta
```

#### Changes list for dam
* Remove salsa scaffolds 556 and 608 from chr1 (small alignment)
* Remove salsa scaffold 815 from chr7 (small alignment)
* Only kept the first 50 mbp of scaffold_539 for chr13 (scaffold_539 shares sequence on chr15)
* Adjusted scaffold_539 for chr15 50484937-126484936
* Salsa scaffold_65 appears to be circular on chr10
* Salsa scaffold_1283 appears to be circular on chr27

```bash
# Adding the missing contigs
perl convert_agp_to_bed.pl dam_scaffolding_tech_mashmap.tab dam_scaffold_salsa.agp dam_scaffold_salsa.agp.bed

bedtools intersect -a dam_assembly_contig_lengths.bed -b dam_scaffold_salsa.agp.bed -v | perl -lane 'print "$F[0]"; system("samtools faidx ../bostaurus_brahma.reformat.fasta $F[0] >> dam_scaffold_salsa_only.missingctg.fasta");'

cat dam_scaffold_salsa_only.ref.fasta dam_scaffold_salsa_only.missingctg.fasta > dam_scaffold_salsa_only.complete.fasta

## Now to generate the sire assembly
perl compare_algns_per_chrassignment.pl -l sire_scaffs.salsa.list -f ../../ars_ucd_123/ARS-UCDv1.0.23.fasta.fai -t 50000 -o sire_scaffs_consensus.salsa.tab
# Oh man, it's worse than the previous one! Still going to try

perl -e '%row = ("salsa" => 3, "phase" => 8); <>; while(<>){chomp; @s = split(/\t/); $scaff = $s[$row{$s[1]} - 1]; @ss = split(/;/, $scaff); $p = 0; for(my $x = 0; $x < scalar(@ss); $x++){$j = $ss[$x]; $p++; $sname = $j =~ s/,.+$//; print "$s[0]\t0\t0\t$p\tA\t$j\t0\t0\t?\t$s[1]\n"; unless($x + 1 >= scalar(@ss)){$p++; print "$s[0]\t0\t0\t$p\tU\t100\tscaffold\tyes\tmap\n";}}}' < sire_scaffs_consensus.salsa.tab > sire_scaffs_consensus.salsa.pagp

perl prepare_master_fasta.pl sire_scaffs_consensus.salsa.pagp sire_scaffold_fastas.tab sire_scaffold_salsa

perl -ne '$_ =~ s/\?/+/g; print $_;' < sire_scaffold_salsa.agp > temp
mv temp sire_scaffold_salsa.agp
```

#### Changes list for Sire
* Removed scaffold_520 from chr1 (small alignment)
* The first 133986306 of scaffold_520 belongs on chr2
* The 133986306-253622838 range of scaffold_520 belongs on chr5
* Only the first 19792045 bases of scaffold_1039 belong on chr8
* NOTE: the middle of scaffold_1032 (34 mb - 45 mb) belong on chrs 28 and 24
* Adjusted the range of scaffold_1039 on chr12 to 19792045-107435571
* Removed scaffold_68 from chr15 (small alignment)
* The first 70494700 of scaffold_68 belongs on chr20
* Removed the first 70494700 of scaffold_68 from chr23
* Removed scaffold_1032 from chr28 (small alignment)

So that's ALLOT more editting needed! The Angus side is somehow more error prone from the Hi-C perspective? Very strange.

```bash
samtools faidx sire_scaffold_salsa.fasta
java -Xmx100g -jar CombineFasta.jar agp2fasta -f sire_scaffold_salsa.fasta -a sire_scaffold_salsa.agp -o sire_scaffold_salsa_only.ref.fasta

sbatch --nodes=1 --ntasks-per-node=5 --mem=25000 -p assemble1 --wrap="../../../binaries/mashmap-Linux64-v2.0/mashmap -r ../bostaurus_angus.reformat.fasta -q f1_sire_salsa.fasta -t 5 -f one-to-one -o f1_sire_salsa.fasta.mashmap"

ls f1_sire*mashmap > sire_mashmap.tab
perl convert_agp_to_bed.pl sire_mashmap.tab sire_scaffold_salsa.agp sire_scaffold_salsa.contigmap.bed

# There were too many agp entries to edit
```