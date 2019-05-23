# Brangus haploid assembly analysis
---
*11/3/2017*

These are my notes on the analysis of the haploid Brangus reference assembly. Paternal haplotypes are derived from an Angus and the maternal is from a Brahman.

## Table of Contents
* [Assembly download, prep and marker analysis](#prep)
	* [Identifying consecutive blocks of missing HD markers](#consecutive)
* [Hi-C scaffolding check, assignment and validation](#scaffcheck)
* [Repeatmasker analysis and putative active retroviral element](#retroviral)
* [Creating full assemblies based on scaffold winners](#scaffwinners)
	* [Changes list for sire](#changes)
	* [Changes list for dam](#changes)
	* [Tracking contig alignments and finding missing sequence](#tracking)
	* [Checking duplication after selective scaffold reduction](#dupcheck)
* [Generating salsa scaffolds and assessing scaffold contiguity](#salsa)
	* [Changes list for dam](#schanges)
	* [Changes list for sire](#schanges)
	* [Salsa only scaffolding statistics](#sstats)
	* [3ddna2 sire only assembly changes](#ustats)
	* [Updated scaffolding statistics](#ustats)
* [QV and other summary statistics](#qv)
* [CNV calls and statistics](#cnv)

<a name="prep"></a>
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
<a name="consecutive"></a>
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

<a name="scaffcheck"></a>
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
<a name="retroviral"></a>
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

<a name="scaffwinners"></a>
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
<a name="changes"></a>
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
<a name="tracking"></a>
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
<a name="dupcheck"></a>
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
<a name="salsa"></a>
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
<a name="schanges"></a>
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

# I modified the script to print out the "last" scaffold association so that bedtools doesn't explode due to unequal column lengths
perl -lane 'print "$F[0]\t1\t$F[1]";' < ../bostaurus_angus.reformat.fasta.fai > sire_assembly_contig_lengths.bed
bedtools intersect -a sire_assembly_contig_lengths.bed -b sire_scaffold_salsa.contigmap.bed -v | perl -lane 'print $F[0];' | sort | uniq | perl -lane 'print $F[0]; system("samtools faidx ../bostaurus_angus.reformat.fasta $F[0] >> sire_scaffold_salsa_only.missingctg.fasta");'

cat sire_scaffold_salsa_only.ref.fasta sire_scaffold_salsa_only.missingctg.fasta > sire_scaffold_salsa_only.complete.fasta
perl -ne '$_ =~ s/\?/+/g; print $_;' < sire_scaffold_3ddna2.agp > temp
mv temp sire_scaffold_3ddna2.agp


```
<a name="sstats"></a>
##### Salsa only scaffolding statistics

|Assembly | OriginalCtgs | UnscaffoldedCtgs | UnscaffoldCtgLen | ScaffoldCtgLen |
| :--- | ---: | ---: | ---: | ---: |
|Sire  | 1747 | 1477 | 87,976,369 | 2,494,659,122 |
|Dam   | 1585 | 1337 | 76,849,908 | 2,604,211,377 |

Serge suggested making just a 3ddnav2 assembly as that scaffolding required fewer manual edits than the salsa assembly.

```bash
sbatch --nodes=1 --ntasks-per-node=25 --mem=25000 -p assemble1 --wrap="../../../binaries/mashmap-Linux64-v2.0/mashmap -r ../bostaurus_angus.reformat.fasta -q f1_sire_3ddna_v2.fasta -t 25 -f one-to-one -o f1_sire_3ddnav2.fasta.mashmap"

perl compare_algns_per_chrassignment.pl -l sire_scaffs.3ddna2.list -f ../../ars_ucd_123/ARS-UCDv1.0.23.fasta.fai -t 50000 -o sire_scaffs_consensus.3ddna2.tab

perl -e '%row = ("3ddna2" => 3, "phase" => 8); <>; while(<>){chomp; @s = split(/\t/); $scaff = $s[$row{$s[1]} - 1]; @ss = split(/;/, $scaff); $p = 0; for(my $x = 0; $x < scalar(@ss); $x++){$j = $ss[$x]; $p++; $sname = $j =~ s/,.+$//; print "$s[0]\t0\t0\t$p\tA\t$j\t0\t0\t?\t$s[1]\n"; unless($x + 1 >= scalar(@ss)){$p++; print "$s[0]\t0\t0\t$p\tU\t100\tscaffold\tyes\tmap\n";}}}' < sire_scaffs_consensus.3ddna2.tab > sire_scaffs_consensus.3ddna2.pagp

perl prepare_master_fasta.pl sire_scaffs_consensus.3ddna2.pagp sire_scaffold_fastas.tab sire_scaffold_3ddna2
samtools faidx sire_scaffold_3ddna2.fasta
ls f1_sire*.mashmap > sire_mashmap.tab

java -Xmx100g -jar CombineFasta.jar agp2fasta -f sire_scaffold_3ddna2.fasta -a sire_scaffold_3ddna2.agp -o sire_scaffold_3ddna2_only.ref.fasta

perl convert_agp_to_bed.pl sire_mashmap.tab sire_scaffold_3ddna2.agp sire_scaffold_3ddna2.contigmap.bed
bedtools intersect -a sire_assembly_contig_lengths.bed -b sire_scaffold_3ddna2.contigmap.bed -v | perl -lane 'print $F[0];' | sort | uniq | perl -lane 'print $F[0]; system("samtools faidx ../bostaurus_angus.reformat.fasta $F[0] >> sire_scaffold_3ddna2_only.missingctg.fasta");'

cat sire_scaffold_3ddna2_only.ref.fasta sire_scaffold_3ddna2_only.missingctg.fasta > angus_3ddna2_v1.complete.fasta
```
<a name="ustats"></a>
#### 3ddna2 sire only assembly changes
* Removed HiC_scaffold_3 from the end of chr1 (small alignment)
* Removed HiC_scaffold_49 from the end of chr24 (small alignment)

##### Updated scaffolding statistics

|Assembly | Scaffolds | OriginalCtgs | UnscaffoldedCtgs | UnscaffoldCtgLen | ScaffoldCtgLen |
| :--- | :--- |---: | ---: | ---: | ---: |
|Sire  | Salsa| 1747 | 1477 | 87,976,369 | 2,494,659,122 |
|Sire  |3ddna2| 1747 | 1198 | 50,677,645 | 2,523,588,061 |
|Dam   | Salsa| 1585 | 1337 | 76,849,908 | 2,604,211,377 |


## Testing gene location assumptions

There is some ambiguity as to the sex chromosome PAR locations. I want to run a read depth analysis to check and confirm the status of some genes on the PAR. Here is the [SRA page](https://www.ncbi.nlm.nih.gov/Traces/study/?WebEnv=NCID_1_70820143_130.14.18.97_5555_1525359272_3835068197_0MetA0_S_HStore&query_key=22).

> assembler2: /mnt/nfs/nfs2/bickhart-users/cattle_asms/angus_x_brahman

```bash
for i in SRX3666413 SRX3666412  SRX3666411  SRX3666410  SRX3666409  SRX3666408  SRX3666407  SRX3666403  SRX3666451  SRX3666450  SRX3666449; do sbatch download_sra_align.sh $i; done
```
<a name="qv"></a>
## QV and other summary statistics

I had set this up a while back. Here are the commands I used to process the data.

```bash
module load java/1.8.0_121
sbatch --nodes=1 --mem=12000 --ntasks-per-node=1 -p short --wrap="bwa index bostaurus_angus.reformat.fasta"
sbatch --nodes=1 --mem=12000 --ntasks-per-node=1 -p short --wrap="bwa index bostaurus_brahma.reformat.fasta"
ls /scinet01/gov/usda/ars/scinet/project/cattle_genome_assemblies/angusxbrahman/illumina/*.gz > illumina_angusxbrahman.fastqs.tab
python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b angus -t illumina_angusxbrahman.fastqs.tab -f bostaurus_angus.reformat.fasta -m -p short
python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b brahman -t illumina_angusxbrahman.fastqs.tab -f bostaurus_brahma.reformat.fasta -m -p short
sbatch --nodes=1 --ntasks-per-node=4 --mem=20000 -p short --wrap="java -Xmx19g -jar ~/rumen_longread_metagenome_assembly/binaries/JaRMS/store/JaRMS.jar call -i angus/angusxbrahman/angusxbrahman.sorted.merged.bam -f bostaurus_angus.reformat.fasta -o bostaurus_angus.jarms -t 4"
sbatch --nodes=1 --ntasks-per-node=4 --mem=20000 -p short --wrap="java -Xmx19g -jar ~/rumen_longread_metagenome_assembly/binaries/JaRMS/store/JaRMS.jar call -i brahman/angusxbrahman/angusxbrahman.sorted.merged.bam -f bostaurus_brahma.reformat.fasta -o bostaurus_brahman.jarms -t 4"

# Restarting with the final, NCBI versions

sbatch --nodes=1 --mem=12000 --ntasks-per-node=1 -p short --wrap="bwa index bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil.fasta"
sbatch --nodes=1 --mem=12000 --ntasks-per-node=1 -p short --wrap="bwa index bostaurus_brahma_bionano_NCBI_full_corrected_gapfill_arrow_fil_withM.fasta"
python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b f_angus -t illumina_angusxbrahman.fastqs.tab -f bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil.fasta -m -p short
python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b f_brahman -t illumina_angusxbrahman.fastqs.tab -f bostaurus_brahma_bionano_NCBI_full_corrected_gapfill_arrow_fil_withM.fasta -m -p short

sbatch freebayes -C 2 -0 -O -q 20 -z 0.10 -E 0 -X -u -p 2 -F 0.75 -f bostaurus_brahma_bionano_NCBI_full_corrected_gapfill_arrow_fil_withM.fasta -v bostaurus_brahman_bionano_NCBI_full_corrected.freebayes.vcf f_brahman/angusxbrahman/angusxbrahman.sorted.merged.bam
sbatch freebayes -C 2 -0 -O -q 20 -z 0.10 -E 0 -X -u -p 2 -F 0.75 -f bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil.fasta -v bostaurus_angus_bionano_NCBI_full_corrected.freebayes.vcf f_angus/angusxbrahman/angusxbrahman.sorted.merged.bam

sbatch calculate_qv.sh bostaurus_brahman_bionano_NCBI_full_corrected.freebayes.vcf f_brahman/angusxbrahman/angusxbrahman.sorted.merged.bam bostaurus_brahman_bionano_NCBI_full_corrected.QV
sbatch calculate_qv.sh bostaurus_angus_bionano_NCBI_full_corrected.freebayes.vcf f_angus/angusxbrahman/angusxbrahman.sorted.merged.bam bostaurus_angus_bionano_NCBI_full_corrected.QV

sbatch --nodes=1 --ntasks-per-node=4 --mem=20000 -p short --wrap="java -Xmx19g -jar ~/rumen_longread_metagenome_assembly/binaries/JaRMS/store/JaRMS.jar call -i f_angus/angusxbrahman/angusxbrahman.sorted.merged.bam -f bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil.fasta -o angus_bionano_NCBI.recall -t 4 -m 10000000"
sbatch --nodes=1 --ntasks-per-node=4 --mem=20000 -p short --wrap="java -Xmx19g -jar ~/rumen_longread_metagenome_assembly/binaries/JaRMS/store/JaRMS.jar call -i f_brahman/angusxbrahman/angusxbrahman.sorted.merged.bam -f bostaurus_brahma_bionano_NCBI_full_corrected_gapfill_arrow_fil_withM.fasta -o brahman_bionano_NCBI.recall -t 4 -m 10000000"

# Running the reads on the ARS_UCDv1.2 reference
python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b ars_ucd -t illumina_angusxbrahman.fastqs.tab -f ../../dominette/ARS-UCD1.2_Btau5.0.1Y/ARS-UCD1.2_Btau5.0.1Y.fa -m -p short


#### FRC align ####
module load frc_align
sbatch --nodes=1 --mem=25000 --ntasks-per-node=3 -p short --wrap="FRC --pe-sam f_angus/angusxbrahman/angusxbrahman.sorted.merged.bam --output f_angus_frc_stats"
sbatch --nodes=1 --mem=25000 --ntasks-per-node=3 -p short --wrap="FRC --pe-sam f_brahman/angusxbrahman/angusxbrahman.sorted.merged.bam --output f_brahman_frc_stats"
sbatch --nodes=1 --mem=25000 --ntasks-per-node=3 -p short --wrap="FRC --pe-sam ars_ucd/angusxbrahman/angusxbrahman.sorted.merged.bam --output ars_ucd_frc_stats"
```
<a name="cnv"></a>
## CNV calls and statistics

Now I need to download the SRA datasets for CNV calling. I'll queue them up as sequential, simultaneous, tasks.

> Ceres: /home/derek.bickharhth/cattle_genome_assemblies/angusxbrahman/public

```bash
perl -lane '@bsegs = split(/_/, $F[2]); print join("\t", @bsegs);' < lloyd_master_pe > beef_panel_simple_accession.tab
module load sratoolkit/2.9.0

cat beef_panel_simple_accession.tab | cut -f1 | xargs -I {} sbatch --nodes=1 --mem=2000 --ntasks-per-node=1 -p short --wrap="fastq-dump.2 -I --gzip --split-files {}"

# Generating the spreadsheet I need for my script
perl -lane '@f = `ls $F[0]*.gz`; chomp(@f); print join("\t", @f) . "\t$F[0]\t$F[1]";' < beef_panel_simple_accession.tab > beef_panel_pipeline_spreadsheet.tab

# Grrr! The SRA read splitting doesn't work with BWA because of non-unique pair names! Fixing...
perl -lane 'print $F[2];' < beef_panel_pipeline_spreadsheet.tab | xargs -I {} sbatch --nodes=1 --mem=3000 --ntasks-per-node=1 -p short --wrap="python3 ~/python_toolchain/sequenceData/fixSRAFastqFiles.py -f {}_1.fastq.gz -r {}_2.fastq.gz -o {}.fmt -l {}.log"

perl -lane 'for($x = 0; $x < 2; $x++){$F[$x] =~ s/\_([12])/\.fmt\.$1/; $F[$x] =~ s/fastq/fq/;} print join("\t", @F);' < beef_panel_pipeline_spreadsheet.fullp.tab > beef_panel_pipeline_spreadsheet.fullp.rfmt.tab
```

Now I can do the alignments. I'm using these assembly fastas on the cluster:

* /home/derek.bickharhth/bostauruscnv/assembly/ARS-UCD1.2_Btau5.0.1Y.fa
* /home/derek.bickharhth/cattle_genome_assemblies/angusxbrahman/asms/bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil.fasta
* /home/derek.bickharhth/cattle_genome_assemblies/angusxbrahman/asms/bostaurus_brahma_bionano_NCBI_full_corrected_gapfill_arrow_fil_withM.fasta

```bash
# Angus
python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b angus_asm -t beef_panel_pipeline_spreadsheet.fullp.rfmt.tab -f /project/cattle_genome_assemblies/angusxbrahman/asms/bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil.fasta -m -p msn

python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b brahman_asm -t beef_panel_pipeline_spreadsheet.fullp.rfmt.tab -f /project/cattle_genome_assemblies/angusxbrahman/asms/bostaurus_brahma_bionano_NCBI_full_corrected_gapfill_arrow_fil_withM.fasta -m -p msn

python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b arsucd_asm -t beef_panel_pipeline_spreadsheet.fullp.rfmt.tab -f /beegfs/project/bostauruscnv/assembly/ARS-UCD1.2_Btau5.0.1Y.fa -m -p msn

mkdir angus_jarms; module load java/1.8.0_121; for i in angus_asm/*/*.merged.bam; do name=`echo $i | cut -d'/' -f2`; echo $name; sbatch --nodes=1 --ntasks-per-node=4 --mem=20000 -p short --wrap="java -Xmx19g -jar ~/rumen_longread_metagenome_assembly/binaries/JaRMS/store/JaRMS.jar call -i $i -f /project/cattle_genome_assemblies/angusxbrahman/asms/bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil.fasta -o $name.JaRMs -t 4 -m 10000000; mv $name.JaRMs* ./angus_jarms/"; done
mkdir angus_lumpy; module load bwa samtools; for i in angus_asm/*/*.merged.bam; do name=`echo $i | cut -d'/' -f2`; echo $name; sbatch --nodes=1 --ntasks-per-node=1 --mem=15000 -p short --wrap="lumpyexpress -B $i -o $name.vcf; mv $name.vcf ./angus_lumpy/;"; done

mkdir brahman_jarms; module load java/1.8.0_121; for i in brahman_asm/*/*.merged.bam; do name=`echo $i | cut -d'/' -f2`; echo $name; sbatch --nodes=1 --ntasks-per-node=4 --mem=20000 -p short --wrap="java -Xmx19g -jar ~/rumen_longread_metagenome_assembly/binaries/JaRMS/store/JaRMS.jar call -i $i -f /project/cattle_genome_assemblies/angusxbrahman/asms/bostaurus_brahma_bionano_NCBI_full_corrected_gapfill_arrow_fil_withM.fasta -o $name.JaRMs -t 4 -m 10000000; mv $name.JaRMs* ./brahman_jarms/"; done
mkdir brahman_lumpy; module load bwa samtools; for i in brahman_asm/*/*.merged.bam; do name=`echo $i | cut -d'/' -f2`; echo $name; sbatch --nodes=1 --ntasks-per-node=1 --mem=15000 -p short --wrap="lumpyexpress -B $i -o $name.vcf; mv $name.vcf ./brahman_lumpy/;"; done

mkdir arsucd_jarms; module load java/1.8.0_121; for i in arsucd_asm/*/*.merged.bam; do name=`echo $i | cut -d'/' -f2`; echo $name; sbatch --nodes=1 --ntasks-per-node=4 --mem=20000 -p short --wrap="java -Xmx19g -jar ~/rumen_longread_metagenome_assembly/binaries/JaRMS/store/JaRMS.jar call -i $i -f /beegfs/project/bostauruscnv/assembly/ARS-UCD1.2_Btau5.0.1Y.fa -o $name.JaRMs -t 4 -m 10000000; mv $name.JaRMs* ./arsucd_jarms/"; done
mkdir arsucd_lumpy; module load bwa samtools; for i in arsucd_asm/*/*.merged.bam; do name=`echo $i | cut -d'/' -f2`; echo $name; sbatch --nodes=1 --ntasks-per-node=1 --mem=15000 -p short --wrap="lumpyexpress -B $i -o $name.vcf; mv $name.vcf ./arsucd_lumpy/;"; done

# OK, now to convert to BEDPE
for i in *_lumpy/*.vcf; do echo $i; sbatch --nodes=1 --ntasks-per-node=1 --mem=5000 -p msn --wrap="/home/derek.bickharhth/lumpy-sv/scripts/vcfToBedpe -i $i -o $i.bedpe"; done

# First, let's try to generate counts for each event just to see if raw counts are less on the "native" assembly for each breed.
# I  need to generate a list for each file, but then I think that I can push the data through my tabfile column counter script
module load perl/5.24.1
for i in Angus Brahman Gelbvieh Hereford RedAngus Shorthorn Simmental; do for j in `seq 1 6`; do echo "$i$j"; done; done > combinations.list

mkdir raw_counts
for i in `cat combinations.list`; do echo $i; perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f angus_lumpy/$i.vcf.bedpe,brahman_lumpy/$i.vcf.bedpe,arsucd_lumpy/$i.vcf.bedpe -c 10 -e '#' -o raw_counts/$i.lumpy.counts.table; done

# I wrote a quick compilation script to generate raw counts for Lloyd
cd raw_counts/
perl calculate_averages.pl Angus
```

And here is the short script I used to tabulate the count information:

```perl
#!/usr/bin/perl
# This is a one-off script designed to process tabFileColumnCounter output from the CNV files

use strict;
my $usage = "perl $0 <prefix of file>\n";
chomp(@ARGV);
unless(scalar(@ARGV) == 1){
        print $usage;
        exit;
}

my @files = `ls $ARGV[0]*`;
chomp(@files);

my %data; # {SVtype}->{asm}->{file} = count
foreach my $f (@files){
        open(my $IN, "< $f");
        for(my $x = 0; $x < 5; $x++){
                <$IN>;
        }
        while(my $line = <$IN>){
                chomp $line;
                my @segs = split(/\t/, $line);
                for(my $x = 0; $x < scalar(@segs); $x += 2){
                        $data{$segs[$x]}->{"asm$x"}->{$f} = $segs[$x + 1] * 1;
                }
        }
        close $IN;
}

print "Asm\tSV\tsum\tavgnum\tstdev\n";
foreach my $svs (sort {$a cmp $b} keys(%data)){
        if($svs eq "" || length($svs) == 0){next;}
        foreach my $asm (sort {$a cmp $b} keys(%{$data{$svs}})){
                my $count = 0;
                my $sum = 0;
                foreach my $file (sort {$a cmp $b} keys(%{$data{$svs}->{$asm}})){
                        $count += 1;
                        $sum += $data{$svs}->{$asm}->{$file};
                }
                if($count == 0){
                        print "$asm\t$svs\t0\t0\t0\n";
                }else{
                        my $avg = $sum / $count;
                        if($count == 1){
                                $avg = sprintf("%.3f", $avg);
                                print "$asm\t$svs\t$sum\t$avg\t0\n";
                                next;
                        }
                        my $ss = 0;
                        foreach my $file (sort {$a cmp $b} keys(%{$data{$svs}->{$asm}})){
                                $ss += ($data{$svs}->{$asm}->{$file} - $avg)**2;
                        }
                        $avg = sprintf("%.3f", $avg);
                        my $stdev = sprintf("%.3f", sqrt($ss / ($count - 1)));
                        print "$asm\t$svs\t$sum\t$avg\t$stdev\n";
                }
        }
}
```

Now to package things up to send to Lloyd.

```bash
mkdir angus_cnvdata
mkdir brahman_cnvdata
mkdir arsucd_cnvdata

cp angus_lumpy/*.bedpe ./angus_cnvdata/
for i in brahman arsucd; do cp $i"_lumpy"/*.bedpe ./$i"_cnvdata"/; done

for i in angus brahman arsucd; do cp $i"_jarms"/*.levels ./$i"_cnvdata"/; done

for i in angus brahman arsucd; do tar -czvf $i"_cnvdata.tar.gz" $i"_cnvdata"; done
```

## VST analysis and gene coordinate conversion

I want to estimate the variance of individual copy number of a gene-by-gene basis and send a list of highly variable genes to Lloyd for analysis.

I will use the Angus and Brahman individuals as an outgroup for this comparison, and calculate the Vst on each respective assembly.

> Ceres: /home/derek.bickharhth/cattle_genome_assemblies/angusxbrahman/public

```bash
# Converting the gtf and NCBI files into bed files for coordinate intersections
module load bedtools/2.25.0

perl -e 'while(<>){chomp; if($_ =~ /^\#/){next;} @s = split(/\t/); @bsegs = split(/\"/, $s[8]); if($s[2] =~ /gene/){print "$s[0]\t$s[3]\t$s[4]\t$bsegs[1]\n";}}' < sire.UOA_angus_1.96.gtf | bedtools sort -i stdin > sire.UOA_angus_1.96.gene.bed
perl -e 'while(<>){chomp; if($_ =~ /^\#/){next;} @s = split(/\t/); @bsegs = split(/\"/, $s[8]); if($s[2] =~ /gene/){print "$s[0]\t$s[3]\t$s[4]\t$bsegs[1]\n";}}' < dam.UOA_brahman_1.96.gtf | bedtools sort -i stdin > dam.UOA_brahman_1.96.gene.bed

dos2unix ncbi_ars_ucd_1.2_refseq.txt
perl -e '<>; while(<>){chomp; @s = split(/\t/); if($s[12] eq "" || $s[13] eq ""){next;} print "$s[10]\t$s[12]\t$s[13]\t$s[5]\n";}' < ncbi_ars_ucd_1.2_refseq.txt | bedtools sort -i stdin > ncbi_ars_ucd_1.2_refseq.bed

echo -e 'dam.UOA_brahman_1.96.gene.bed\tensembl' > dam.UOA_brahman_1.96.gene.list
echo -e 'sire.UOA_angus_1.96.gene.bed\tensembl' > sire.UOA_angus_1.96.gene.list
echo -e 'ncbi_ars_ucd_1.2_refseq.bed\trefseq' > ncbi_ars_ucd_1.2_refseq.list

module load java/1.8.0_121
sbatch annotateScript.sh angus sire.UOA_angus_1.96.gene.list angus.cnvrs
sbatch annotateScript.sh brahman dam.UOA_brahman_1.96.gene.list brahman.cnvrs
sbatch annotateScript.sh arsucd ncbi_ars_ucd_1.2_refseq.list arsucd.cnvrs

head -n 1 brahman.cnvrs_windows_ensembl.tab | perl -lane 'for($x = 6; $x < scalar(@F); $x++){print "$F[$x]\t1";}' > pop_list_base.tab

cp pop_list_base.tab brahman_pop_list.tab
cp pop_list_base.tab angus_pop_list.tab
cp pop_list_base.tab arsucd_pop_list.tab

python3 ~/python_toolchain/sequenceData/calculateVstDifferences.py -p brahman_pop_list.tab -c brahman.cnvrs_windows_ensembl.tab -o brahman.cnvrs_windows_vst_genes.bed
perl calculate_vst_differences_cn.pl -c brahman.cnvrs_windows_ensembl.tab -p brahman_pop_list.tab -o brahman.vst_test.bed

# I just updated the script to print out a melted file for plotting VST stats
# I also added a minimum filter for differences in CN count between averages
python3 ~/python_toolchain/sequenceData/calculateVstDifferences.py -p brahman_pop_list.tab -c brahman.cnvrs_windows_ensembl.tab -o brahman.cnvrs_windows_vst_genes.bed -m brahman.cnvrs_windows_vst_genes.melt -g 1.5

python3 ~/python_toolchain/sequenceData/calculateVstDifferences.py -p angus_pop_list.tab -c angus.cnvrs_windows_ensembl.tab -o angus.cnvrs_windows_vst_genes.bed -m angus.cnvrs_windows_vst_genes.melt -g 1.0
Melt    ENSBIXG00000012150      0.21302147685415118     2.8504814814814807
Melt    ENSBIXG00000001131      0.22260541603045392     4.582777777777775
Dealt with 805 null fields

python3 ~/python_toolchain/sequenceData/calculateVstDifferences.py -p arsucd_pop_list.tab -c arsucd.cnvrs_windows_refseq.tab -o arsucd.cnvrs_windows_vst_genes.bed -m arsucd.cnvrs_windows_vst_genes.melt -g 1.0
Melt    LOC112444932    0.25849345715440125     7.64288
Melt    LOC112444962    0.3087315274489216      10.855250909090909
Melt    LOC783362       0.4513821996097366      2.614326599326601
Melt    LOC112445156    0.40810151414518686     3.878969696969696
Melt    LOC112445196    0.3437332492268294      3.1218956228956225
Melt    TACR1   0.3711020323483955      1.1097407407407403
Melt    PPM1B   0.2389241763285037      1.3348316498316501
Melt    LOC100848847    0.28723153448177546     2.887454545454545
Melt    LOC781990       0.4349602508760344      1.3241548821548819
Melt    LOC112441660    0.4062663573458778      1.2487272727272725
Melt    LOC112442670    0.21099051517668033     1.3915999999999995
Melt    GALNT3  0.3525093476571583      1.0180290909090908
Melt    LOC112448013    0.33571454830575703     1.3718145454545456
Melt    TM2D3   0.3619254154831361      2.2063272727272722
Dealt with 13 null fields

# The ARS-UCD data is misbehaving because of my selection of both brahman and angus as the VST grouping. Let's make things completely consistent and use the same populations in each
# the Brahman file will remain the same, but the Angus and ARS-UCD files will be replaced
python3 ~/python_toolchain/sequenceData/calculateVstDifferences.py -p brahman_pop_list.tab -c angus.cnvrs_windows_ensembl.tab -o angus.cnvrs_windows_vst_genes.tvi.bed -m angus.cnvrs_windows_vst_genes.tvi.melt -g 1.5

python3 ~/python_toolchain/sequenceData/calculateVstDifferences.py -p brahman_pop_list.tab -c arsucd.cnvrs_windows_refseq.tab -o arsucd.cnvrs_windows_vst_genes.tvi.bed -m arsucd.cnvrs_windows_vst_genes.tvi.melt -g 1.5
```

And here is the script I'm using to call the variants.

```bash
#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=10000
#SBATCH -p msn
# $1 = type
# $2 = annotation_db file
# $3 = output basename

lumpy=${1}"_lumpy"
jarms=${1}"_jarms"

cndata=${1}".cnlist"
cnvs=${1}".cnvs"

# Process bedpe file
for i in `ls $lumpy/*.bedpe`
do
        grep -v 'BND' $i | perl -e '<>; while(<>){chomp; @s = split(/\t/); print "$s[0]\t$s[1]\t$s[5]\n";' > $i.bed
        python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f $i -c 10 -d '\t' > $i.stats
done

# Create cnv and cn lists
for i in `ls $lumpy/*.bed`
do
        sample=`basename $i | cut -d'.' -f1`
        echo -e "$i\t$sample"
done > $cnvs

for i in `ls $jarms/*.bed.levels`
do
        sample=`basename $i | cut -d'.' -f1`
        echo -e "$i\t$sample"
done > $cndata

java -Xmx10g -jar ~/rumen_longread_metagenome_assembly/binaries/AnnotateUsingGenomicInfo/store/AnnotateUsingGenomicInfo.jar -d $2 -i $cnvs -c $cndata -o $3 -t
```

Now let's plot some of the VST-melted data.

> F:/SharedFolders/brangus_assembly/vst_analysis

```R
setwd("F:/SharedFolders/brangus_assembly/vst_analysis/")
library(ggplot2)
library(dplyr)
library(readr)
angusasm <- read.delim("angus.cnvrs_windows_vst_genes.melt", header=TRUE, sep = "\t")
angusasm$Pop <- as.factor(angusasm$Pop)
ggplot(angusasm, aes(y=CN, x=Gene, fill=Pop)) + geom_boxplot() + scale_fill_brewer(palette="Dark2") + theme_bw() + coord_flip() + geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.75)) + labs(title="Angus ASM Angus WGS CN differences", x="Gene (Vst)", y="Copy Number")
dev.copy2pdf(file="angus_asm_top_vst_genes.pdf", useDingbats=FALSE)

brahmanasm <- read.delim("brahman.cnvrs_windows_vst_genes.melt", header=TRUE, sep="\t")
brahmanasm$Pop <- as.factor(brahmanasm$Pop)
ggplot(brahmanasm, aes(y=CN, x=Gene, fill=Pop)) + geom_boxplot() + scale_fill_brewer(palette="Dark2") + theme_bw() + coord_flip() + geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.75), dotsize=0.25) + labs(title="Brahman ASM Brahman WGS CN differences", x="Gene (Vst)", y="Copy Number")
dev.copy2pdf(file="brahman_asm_top_vst_genes.pdf", useDingbats=FALSE)

arsucd <- read.delim("arsucd.cnvrs_windows_vst_genes.melt", header=TRUE, sep="\t")
arsucd$Pop <- as.factor(arsucd$Pop)
ggplot(arsucd, aes(y=CN, x=Gene, fill=Pop)) + geom_boxplot() + scale_fill_brewer(palette="Dark2") + theme_bw() + coord_flip() + geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.75), dotsize=0.25) + labs(title="ARSUCD1.2 ASM Brahman and Angus WGS CN differences", x="Gene (Vst)", y="Copy Number")
dev.copy2pdf(file="arsucd_asm_top_vst_genes.pdf", useDingbats=FALSE)

# Now to calculate the stats on the larger datasets
# Angus first
data <- read.delim("angus.cnvrs_windows_vst_genes.bed", header=FALSE)
colnames(data) <- c("chr", "start", "end", "Vst", "Gene")
data <- data %>% mutate(Len = end - start)
data <- data %>% mutate(Outlier = ifelse(Vst > 0.35, as.character(Gene), ""))
ggplot(data, aes(x=chr, y=Vst, color=Vst)) + geom_jitter() + scale_color_gradient(low="blue", high="red") + theme_bw() + geom_text(aes(label=Outlier), na.rm = TRUE, hjust = -0.3, size = 1.5) + labs(title="Angus ASM Vst Gene distribution", x="Chromosome", y="Vst")
dev.copy2pdf(file="angus_chr_vst_plot.pdf", useDingbats=FALSE)

write_tsv(print(data %>% group_by(chr) %>% summarize(numGenes = n(), geneLen10 = quantile(Len, 0.1), geneLenAvg = mean(Len), geneLen90 = quantile(Len, 0.9), VstAvg = mean(Vst), VstStd = sd(Vst), VstMax = max(Vst)), n = Inf), path="angus_vst_summary_stats.tab", quote_escape="none")

# Now Brahman
data <- read.delim("brahman.cnvrs_windows_vst_genes.bed", header=FALSE)
colnames(data) <- c("chr", "start", "end", "Vst", "Gene")
data <- data %>% mutate(Len = end - start, Outlier = ifelse(Vst > 0.5, as.character(Gene), ""))
ggplot(data, aes(x=chr, y=Vst, color=Vst)) + geom_jitter() + scale_color_gradient(low="blue", high="red") + theme_bw() + geom_text(aes(label=Outlier), na.rm = TRUE, hjust = -0.3, size = 1.5) + labs(title="Brahman ASM Vst Gene distribution", x="Chromosome", y="Vst")
dev.copy2pdf(file="brahman_chr_vst_plot.pdf", useDingbats=FALSE)

write_tsv(print(data %>% group_by(chr) %>% summarize(numGenes = n(), geneLen10 = quantile(Len, 0.1), geneLenAvg = mean(Len), geneLen90 = quantile(Len, 0.9), VstAvg = mean(Vst), VstStd = sd(Vst), VstMax = max(Vst)), n = Inf), path="brahman_vst_summary_stats.tab", quote_escape="none")

# Finally, ARS-UCD
data <- read.delim("arsucd.cnvrs_windows_vst_genes.bed", header=FALSE)
colnames(data) <- c("chr", "start", "end", "Vst", "Gene")
data <- data %>% mutate(Len = end - start, Outlier = ifelse(Vst > 0.3, as.character(Gene), ""))
ggplot(data, aes(x=chr, y=Vst, color=Vst)) + geom_jitter() + scale_color_gradient(low="blue", high="red") + theme_bw() + geom_text(aes(label=Outlier), na.rm = TRUE, position = position_dodge(width = 0.5), hjust = -0.2, size = 1.5) + labs(title="ARSUCD ASM Vst Gene distribution", x="Chromosome", y="Vst")
dev.copy2pdf(file="arsucd_chr_vst_plot.pdf", useDingbats=FALSE)

write_tsv(print(data %>% group_by(chr) %>% summarize(numGenes = n(), geneLen10 = quantile(Len, 0.1), geneLenAvg = mean(Len), geneLen90 = quantile(Len, 0.9), VstAvg = mean(Vst), VstStd = sd(Vst), VstMax = max(Vst)), n = Inf), path="arsucd_vst_summary_stats.tab", quote_escape="none")

# Generating plots for the taurine vs indicine analysis
data <- read.delim("angus.cnvrs_windows_vst_genes.tvi.bed", header=FALSE)
colnames(data) <- c("chr", "start", "end", "Vst", "Gene")
data <- data %>% mutate(Len = end - start, Outlier = ifelse(Vst > 0.5, as.character(Gene), ""))
ggplot(data, aes(x=chr, y=Vst, color=Vst)) + geom_jitter() + scale_color_gradient(low="blue", high="red") + theme_bw() + geom_text(aes(label=Outlier), na.rm = TRUE, hjust = -0.3, size = 1.5) + labs(title="Angus ASM Vst Gene distribution", x="Chromosome", y="Vst")
dev.copy2pdf(file="angus_chr_vst_plot_tvi.pdf", useDingbats=FALSE)

write_tsv(print(data %>% group_by(chr) %>% summarize(numGenes = n(), geneLen10 = quantile(Len, 0.1), geneLenAvg = mean(Len), geneLen90 = quantile(Len, 0.9), VstAvg = mean(Vst), VstStd = sd(Vst), VstMax = max(Vst)), n = Inf), path="angus_vst_summary_stats_tvi.tab", quote_escape="none")


data <- read.delim("arsucd.cnvrs_windows_vst_genes.tvi.bed", header=FALSE)
colnames(data) <- c("chr", "start", "end", "Vst", "Gene")
data <- data %>% mutate(Len = end - start, Outlier = ifelse(Vst > 0.5, as.character(Gene), ""))
ggplot(data, aes(x=chr, y=Vst, color=Vst)) + geom_jitter() + scale_color_gradient(low="blue", high="red") + theme_bw() + geom_text(aes(label=Outlier), na.rm = TRUE, hjust = -0.3, size = 1.5) + labs(title="ARS-UCDv1.2 ASM Vst Gene distribution", x="Chromosome", y="Vst")
dev.copy2pdf(file="arsucd_chr_vst_plot_tvi.pdf", useDingbats=FALSE)

write_tsv(print(data %>% group_by(chr) %>% summarize(numGenes = n(), geneLen10 = quantile(Len, 0.1), geneLenAvg = mean(Len), geneLen90 = quantile(Len, 0.9), VstAvg = mean(Vst), VstStd = sd(Vst), VstMax = max(Vst)), n = Inf), path="arsucd_vst_summary_stats_tvi.tab", quote_escape="none")


angusasm <- read.delim("angus.cnvrs_windows_vst_genes.tvi.melt", header=TRUE, sep = "\t")
angusasm$Pop <- as.factor(angusasm$Pop)
ggplot(angusasm, aes(y=CN, x=Gene, fill=Pop)) + geom_boxplot() + scale_fill_brewer(palette="Dark2") + theme_bw() + coord_flip() + geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.75), dotsize=0.25) + labs(title="Angus ASM Angus WGS CN differences", x="Gene (Vst)", y="Copy Number")
dev.copy2pdf(file="angus_asm_top_vst_genes.tvi.pdf", useDingbats=FALSE)

arsucd <- read.delim("arsucd.cnvrs_windows_vst_genes.tvi.melt", header=TRUE, sep="\t")
arsucd$Pop <- as.factor(arsucd$Pop)
ggplot(arsucd, aes(y=CN, x=Gene, fill=Pop)) + geom_boxplot() + scale_fill_brewer(palette="Dark2") + theme_bw() + coord_flip() + geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.75), dotsize=0.25) + labs(title="ARSUCD1.2 ASM Brahman and Angus WGS CN differences", x="Gene (Vst)", y="Copy Number")
dev.copy2pdf(file="arsucd_asm_top_vst_genes.tvi.pdf", useDingbats=FALSE)
```

## Coordinate liftover and CNV comparison analysis

For future reference, here are the liftover chain file creation scripts we're going to use for the analysis:

#### generateLiftoverChain.sh

```bash
#!/usr/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100000

#$1 = reference fasta
#$2 = query fasta
#$3 = work folder

module load minimap2/2.6

converter=/beegfs/project/rumen_longread_metagenome_assembly/binaries/kentUtils/bin/linux.x86_64/sam2psl.py
kentFolder=/beegfs/project/rumen_longread_metagenome_assembly/binaries/kentUtils/bin/linux.x86_64/

name1=`basename $1 | cut -d'.' -f1`
name2=`basename $2 | cut -d'.' -f1`
twoBit1=${3}/${name1}.blat.2bit
twoBit2=${3}/${name2}.blat.2bit


minimap2 -k 19 -w 19 -d $1.mmi $1


echo "aligning ${name2} to ${name1}"

minimap2 -ax asm5 $1.mmi $2 > ${3}/${name2}_to_${name1}.sam

echo "converting sam to psl"

python $converter -i ${3}/${name2}_to_${name1}.sam -o ${3}/${name2}_to_${name1}_mmap.psl

echo "done aligning"

# Generate two bit files and info
echo "Generating TwoBit indicies"

${kentFolder}/faToTwoBit $1 $twoBit1
${kentFolder}/twoBitInfo $twoBit1 ${twoBit1}.info
${kentFolder}/faToTwoBit $2 $twoBit2
${kentFolder}/twoBitInfo $twoBit2 ${twoBit2}.info

echo "AxtChain file creation"
${kentFolder}/axtChain -linearGap=medium -psl ${3}/${name2}_to_${name1}_mmap.psl ${twoBit2} ${twoBit1} ${3}/${name2}_to_${name1}.chain

echo "Sorting Chain file"
${kentFolder}/chainSort ${3}/${name2}_to_${name1}.chain ${3}/${name2}_to_${name1}.sorted.chain

echo "Creating Net file"
${kentFolder}/chainNet ${3}/${name2}_to_${name1}.sorted.chain ${twoBit2}.info ${twoBit1}.info ${3}/${name2}_to_${name1}.net /dev/null

echo "Creating liftover file ${3}/${name2}_to_${name1}.liftover.chain"
${kentFolder}/netChainSubset ${3}/${name2}_to_${name1}.net ${3}/${name2}_to_${name1}.sorted.chain ${3}/${name2}_to_${name1}.liftover.chain
```

Now to start the liftover process by generating the conversion files.

> Ceres: /home/derek.bickharhth/cattle_genome_assemblies/angusxbrahman/public

```bash
mkdir liftoverChains; sbatch -p msn ~/bin/generateLiftoverChain.sh /beegfs/project/cattle_genome_assemblies/dominette/ARS-UCD1.2_Btau5.0.1Y/ARS-UCD1.2_Btau5.0.1Y.fa /beegfs/project/cattle_genome_assemblies/angusxbrahman/asms/bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil.fasta liftoverChains

sbatch -p msn ~/bin/generateLiftoverChain.sh /beegfs/project/cattle_genome_assemblies/dominette/ARS-UCD1.2_Btau5.0.1Y/ARS-UCD1.2_Btau5.0.1Y.fa /beegfs/project/cattle_genome_assemblies/angusxbrahman/asms/bostaurus_brahma_bionano_NCBI_full_corrected_gapfill_arrow_fil_withM.fasta liftoverChains/

# And reversed
sbatch -p msn ~/bin/generateLiftoverChain.sh /beegfs/project/cattle_genome_assemblies/angusxbrahman/asms/bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil.fasta /beegfs/project/cattle_genome_assemblies/dominette/ARS-UCD1.2_Btau5.0.1Y/ARS-UCD1.2_Btau5.0.1Y.fa liftoverChains

sbatch -p msn ~/bin/generateLiftoverChain.sh /beegfs/project/cattle_genome_assemblies/angusxbrahman/asms/bostaurus_brahma_bionano_NCBI_full_corrected_gapfill_arrow_fil_withM.fasta /beegfs/project/cattle_genome_assemblies/dominette/ARS-UCD1.2_Btau5.0.1Y/ARS-UCD1.2_Btau5.0.1Y.fa liftoverChains/


# OK, now to convert my annotation files into bed files
ls *.cnvrs_regions.tab
angus.cnvrs_regions.tab  arsucd.cnvrs_regions.tab  brahman.cnvrs_regions.tab

for i in *.cnvrs_regions.tab; do name=`echo $i | cut -d'.' -f1`; echo $name; perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); <IN>; while(<IN>){chomp; @s = split(/\t/); if($s[2] < 0 ){$s[2] = 1;} print "$s[1]\t$s[2]\t$s[3]\t$ARGV[1]\t$s[0]\t$s[4]\n";} close IN;' $i $name > $name.cnvrs_regions.bed; done

~/rumen_longread_metagenome_assembly/binaries/kentUtils/bin/linux.x86_64/liftOver angus.cnvrs_regions.bed liftoverChains/bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil_to_ARS-UCD1.liftover.chain angus.cnvrs_regions.arsucd.bed angus.cnvrs_regions.arsucd.unmapped
~/rumen_longread_metagenome_assembly/binaries/kentUtils/bin/linux.x86_64/liftOver brahman.cnvrs_regions.bed liftoverChains/bostaurus_brahma_bionano_NCBI_full_corrected_gapfill_arrow_fil_withM_to_ARS-UCD1.liftover.chain brahman.cnvrs_regions.arsucd.bed brahman.cnvrs_regions.arsucd.unmapped

```

#### Results of the liftover

| Assembly | StartingCNVRs | LiftedOver | Unmapped |
| :------- | ------------: | ---------: | -------: |
| Angus    |       1416    |    1313    |   103    |
| Brahman  |       2147    |    2021    |   126    |

Now to try to create upSet plots using the three assemblies. We're going to try a raw overlap first and see which regions are unique only to ARS-UCDv1.2. In our ideal model, homozygous SVs should be present in all three. Hets should be in ARS-UCDV1.2 and the parental assembly. SVs only in the parental genome are likely compression errors, and SVs only in ARS-UCDV1.2 are likely the same (or misassemblies!)

```bash
module load bedtools/2.25.0
# Let's convert them all to a format that Bedtools will recognize

for i in angus.cnvrs_regions.arsucd.bed brahman.cnvrs_regions.arsucd.bed arsucd.cnvrs_regions.bed; do name=`echo $i | cut -d'.' -f1`; echo $name;

perl -lane 'print "$F[0]\t$F[1]\t$F[2]\t$F[3]";' < angus.cnvrs_regions.arsucd.bed | bedtools sort -i stdin > angus.cnvrs_regions.arsucd.sort.bed
perl -lane 'print "$F[0]\t$F[1]\t$F[2]\t$F[3]";' < brahman.cnvrs_regions.arsucd.bed | bedtools sort -i stdin > brahman.cnvrs_regions.arsucd.sort.bed
perl -lane 'print "$F[0]\t$F[1]\t$F[2]\t$F[3]";' < arsucd.cnvrs_regions.bed | bedtools sort -i stdin > arsucd.cnvrs_regions.sort.bed

cat angus.cnvrs_regions.arsucd.sort.bed brahman.cnvrs_regions.arsucd.sort.bed arsucd.cnvrs_regions.sort.bed | bedtools sort -i stdin | bedtools merge -c 4 -o collapse -delim ';' -i stdin > combined_merger.merge.bed

# A couple of chromosomes have huge overlap regions. Still, let's count the naieve overlap
perl -lane '%h = (); foreach $x (split(/\;/, $F[3])){$h{$x} = 1;} $e = join(";", sort{$a cmp $b} keys(%h)); print "$F[0]\t$F[1]\t$F[2]\t$e";' < combined_merger.merge.bed > combined_merger.merge.distinct.bed

python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f combined_merger.merge.distinct.bed -c 3 -d '\t' -m
```

|Entry                | Value|
|:--------------------|-----:|
|brahman              |  1733|
|arsucd               |  1183|
|angus                |  1028|
|angus;arsucd;brahman |    25|
|arsucd;brahman       |    24|
|angus;arsucd         |    16|
|angus;brahman        |    14|


```bash
# Ah, I see the problem. Some of the regions are lumpy errors reporting megabase sized CNVs
mkdir cnvoverlap
cp *.sort.bed ./cnvoverlap/
cd cnvoverlap/

for i in *.sort.bed; do echo $i; perl -lane 'print ($F[2] - $F[1]);' < $i | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl ; done
angus.cnvrs_regions.arsucd.sort.bed
total   1313
Sum:    181469943
Minimum 291
Maximum 88557271
Average 138210.162224
Median  1335
Standard Deviation      2817030.832286
Mode(Highest Distributed Value) 654

arsucd.cnvrs_regions.sort.bed
total   1342
Sum:    454548690
Minimum 168
Maximum 86128770
Average 338709.903130
Median  1402
Standard Deviation      4505098.254899
Mode(Highest Distributed Value) 903

brahman.cnvrs_regions.arsucd.sort.bed
total   2021
Sum:    336551683
Minimum 90
Maximum 107020199
Average 166527.304800
Median  1316
Standard Deviation      3217020.315354
Mode(Highest Distributed Value) 1319

# OK, let's remove all CNVRs larger than 1 megabase
for i in *.sort.bed; do name=`echo $i | cut -d '.' -f1`; perl -lane 'if($F[2] - $F[1] > 1000000){next;}else{print $_;}' < $i > $name.cnvrs_filt.sort.bed; done

cat *.cnvrs_filt.sort.bed | bedtools sort -i stdin | bedtools merge -c 4 -o distinct -delim ';' -i stdin > combined_merger_filt.merge.bed
python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f combined_merger_filt.merge.bed -c 3 -d '\t' -m
```

|Entry                | Value|
|:--------------------|-----:|
|brahman              |  1953|
|arsucd               |  1277|
|angus                |  1260|
|arsucd;brahman       |    22|
|angus;brahman        |    19|
|angus;arsucd         |    15|
|angus;arsucd;brahman |     6|

OK, I think that this is more reasonable. Let's plot this for R summary figures

```bash
perl -e '$cnv =  1; while(<>){chomp; @s = split(/\t/); @bsegs = split(/\;/, $s[3]); foreach my $b (@bsegs){open(OUT, ">> $b.upset.list"); print OUT "$cnv\n"; close OUT;} $cnv++;}' < combined_merger_filt.merge.bed
```

Now for the upSet plot.

> F:/SharedFolders/brangus_assembly/vst_analysis/

```R
library(UpSetR)

angus <- read.delim("angus.upset.list", header=FALSE)
brahman <- read.delim("brahman.upset.list", header=FALSE)
arsucd <- read.delim("arsucd.upset.list", header=FALSE)

angus <- as.vector(angus$V1)
brahman <- as.vector(brahman$V1)
arsucd <- as.vector(arsucd$V1)

input <- list(angus = angus, brahman = brahman, arsucd = arsucd)
upset(fromList(input), order.by = "freq")
dev.copy2pdf(file="upset_cnvr_intersection_plot.pdf", useDingbats=FALSE)
```