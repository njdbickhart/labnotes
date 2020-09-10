# Gaur Genome Assembly analysis
---
*4/28/2020*

These are my notes on running analysis on the Gaur Genome Assembly ahead of publication.

## Table of Contents

## FRC Align

> Ceres: /project/gaur_genome_assembly/HenryDoorly_FRCalign

```bash
# Creating aligned bams for FRC analysis
cat /project/gaur_genome_assembly/HenryDoorly_gaur_illumina_data/input.fofn > gaur_reads_spreadsheet.tab

# removing those horrible pipes in fasta headers
perl -e 'while(<>){$_ =~ s/\|/_/g; print $_;}' < ../HenryDoorly_gaur_pacbio/UOA_Gaur_PB_HiC_cattle.fasta > UOA_Gaur_PB_HiC_cattle.reformat.fasta

sbatch --nodes=1 --mem=5000 --ntasks-per-node=1 -p priority -q msn --wrap="module load samtools; samtools faidx UOA_Gaur_PB_HiC_cattle.reformat.fasta"
sbatch --nodes=1 --mem=15000 --ntasks-per-node=2 -p priority -q msn --wrap="module load bwa; bwa index UOA_Gaur_PB_HiC_cattle.reformat.fasta"

# Aligning reads to the reference
sbatch --nodes=1 --mem=10000 --ntasks-per-node=1 --dependency=afterok:2868201 -p priority -q msn --wrap="python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b gaur_scaffolds -t gaur_reads_spreadsheet.tab -f /project/gaur_genome_assembly/HenryDoorly_FRCalign/UOA_Gaur_PB_HiC_cattle.reformat.fasta -q msn -p priority -m"
```

## PBJelly fiasco

PBjelly did it again! It looks like it added 7 megabases of sequence. I'm highly skeptical of this, so let's check it out and make sure the program isn't just adding random sequence here. 

> Ceres: /project/gaur_genome_assembly/HenryDoorly_gaur_asms

```bash
module load java/1.8.0_121 samtools bedtools/2.25.0
# What a mess... I need to align the old scaffolded assembly against the new to try to infer the gap coordinates
minimap2 -x asm5 UOA_Gaur_PB_HiC_Arrow.fasta ARS-UOA_Gaur_PB_Arrow2_HiC_PBJelly_Arrow.fasta > arrow_to_pbjelly_comparison.paf

# Now to hack out the unaligned regions using a combination of bedtools and perl
# This creates a bed with just the major chromosomes of the final assembly
perl -lane '@bsegs = split(/\|/, $F[0]); if(length($bsegs[0]) <= 2){print "$bsegs[0]\t$F[2]\t$F[3]\t$F[5]\t$F[7]\t$F[8]";}' < arrow_to_pbjelly_comparison.paf | bedtools sort -i stdin > arrow_to_pbjelly_comparison.chrs.bed

# I had to manually edit the file to remove obvious repeat alignments that were internal to other alignments
perl -lane 'print "$F[0]\t$F[1]\t$F[2]";' < arrow_to_pbjelly_comparison.chrs.bed | bedtools merge -i stdin > arrow_to_pbjelly_comparison.algngaps.bed

# Getting the alignment based gaps
perl -lane '@bsegs = split(/\|/, $F[0]); if(length($bsegs[0]) <= 4){print "$bsegs[0]\t1\t$F[1]";}' < ARS-UOA_Gaur_PB_Arrow2_HiC_PBJelly_Arrow.fasta.fai | bedtools sort -i stdin > ARS-UOA_Gaur_PB_Arrow2_HiC_PBJelly_Arrow.chrlens.bed

bedtools subtract -b arrow_to_pbjelly_comparison.algngaps.bed -a ARS-UOA_Gaur_PB_Arrow2_HiC_PBJelly_Arrow.chrlens.bed > arrow_to_pbjelly_comparison.algngaps.subt.bed


# Identifying gaps in the gap-filled assembly
java -Xmx10g -jar ~/rumen_longread_metagenome_assembly/binaries/GetMaskBedFasta/store/GetMaskBedFasta.jar -f ARS-UOA_Gaur_PB_Arrow2_HiC_PBJelly_Arrow.fasta -o ARS-UOA_Gaur_PB_Arrow2_HiC_PBJelly_Arrow.gaps.bed -s ARS-UOA_Gaur_PB_Arrow2_HiC_PBJelly_Arrow.gaps.stats
perl -lane '@bsegs = split(/\|/, $F[0]); print "$bsegs[0]\t$F[1]\t$F[2]";' < ARS-UOA_Gaur_PB_Arrow2_HiC_PBJelly_Arrow.gaps.bed | bedtools sort -i stdin > ARS-UOA_Gaur_PB_Arrow2_HiC_PBJelly_Arrow.gaps.rfmt.bed

bedtools intersect -a arrow_to_pbjelly_comparison.algngaps.subt.bed -b ARS-UOA_Gaur_PB_Arrow2_HiC_PBJelly_Arrow.gaps.rfmt.bed > arrow_to_pbjelly_comparison.algngaps.subt.nogaps.bed

cat arrow_to_pbjelly_comparison.algngaps.subt.nogaps.bed | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/bed_length_sum.pl
        Interval Numbers:       127
        Total Length:           5236
        Length Average:         41.2283464566929
        Length Median:          24
        Length Stdev:           31.0752677759104
        Smallest Length:        24
        Largest Length:         99

# Damn it!
# Didn't work
```

New tactic: use the pbjelly json information against it

> Ceres: /project/gaur_genome_assembly/HenryDoorly_gaur_PBJelly

```bash
python3 makeSenseOfJelly.py gap_fill_status.txt assembly gap_fill_status_table.tab


```

To parse it all in R

```R
data <- read.delim("gap_fill_status_table.tab", header=TRUE)
temp <- data[,c(1,2,3,4,5,6,7,8)]
temp <- na_if(temp, "Nope")
temp <- na_if(temp, "None")

temp$fillBases <- as.numeric(temp$fillBases)
temp$contribSeqs <- as.numeric(temp$contribSeqs)
temp$spanCount <- as.numeric(temp$spanCount)

temp %>% group_by(GapCat) %>% summarize(n = n(), avgFill = mean(fillBases, na.rm=TRUE), sumFill = sum(fillBases, na.rm=TRUE), avgContrib = mean(contribSeqs, na.rm=TRUE), avgSpanCnt = mean(spanCount, na.rm=TRUE))
  GapCat            n avgFill sumFill avgContrib avgSpanCnt
  <fct>         <int>   <dbl>   <dbl>      <dbl>      <dbl>
1 doubleextend    159    790.  125544       2.91       1
2 filled          269    502.  135095      32.4        5.51
3 minreadfail    1171      1     1171       1          1
4 nofillmetrics    13    NaN        0     NaN        NaN
5 overfilled      638    503.  321082      30.9        1
6 singleextend    462    679.  313852      16.8        1
```

## Pilon correction

> Ceres: /project/gaur_genome_assembly/HenryDoorly_gaur_pilon

```bash
module load bwa samtools
cat ../HenryDoorly_gaur_illumina_data/input.fofn > spreadsheet.tab

perl -ne 'if($_ =~ />/){chomp($_); @bsegs = split(/\|/, $_); print "$bsegs[0]\n";}else{print $_;}' < /project/gaur_genome_assembly/HenryDoorly_gaur_asms/UOA_Gaur_PB_HiC_cattle.fasta > ARS_UOA_Gaur_PB_HiC_cattle.rfmt.fasta

sbatch --nodes=1 --mem=15000 --ntasks-per-node=1 -p priority -q msn -t 1-0 --wrap="bwa index ARS_UOA_Gaur_PB_HiC_cattle.rfmt.fasta"

python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b henrydoorly -t spreadsheet.tab -f /project/gaur_genome_assembly/HenryDoorly_gaur_pilon/ARS_UOA_Gaur_PB_HiC_cattle.rfmt.fasta -m -p priority -q msn -e '2-0'

sbatch --dependency=afterok:3006206 --nodes=1 --mem=12000 --ntasks-per-node=2 -p priority -q msn -t 1-0 --wrap="python3 ~/python_toolchain/sequenceData/slurmPilonCorrectionPipeline.py -g /project/gaur_genome_assembly/HenryDoorly_gaur_pilon/ARS_UOA_Gaur_PB_HiC_cattle.rfmt.fasta -f henrydoorly/gaur/gaur.sorted.merged.bam -o pilonone -p priority -q msn -m -t '2-0'"

## The redo

sbatch --nodes=1 --mem=15000 --ntasks-per-node=1 -p priority -q msn -t 1-0 --wrap="bwa index /project/gaur_genome_assembly/ARS-UOA_Gaur_PB_HiC_Arrow.fasta"

sbatch --dependency=afterok:3009415 --nodes=1 --mem=10000 --ntasks-per-node=1 -p priority -q msn -t 1-0 --wrap="python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b henrydoorly -t spreadsheet.tab -f /project/gaur_genome_assembly/ARS-UOA_Gaur_PB_HiC_Arrow.fasta -m -p priority -q msn -e '2-0'"

sbatch --nodes=1 --mem=12000 --ntasks-per-node=2 -p priority -q msn -t 1-0 --wrap="python3 ~/python_toolchain/sequenceData/slurmPilonCorrectionPipeline.py -g /project/gaur_genome_assembly/ARS-UOA_Gaur_PB_HiC_Arrow.fasta -f henrydoorly/gaur/gaur.sorted.merged.bam -o pilonone -p priority -q msn -m -t '2-0'"

# I had to modify the scripts to increase the memory allocation and use a different version of pilon
for i in pilonone/scripts/*.sh; do echo $i; perl -lane 'if($_ =~ /--mem/){print "\#SBATCH --mem=150000";}elsif($_ =~ /module/){print "module load pilon/1.23 samtools";}elsif($_ =~ /^echo/ || $_ =~ /^java/){$_ =~ s/-Xmx10G/-Xmx150G/; $_ =~ s/pilon-1.22.jar/pilon-1.23.jar/; $_ =~ s/ARS-UOA_Gaur_PB_HiC_Arrow.fasta/pilon_temp.fa/; print $_;}else{print $_;}' < $i > temp; mv temp $i; done

sbatch --nodes=1 --mem=15000 --ntasks-per-node=1 -p priority -q msn -t 1-0 --wrap="module load bwa samtools; cat pilonone/*.fasta > ARS-UOA_Gaur_PB_HiC_Arrow_p1.fasta; bwa index ARS-UOA_Gaur_PB_HiC_Arrow_p1.fasta; samtools faidx ARS-UOA_Gaur_PB_HiC_Arrow_p1.fasta"

# Now time for the round 2 alignments
sbatch --nodes=1 --mem=10000 --ntasks-per-node=1 -p priority -q msn -t 1-0 --wrap="python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b p1hdoorly -t spreadsheet.tab -f /project/gaur_genome_assembly/HenryDoorly_gaur_pilon/ARS-UOA_Gaur_PB_HiC_Arrow_p1.fasta -m -p priority -q msn -e '2-0'"

# And the pilon alignment
sbatch --nodes=1 --mem=12000 --ntasks-per-node=2 -p priority -q msn -t 1-0 --wrap="python3 ~/python_toolchain/sequenceData/slurmPilonCorrectionPipeline.py -g /project/gaur_genome_assembly/HenryDoorly_gaur_pilon/ARS-UOA_Gaur_PB_HiC_Arrow_p1.fasta -f p1hdoorly/gaur/gaur.sorted.merged.bam -o pilontwo -p priority -q msn -m -t '2-0' -e 150000"

# Finally, reformatting the fasta for submission
python3 ~/python_toolchain/assembly/reformatPilonOutput.py -f pilontwo -o ARS-UOA_Gaur_PB_HiC_Arrow_p2.fasta -t 8
```

## Assembly QV estimation

I've written a pipeline that should help to automate this. Let's try queuing it up and crossing our fingers!

> Ceres: /lustre/project/gaur_genome_assembly/HenryDoorly_gaur_qv

```bash
module load miniconda/3.6

sbatch --nodes=1 --mem=5000 --ntasks-per-node=2 -p priority -q msn snakemake --cluster-config ~/python_toolchain/snakeMake/assemblyValidation/cluster.json --cluster "sbatch --nodes={cluster.nodes} --ntasks-per-node={cluster.ntasks-per-node} --mem={cluster.mem} --partition={cluster.partition} -q {cluster.qos} -o {cluster.stdout}" -p --jobs 250 -s ~/python_toolchain/snakeMake/assemblyValidation/assemblyValidation --use-conda

# This ideogram plotting script works:
python ~/python_toolchain/snakeMake/assemblyValidation/scripts/ideogram_feature_plot.py -b mapped/merged.bam -f calls/merged_frc.txt_Features.txt -o final/ideogram_frc_plot.pdf
```

## Bison assembly modifications

This is a separate, but smaller project that needs my attention. I am going to modify the Bison sire haplotig assembly based on Jonas' instructions.

> Ceres: /lustre/project/cattle_genome_assemblies/bison_x_simmental/bison_sire_fixes

```bash
module load samtools
ln -s ../Salsa_Bison_sire/scaffolds/scaffolds_FINAL.fasta sire_salsa_scaffolds.fasta

samtools faidx sire_salsa_scaffolds.fasta
perl -lane 'print "$F[0]\t1\t$F[1]\t$F[0]\t1\t+";' < sire_salsa_scaffolds.fasta.fai > sire_salsa_scaffolds_plan.bed

# Now to do manual edits before running my program!
sbatch --nodes=1 --mem=200000 --ntasks-per-node=10 -iority -q msn -t 1-0 --wrap="java -Xmx200g -jar /lustre/project/rumen_longread_metagenome_assembly/binaries/CombineFasta/store/CombineFasta.jar agp2fasta -f sire_salsa_scaffolds.fasta -b sire_salsa_scaffolds_plan.bed -o sire_salsa_jonas_fixes.fasta"

minimap2 -x asm5 sire_salsa_jonas_fixes.fasta.gz sire_salsa_scaffolds.fasta > sire_sire_checks.paf
```