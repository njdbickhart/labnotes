# Yaklander genome assembly analysis
---
*3/25/2019*

These are my notes on running some analysis on the Yaklander (Yak x Highland cross) triobinned assembly.

## Table of Contents


## Preparing the assembly fastas

I just need to gather the materials I need for the analysis.

> Ceres: /home/derek.bickharhth/cattle_genome_assemblies/yaklander

```bash
module load bwa
sbatch --nodes=1 --mem=14000 --ntasks-per-node=1 -p short --wrap="bwa index yaklander.dam.gapfilled.arrow2.fasta"
sbatch --nodes=1 --mem=14000 --ntasks-per-node=1 -p short --wrap="bwa index yaklander.sire.gapfilled.arrow2.fasta"
```


## Running the repeat analysis

OK, now that everything is ready, let's queue up RepeatMasker and generate the files I need. I've already done repeatmasking on the ARS-UCDv1.2 assembly and the UMD3.1 assembly, so those comparative repeat lengths can be added to the plots if needed for comparison.

> Ceres: /home/derek.bickharhth/cattle_genome_assemblies/yaklander
 
```bash
module unload perl/5.24.1
sbatch --nodes=1 --mem=40000 --ntasks-per-node=25 -p medium --wrap="/beegfs/project/rumen_longread_metagenome_assembly/binaries/RepeatMasker/RepeatMasker -pa 25 -q -species cow -no_is -gff yaklander.dam.gapfilled.arrow2.fasta"
sbatch --nodes=1 --mem=40000 --ntasks-per-node=25 -p medium --wrap="/beegfs/project/rumen_longread_metagenome_assembly/binaries/RepeatMasker/RepeatMasker -pa 25 -q -species cow -no_is -gff yaklander.sire.gapfilled.arrow2.fasta"

# I've combined the gap and repeat analysis into one rolling script
sbatch generate_repeat_counts.sh yaklander.dam.gapfilled.arrow2.fasta yakdam ../dominette/ARS-UCD1.2_Btau5.0.1Y/ARS-UCD1.2_Btau5.0.1Y.fa
sbatch generate_repeat_counts.sh yaklander.sire.gapfilled.arrow2.fasta yaksire ../dominette/ARS-UCD1.2_Btau5.0.1Y/ARS-UCD1.2_Btau5.0.1Y.fa

# Grrr! The "|arrow" delimiters strike again! I need to remove them 
module load bwa samtools
perl -e 'while(<>){if($_ =~ /^>/){$_ =~ s/\|arrow//g;} print $_;}' < yaklander.dam.gapfilled.arrow2.fasta > yaklander.dam.gapfilled.rfmt.fa
perl -e 'while(<>){if($_ =~ /^>/){$_ =~ s/\|arrow//g;} print $_;}' < yaklander.sire.gapfilled.arrow2.fasta > yaklander.sire.gapfilled.rfmt.fa
sbatch --nodes=1 --mem=12000 --ntasks-per-node=1 -p msn --wrap="bwa index yaklander.dam.gapfilled.rfmt.fa; samtools faidx yaklander.dam.gapfilled.rfmt.fa"
sbatch --nodes=1 --mem=12000 --ntasks-per-node=1 -p msn --wrap="bwa index yaklander.sire.gapfilled.rfmt.fa; samtools faidx yaklander.sire.gapfilled.rfmt.fa"

sbatch generate_repeat_counts.sh yaklander.dam.gapfilled.arrow2.fasta yakdam ../dominette/ARS-UCD1.2_Btau5.0.1Y/ARS-UCD1.2_Btau5.0.1Y.fa
sbatch generate_repeat_counts.sh yaklander.sire.gapfilled.arrow2.fasta yaksire ../dominette/ARS-UCD1.2_Btau5.0.1Y/ARS-UCD1.2_Btau5.0.1Y.fa

sbatch generate_repeat_counts.sh yaklander.dam.gapfilled.rfmt.fa yakdam ../dominette/ARS-UCD1.2_Btau5.0.1Y/ARS-UCD1.2_Btau5.0.1Y.fa
sbatch --dependency=afterok:659870 generate_repeat_counts.sh yaklander.sire.gapfilled.rfmt.fa yaksire ../dominette/ARS-UCD1.2_Btau5.0.1Y/ARS-UCD1.2_Btau5.0.1Y.fa


# Ack! I screwed up my script and read in the fasta file instead of the repeatmasker output!
perl -e '<>; <>; <>; while(<>){ $_ =~ s/^\s+//; @s = split(/\s+/); $orient = ($s[8] eq "+")? "+" : "-"; $qlen = $s[12] - $s[11]; $s[4] =~ s/\|arrow//g; print "$s[4]\t$s[5]\t$s[6]\t$orient\t$s[9]\t$s[10]\t$qlen\n";}' < yaklander.dam.gapfilled.arrow2.fasta.out > yakdam.repeat.bed
python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f yakdam.repeat.bed -c 5 -d '\t' -m > yakdam.rep.count.md
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); while(<IN>){chomp; @F = split(/\t/); $F[5] =~ s/\?//g; $F[5] =~ tr/\//_/; print "$F[5]\t$F[6]\t$ARGV[1]\n";} close IN;' yakdam.repeat.bed yakdam > yakdam.repeat.lens

perl -e '<>; <>; <>; while(<>){ $_ =~ s/^\s+//; @s = split(/\s+/); $orient = ($s[8] eq "+")? "+" : "-"; $qlen = $s[12] - $s[11]; $s[4] =~ s/\|arrow//g; print "$s[4]\t$s[5]\t$s[6]\t$orient\t$s[9]\t$s[10]\t$qlen\n";}' < yaklander.sire.gapfilled.arrow2.fasta.out > yaksire.repeat.bed
python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f yaksire.repeat.bed -c 5 -d '\t' -m > yaksire.rep.count.md
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); while(<IN>){chomp; @F = split(/\t/); $F[5] =~ s/\?//g; $F[5] =~ tr/\//_/; print "$F[5]\t$F[6]\t$ARGV[1]\n";} close IN;' yaksire.repeat.bed yaksire > yaksire.repeat.lens
```