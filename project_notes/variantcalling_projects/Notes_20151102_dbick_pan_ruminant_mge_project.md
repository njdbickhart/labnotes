# Ruminant variant genome-wide comparison project
---
*11/2/2015*

These are my notes on the pan ruminant comparative genomics project I'm hoping to present at PAG. 

Let's start by first generating the BAMs I need. Buffalo will take some time so let's get that started.

> Blade2: /mnt/iscsi/md3200i_4/dbickhart/Buffalo/new_umd3

```bash
# Generating the spreadsheet
ls ../../../schroeder/Project_WaterBuffalo/full/ | perl -e '%h; while(<>){chomp; @s = split(/_/); $h{$s[0]}->{$s[4]}->{$s[3]} = $_;} foreach $an (keys(%h)){foreach $num (keys(%{$h{$an}})){ $fq1 = $h{$an}->{$num}->{"R1"}; $fq2 = $h{$an}->{$num}->{"R2"}; print "$fq1\t$fq2\t$an\t$an\n";}}' > new_umd3_spreadsheet.tab
```

Damn, I can't install my version of Perl on the system! A prior local-lib installation is interfering with everything and I don't have time to correct it. I'm going to have to make due with the cow4 alignments I guess!

OH! Good news! Thanks to my checking my old notes, Steve had apparently transferred the old Buffalo data to my QNAP drive! Nice! Let's get the spreadsheets ready on Blade14 instead (where I know my pipeline can run).

> Blade14: /home/dbickhart/vnx_gliu_7/ruminant_project

```bash
# Getting the buffalo spreadsheet ready
ls /mnt/cifs/bickhart-qnap/buffalo_sequence/Sample_*/*.gz | perl -e 'use File::Basename; %h; while(<>){chomp; $b = basename($_); @s = split(/_/, $b); $h{$s[0]}->{$s[4]}->{$s[3]} = $_;} foreach $an (keys(%h)){foreach $num (keys(%{$h{$an}})){ $fq1 = $h{$an}->{$num}->{"R1"}; $fq2 = $h{$an}->{$num}->{"R2"}; print "$fq1\t$fq2\t$an\t$an\n";}}' > buffalo_fq_spreadsheet.tab

# Now the goat
ls /mnt/nfs/nfs2/GoatData/African-Goat-Resequencing/NextSeq500/150820_NS500432_0014_AH15AFBGXX/H15AFBGXX/African-Goat-ReSeq/G0*/*.gz | perl -e 'use File::Basename; %h; while(<>){chomp; $b = basename($_); @s = split(/_/, $b); $h{$s[0]}->{$s[4]}->{$s[3]} = $_;} foreach $an (keys(%h)){foreach $num (keys(%{$h{$an}})){ $fq1 = $h{$an}->{$num}->{"R1"}; $fq2 = $h{$an}->{$num}->{"R2"}; print "$fq1\t$fq2\t$an\t$an\n";}}' > goat_spreadsheet.tab

# Just some minor housecleaning -- I don't want to waste time with EFF and SNP calling!
cat buffalo_fq_spreadsheet.tab goat_spreadsheet.tab > combined_spreadsheet.tab
cp ~/.mergedpipeline.cnfg ./temp_pipeline.cnfg


perl ~/perl_toolchain/sequence_data_pipeline/runMergedBamPipeline.pl --fastqs combined_spreadsheet.tab --output goat_buff_bams --reference ../reference/umd3_kary_unmask_ngap.fa --config temp_pipeline.cnfg --threads 20