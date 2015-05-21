# African cattle sample processing
---

*5/18/2015*

These are my notes on generating statistics and variant calls on the african cattle samples. I will be referring to my [notes from my older note style](https://github.com/njdbickhart/labnotes/blob/master/textpad_archive/Lab_notes_20141022_dbick_running_african_waterbuff.note) when appropriate. Sadly, it looks like I had to take notes on the server and that those notes were lost. 

Still, the SNP calling and other methods are pretty standard, so I should be able to reproduce them easily.

## Getting BAM stats

This should be pretty easy to do. I can pull the X coverage from the BAMs by using my perl script:

> Blade14: /mnt/iscsi/vnx_gliu_7/african_cattle

```bash
# I need to get the names of the bams in comma delimited form
for i in ./*/*merged.sorted.bam; do echo -n "$i,"; done; echo
	./Ndama177/Ndama177.rehead.merged.sorted.bam,./Ndama180/Ndama180.merged.sorted.bam,./Ndama186/Ndama186.merged.sorted.bam,./Ndama188/Ndama188.merged.sorted.bam,./Ndama190/Ndama190.merged.sorted.bam,./Ndama199/Ndama199.merged.sorted.bam,./NS1-01/NS1-01.merged.sorted.bam,./NS1-02/NS1-02.merged.sorted.bam,./NS1-04/NS1-04.merged.sorted.bam,./NS1-05/NS1-05.merged.sorted.bam,./NS1-07/NS1-07.merged.sorted.bam,./NS1-09/NS1-09.merged.sorted.bam 

perl ~/perl_toolchain/sequence_data_scripts/getBamStats.pl -b ./Ndama177/Ndama177.rehead.merged.sorted.bam,./Ndama180/Ndama180.merged.sorted.bam,./Ndama186/Ndama186.merged.sorted.bam,./Ndama188/Ndama188.merged.sorted.bam,./Ndama190/Ndama190.merged.sorted.bam,./Ndama199/Ndama199.merged.sorted.bam,./NS1-01/NS1-01.merged.sorted.bam,./NS1-02/NS1-02.merged.sorted.bam,./NS1-04/NS1-04.merged.sorted.bam,./NS1-05/NS1-05.merged.sorted.bam,./NS1-07/NS1-07.merged.sorted.bam,./NS1-09/NS1-09.merged.sorted.bam -g 2800000000
```

The stats:

BamName |RawXCov| MapXCov| AvgRawChrcov|    AvgMapChrcov
|:--- | :---: | :---: | :---:| :---:|
./Ndama177/Ndama177.rehead.merged.sorted.bam|    5.08554830661171 |       5.04423094543229|        4.9872021856494| 4.95973458460606
./Ndama180/Ndama180.merged.sorted.bam|   3.43362628419845  |      3.33276396957675  |      3.38720573548025 |       3.29383298185112
./Ndama186/Ndama186.merged.sorted.bam |  3.94761317018214   |     3.91129512426424   |     3.90226247462167  |      3.87776681944447
./Ndama188/Ndama188.merged.sorted.bam |  3.54640113469154   |     3.4662634888881 |3.48522278516829       | 3.4145177870431
./Ndama190/Ndama190.merged.sorted.bam |  2.85545393694522   |     2.82615565691045 |       2.73915620565214  |      2.72317887625559
./Ndama199/Ndama199.merged.sorted.bam |  4.34291864241651   |     4.30902938842751 |       4.16805178636062  |      4.15318189746619
./NS1-01/NS1-01.merged.sorted.bam     |  4.29751395859412   |     4.23579733387879 |       4.27254320629197  |      4.22422124673683
./NS1-02/NS1-02.merged.sorted.bam     |  2.12785429219184   |     2.09493544663026 |       2.09302372587782  |      2.0803988709479
./NS1-04/NS1-04.merged.sorted.bam     |  2.55297422157921   |     2.53565774892409 |       2.49533081901342  |      2.48357798953967
./NS1-05/NS1-05.merged.sorted.bam     |  2.55480712670914   |     2.53728295340023 |       2.52324810704291  |      2.51095001113714
./NS1-07/NS1-07.merged.sorted.bam     |  3.17198369816393   |     3.14750402015737 |       3.11683666891432  |      3.0984912034217
./NS1-09/NS1-09.merged.sorted.bam     |  2.69971340498915   |     2.67534076569672 |       2.6780372123013| 2.66157007183874

George just told me about samples that I didn't know we had access to. I'm going to run my pipeline on the samples and then merge the bams individually. The sample fastqs are located here:

> Blade14: /mnt/nfs/nfs2/gliu/HiSeq/150218_SN644_0231_AC5FVPACXX/C5FVPACXX_unaligned/Project_G.Liu.African-Cattle/

I'm going to collect the fastqs and run them with my pipeline script.

> Blade14: /mnt/iscsi/vnx_gliu_7/african_cattle

```bash
ls /mnt/nfs/nfs2/gliu/HiSeq/150218_SN644_0231_AC5FVPACXX/C5FVPACXX_unaligned/Project_G.Liu.African-Cattle/*/*.gz > C5FVPACXX_flowcell_fastqlist.txt
perl ~/bin/convert_steve_fq_to_spreadsheet.pl C5FVPACXX_flowcell_fastqlist.txt
perl -lane 'print "$F[0]\t$F[1]\t$F[-1].lib1\t$F[-1]";' < temp_test_starter.txt > C5FVPACXX_flowcell_spreadsheet_starter.tab

# OK, that accounts for the new format that I need for my pipeline. Time to see if it actually works!
perl ~/perl_toolchain/sequence_data_pipeline/runMergedBamPipeline.pl --fastqs C5FVPACXX_flowcell_spreatarter.tab --output C5FVPACXX --reference ../reference/umd3_kary_unmask_ngap.fa --coords ../reference/samtools_chr_segs.txt --threads 10
```

*5/21/2015*

--

OK, due to an issue with the pipeline, the merger didn't go ahead because I changed the version of samtools while it was running. Instead, let's try to use the log file to replicate the steps for the merger.

> Blade14: /mnt/iscsi/vnx_gliu_7/african_cattle

```bash
# This pulls the appropriate merger lines from the log file and uses them to generate new commands
grep 'sammerge' C5FVPACXX/MergedBamPipeline.1431978899.log | perl -ne '$_ =~ s/.+ \| sammerge \- //g; $_ =~ s/-h 1/-p/; chomp;  print "$_\n"; system($_);'

# That should merge the bams, then index them
```