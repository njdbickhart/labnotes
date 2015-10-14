# JaRMS testing in advance of other projects
---
*10/05/2015*

The goal is to test out JaRMS and see if it can reliably reproduce the analysis performed by CNVnator.

I'm going to generate a "test genome" of only two chromosomes with a handful of pre-determined CNVs. Basically, I'll steal cattle chromosomes 29 and 28, use my RAPTR-SV simulation scripts, and then start processing the data.

> pwd: /home/dbickhart/share/side_projects/jarms_testdata

```bash
# Grep out the chromosomes I need
samtools faidx ../../umd3_data/umd3_kary_unmask_ngap.fa

samtools faidx ../../umd3_data/umd3_kary_unmask_ngap.fa chr28:1-46312546 > chr28.fa
samtools faidx ../../umd3_data/umd3_kary_unmask_ngap.fa chr29:1-51505224 > chr29.fa

# Now I've gotta create a "non-repeats" bed file
echo -e "chr28\t1\t46312546" > chr28.len
echo -e "chr29\t1\t51505224" > chr29.len

subtractBed -a chr28.len -b ../../umd3_data/umd3_repeat_simple_repbase_named_merged.bed > chr28.non_repeats.bed
subtractBed -a chr29.len -b ../../umd3_data/umd3_repeat_simple_repbase_named_merged.bed > chr29.non_repeats.bed

# OK, let's introduce the CNVs!
perl ../../RAPTR-SV/simulations/generateTANDAndDelBedRegions.pl chr28.non_repeats.bed chr28.initialsvs.bed 10
perl ../../RAPTR-SV/simulations/generateTANDAndDelBedRegions.pl chr29.non_repeats.bed chr29.initialsvs.bed 10

perl ../../RAPTR-SV/simulations/introduceTandAndDupsLONG.pl -r chr28.fa -i chr28.initialsvs.bed -o chr28.initialsvs.fa
perl ../../RAPTR-SV/simulations/introduceTandAndDupsLONG.pl -r chr29.fa -i chr29.initialsvs.bed -o chr29.initialsvs.fa

# Now, to combine both original fastas for the reference genome, and then combine both SV fastas for the reads
cat chr28.fa chr29.fa > unaltered_sim_refgenome.fa
cat chr28.initialsvs.fa chr29.initialsvs.fa > altered_sv_template.fa

# Creating the fastq files. Aiming for 10X coverage of the two chromosomes
wgsim -1 100 -2 100 -N 10000000 altered_sv_template.fa altered_sv_read_1.fq altered_sv_read_2.fq > SNP_INDEL_locations.tab

bwa index unaltered_sim_refgenome.fa

bwa mem -R "@RG\tID:test1\tLB:test1\tSM:test1" unaltered_sim_refgenome.fa altered_sv_read_1.fq altered_sv_read_2.fq | samtools view -bS - | samtools sort -T test.sort -o test1.altered_sv.bam -
```

After a long time attempting to install CNVnator and Root on Ubuntu, I got it working. Some notes:

* Must use version 5 of root
* Change MakeFile to add the commands listed [here](http://askubuntu.com/questions/521706/error-adding-symbols-dso-missing-from-command-line-while-compiling-g13-driver)
* Must install the following packages in ubuntu: sudo apt-get install git dpkg-dev make g++ gcc binutils libx11-dev libxpm-dev libxft-dev libxext-dev libusb-1.0.0-dev
* Must set ROOTSYS environmental variable
* Must compile Samtools directory in the CNVnator src directory

OK, after all that, let's create the calls.

```bash
# CNVnator callset
~/CNVnator_v0.3/src/cnvnator -root test.root -genome unaltered_sim_refgenome.fa -tree test1.altered_sv.bam
~/CNVnator_v0.3/src/cnvnator -root test.root -genome unaltered_sim_refgenome.fa -his 100
~/CNVnator_v0.3/src/cnvnator -root test.root  -stat 100
#NOTE: it keeps telling me that "fit is empty" but I suspect it is because it is looking for X and Y
~/CNVnator_v0.3/src/cnvnator -root test.root  -eval 100
~/CNVnator_v0.3/src/cnvnator -root test.root -partition 100
~/CNVnator_v0.3/src/cnvnator -root test.root -call 100 > test1.altered.cnvnator.calls

# That took, literally 10 minutes.
wc -l test1.altered.cnvnator.calls
	580 test1.altered.cnvnator.calls	<- total calls
```

Just for bad luck, let's test JaRMS out.

```bash
# Found out I've gotta index my bam first!
samtools index test1.altered_sv.bam
# TODO: implement this in the program

~/jdk1.8.0_60/bin/java -jar ../../programs_source/netbeans_workspace/JaRMS/dist/JaRMS.jar call -i test1.altered_sv.bam -f unaltered_sim_refgenome.fa -o test1.altered.jarms.calls -t 2
[CALLMODE] Error processing bam file!
Exception in thread "main" java.lang.NullPointerException
        at jarms.modes.CallMode.lambda$run$19(CallMode.java:78)
        at java.util.stream.ReferencePipeline$3$1.accept(ReferencePipeline.java:193)
        at java.util.concurrent.ConcurrentHashMap$KeySpliterator.forEachRemaining(ConcurrentHashMap.java:3527)
        at java.util.stream.AbstractPipeline.copyInto(AbstractPipeline.java:481)
        at java.util.stream.AbstractPipeline.wrapAndCopyInto(AbstractPipeline.java:471)
        at java.util.stream.ReduceOps$ReduceOp.evaluateSequential(ReduceOps.java:708)
        at java.util.stream.AbstractPipeline.evaluate(AbstractPipeline.java:234)
        at java.util.stream.ReferencePipeline.reduce(ReferencePipeline.java:474)
        at jarms.modes.CallMode.run(CallMode.java:79)
        at jarms.JaRMS.main(JaRMS.java:51)

```

I've resolved issues with the program and it produces results, but it lacks duplication calls! 

```bash
# CNVnator
perl ../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f test1.altered.cnvnator.calls -c 0 -m
```

|Entry       | Count|
|:-----------|-----:|
|deletion    |   574|
|duplication |     6|

```bash
# JaRMS prototype
perl ../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f test1.altered.jarms.callscnvs.bed -c 3 -m
```

|Entry    | Count|
|:--------|-----:|
|deletion |   163|

I suspect that there may be some problems with:

* My read depth estimation
* The "levels" calculation in the mean shift implementation
* Some statistical test for average/stdev calculation

I'm going to test out the first and then check each one in turn until I get the right results. First, lets get an unbiased estimate of read depth from BedTools:

```bash
# Making the bed files
perl -e '$maxlen = 46312546; for($i = 0; $i < $maxlen; $i += 500){ $e = $i + 499; if($e > $maxlen){$e = $maxlen;} print "chr28\t$i\t$e\n";}' > chr28.500bp.wins.bed
perl -e '$maxlen = 51505224; for($i = 0; $i < $maxlen; $i += 500){ $e = $i + 499; if($e > $maxlen){$e = $maxlen;} print "chr29\t$i\t$e\n";}' > chr29.500bp.wins.bed

cat chr28.500bp.wins.bed chr29.500bp.wins.bed > unaltered.500bp.wins.bed

coverageBed -a unaltered.500bp.wins.bed -b test1.altered_sv.bam -hist -sorted > unaltered.500bp.wins.rawcov.tab

cat unaltered.500bp.wins.rawcov.tab | cut -f4 | statStd.pl
total   195637
Minimum 0
Maximum 276
Average 122.257318
Median  123
Standard Deviation      17.419157
Mode(Highest Distributed Value) 123
```

Only slightly different from my calculations. Let's try to pull specific regions for analysis, starting with a known duplication.

```bash
perl -lane 'if($F[0] eq "chr28" && $F[2] > 7960546 && $F[1] < 7962887){print $_;}' < unaltered.500bp.wins.rawcov.tab
chr28   7960500 7960999 229     499     499     1.0000000
chr28   7961000 7961499 253     499     499     1.0000000
chr28   7961500 7961999 260     499     499     1.0000000
chr28   7962000 7962499 253     499     499     1.0000000
chr28   7962500 7962999 217     499     499     1.0000000

# OK, good reference points! Double the coverage here.
# Let's get a deletion too

perl -lane 'if($F[0] eq "chr28" && $F[2] > 35296378 && $F[1] < 35297411){print $_;}' < unaltered.500bp.wins.rawcov.tab
chr28   35296000        35296499        87      380     499     0.7615231
chr28   35296500        35296999        0       0       499     0.0000000
chr28   35297000        35297499        27      91      499     0.1823647

# OK, at least 2 500bp windows that are lower than the average in an appreciable way.
```

Now that I have the target windows, let's check them in my IDE.

| chr | start | end | reads | gccorrect | level |
| :--- | :--- | :--- | ---: | ---: | ---: |
chr28|   7960500 | 7960999 | 200 | 199 | 207
chr28 |  7961000 | 7961499 | 208 | 207 | 207
chr28 |  7961500 | 7961999 | 225 | 222 | 207
chr28 |  7962000 | 7962499 | 208 | 207 | 207
chr28 |  7962500 | 7962999 | 294 | 193 | 207
chr28 |  35296000    |    35296499 | 77 | 76 | 99
chr28 |  35296500    |    35296999 | 0  | 0  |  37
chr28 |  35297000    |    35297499 | 23 | 22 | 37

Everything looks fine! Even better, excellent! But the duplication isn't being called somehow. Let's check inside the routine. 

OK, it looks like there's modification of start and end coordinates within a subroutine. Because Java passes by identity, there's an issue here! I changed the Ttest to give the CDF rather than the density value -- let's see if that works instead!

I'm really confused with the "GetMeanSigma" function. I'm not sure how the values were fitted there. I've recompiled CNVnator to get those values.

```bash
~/CNVnator_v0.3/src/cnvnator -root test2.root -genome unaltered_sim_refgenome.fa -tree test1.altered_sv.bam
~/CNVnator_v0.3/src/cnvnator -root test2.root -genome unaltered_sim_refgenome.fa -his 500
~/CNVnator_v0.3/src/cnvnator -root test2.root  -stat 500
getMeanSigma:   mean:   102.231 sigma:  14.9362 After corr:     mean:   103.44  sigma:  11.5275
getMeanSigma, round2:   103.433 11.313  minfit: 80.3846 maxfit: 126.495
Average RD per bin (1-22) is 103.433 +- 11.313 (before GC correction)
...
Correcting counts by GC-content for 'chr28' ...
getMeanSigma:   mean:   102.231 sigma:  14.9362 After corr:     mean:   103.44  sigma:  11.5275
getMeanSigma, round2:   103.433 11.313  minfit: 80.3846 maxfit: 126.495
Zero value of GC average.
Bin 7146 with center 3.57275e+06 is not corrected.
Correcting counts by GC-content for 'chr29' ...
getMeanSigma:   mean:   102.231 sigma:  14.9362 After corr:     mean:   103.44  sigma:  11.5275
getMeanSigma, round2:   103.433 11.313  minfit: 80.3846 maxfit: 126.495
Zero value of GC average.
Bin 50448 with center 2.52238e+07 is not corrected.
Zero value of GC average.
Bin 92731 with center 4.63652e+07 is not corrected.
Zero value of GC average.
Bin 92903 with center 4.64512e+07 is not corrected.
Zero value of GC average.
Bin 95317 with center 4.76582e+07 is not corrected.
Making statistics for chr28 after GC correction ...
Making statistics for chr29 after GC correction ...
getMeanSigma:   mean:   103.163 sigma:  15.098  After corr:     mean:   104.443 sigma:  11.6203
getMeanSigma, round2:   104.424 11.3559 minfit: 81.2025 maxfit: 127.684
Average RD per bin (1-22) is 104.424 +- 11.3559 (after GC correction)
...


~/CNVnator_v0.3/src/cnvnator -root test2.root  -eval 500
~/CNVnator_v0.3/src/cnvnator -root test2.root  -partition 500
~/CNVnator_v0.3/src/cnvnator -root test2.root -call 500 > test1.altered.cnvnator.500.calls

# I removed the sex chromosome fit error listings
```


```c++
// Root commands
TFile file("test2.root")
file.ls()
TDirectory *d = (TDirectory*)file.Get("bin_500")
d->cd()
d->ls()
TH1 *his = (TH1*)d->Get("rd_p_GC_500")
TF1 *fg = new TF1("fb", "gaus(0)", 0, 10)
double mean = his->GetMean()
mean
(double)1.03162798582732051e+02
double sigma = his->GetRMS()
sigma
(double)1.50979568939763453e+01
double constant = his->GetEntries() * 0.4/sigma
constant
(double)5.18313839081243077e+03

fg->SetParameters(constant, mean, 0.5*sigma)
fg->SetParLimits(1, 0, 3 * mean)
fg->SetParLimits(2, 0, 3 * sigma)
his->Fit(fg, "qN")
 his->Draw()
double min_fit = fg->GetParameter(1) - 2.0*fg->GetParameter(2)
double max_fit = fg->GetParameter(1) + 2.0*fg->GetParameter(2)
his->Fit(fg, "qN", "", min_fit, max_fit)
```

OK! I figured it out! I used the Apache commons Gaussian fitter function after "binning" the signal into rounded integer read count bins. Now I'm getting similar standard deviation values in the mean shift algorithm and during SV segmentation. This is also improving my e value estimates and is catching more true positives. 

Let's test how many we've gotten over the base CNVnator algorithm.

```bash
grep 'chr28' test1.altered.jarms.callscnvs.bed | cut -f1,2,3,4 | intersectBed -a chr28.initialsvs.bed -b stdin
chr28   7962001 7962887 TAND
chr28   13321161        13322000        TAND
chr28   38523501        38524500        TAND

perl -lane '($c, $s, $e) = $F[1] =~ /(chr.+):(\d+)-(\d+)/; if($c eq "chr28"){print "$c\t$s\t$e\t$F[0]";}' < test1.altered.cnvnator.500.calls | intersectBed -a chr28.initialsvs.bed -b stdin
chr28   35296501        35297411        DEL
chr28   44928150        44929771        DEL

# Interestingly, if I use the 100 bin wins, I get this:
perl -lane '($c, $s, $e) = $F[1] =~ /(chr.+):(\d+)-(\d+)/; if($c eq "chr28"){print "$c\t$s\t$e\t$F[0]";}' < test1.altered.cnvnator.calls | intersectBed -a chr28.initialsvs.bed -b stdin
chr28   7960546 7962887 TAND
chr28   17103516        17104500        TAND
chr28   35296401        35297400        DEL
chr28   38522662        38524563        TAND
chr28   44928201        44929771        DEL

grep 'chr29' test1.altered.jarms.callscnvs.bed | cut -f1,2,3,4 | intersectBed -a chr29.initialsvs.bed -b stdin
chr29   1749117 1750206 DEL
chr29   20259148        20262500        DEL
chr29   21823225        21826779        DEL
chr29   22600501        22601457        TAND
chr29   23690522        23697000        DEL
chr29   40319001        40322301        DEL

perl -lane '($c, $s, $e) = $F[1] =~ /(chr.+):(\d+)-(\d+)/; if($c eq "chr29"){print "$c\t$s\t$e\t$F[0]";}' < test1.altered.cnvnator.500.calls | intersectBed -a chr29.initialsvs.bed -b stdin
chr29   1749117 1750206 DEL
chr29   20259148        20262500        DEL
chr29   21823225        21826779        DEL
chr29   23690522        23697000        DEL
chr29   40318824        40322301        DEL

# And the 100 bp bins
perl -lane '($c, $s, $e) = $F[1] =~ /(chr.+):(\d+)-(\d+)/; if($c eq "chr29"){print "$c\t$s\t$e\t$F[0]";}' < test1.altered.cnvnator.calls | intersectBed -a chr29.initialsvs.bed -b stdin
chr29   1749117 1750206 DEL
chr29   17213201        17213600        DEL
chr29   19499451        19500286        DEL
chr29   19623001        19624000        DEL
chr29   20259201        20262500        DEL
chr29   21823225        21826779        DEL
chr29   22599062        22601400        TAND
chr29   23690522        23697180        DEL
chr29   34373733        34375000        TAND
chr29   36237119        36238580        TAND
chr29   40318824        40322300        DEL

wc -l test1.altered.cnvnator.500.calls test1.altered.cnvnator.calls test1.altered.jarms.callscnvs.bed
  426 test1.altered.cnvnator.500.calls	<- this had 8 lines dedicated to my rewrite of the source code.
  580 test1.altered.cnvnator.calls
   53 test1.altered.jarms.callscnvs.bed
 1059 total
```

So, chr29 was the best for my implementation, and I had weird results with the chr28 comparisons. I think that my gc correction may be causing the difference here.

Some c modifications to debug further:

```cpp
// testtwo
if(ret < 0.05){
        cerr<<"ttwo:\t"<<ret<<"\t"<<preret<<"\t"<<ndf<<"\t"<<scale<<"\t";
        cerr<<m1<<"\t"<<s1<<"\t"<<n1<<"\t"<<m2<<"\t"<<s2<<"\t"<<n2<<endl;
  }

//evalue
if (ret < 0.05){
        cerr<<"Eval:\t"<<ret<<"\t"<<start<<"\t"<<end<<"\t"<<mean<<"\t"<<sigma<<"\t"<<aver;
        for(int i = start; i <= end; i++){
                cerr<<"\t"<<rd[i];
        }
        cerr<<endl;
  }
```

I actually found out that I had not set the "potential deletion" region correctly! Oops! Fixed it and the data looks better.

