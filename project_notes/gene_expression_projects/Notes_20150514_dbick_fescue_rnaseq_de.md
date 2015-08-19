# Processing Fescue RNA-seq results
---
*5/14/2015*

## Table of Contents

* [Initial Information](#initial)
* [Identifying differentially expressed transcripts](#rsem)
* [Generating the results](#results)
* [Running the Analysis](#running)
* [Interpreting the results](#interpret)
* [Group Specific Results](#specific)
* [Dissecting more useful information](#dissect)

<a name="initial"></a>
## Initial information on files

The files are all located at this address:

> Blade14: /mnt/cifs/bickhart-qnap/RBaldwin/EA12095_20130418/

And here are my previous notes on the subject. 

I received an email from Randy regarding the fields of the file that I had discovered previously. Here is the file:

> FILE:  /mnt/cifs/bickhart-qnap/RBaldwin/EA12095_20130418/20130430-H22821.Analysis/3013_1_26_11_L8.LB3.rsem.genes.results.gz 
CONTENTS:
1433E_BOVIN:ENSBTAG00000005664      10301.02              0.000405000950203639  NM_174491:ENSBTAT00000007442
1433S_BOVIN:ENSBTAG00000009223      725.00   3.75504296409274e-05   NM_001075912:ENSBTAT00000012154
1433T_BOVIN:ENSBTAG00000002108      7197.26 0.00026155699827383     NM_001078127:ENSBTAT00000032851

And here is the email response from the technician:

> Hi Randy,

>This is what I found out:

>The genes.results file is formatted like this:
Gene ID
Estimated count (NOT normalized)
Tau (a type of normalization performed by RSEM)
Transcript ID(s) associated with the gene

>Let me know if you need anything else.

>Amy

So, RSEM is the magic word here. I found the RSEM package online at this [address](http://deweylab.biostat.wisc.edu/rsem/README.html). I suspect that I can install the software and use it to get a better sense of the current analysis and to generate figures/tables for Randy to use. 

It looks like RSEM is a combination of perl scripts, C++ hardpoints and R packages. Let's get the [github repository](https://github.com/deweylab/RSEM) for RSEM installed so i can start tooling around with the packages.

> Blade14: ~/

```bash
git clone https://github.com/deweylab/RSEM

cd RSEM
make
make ebseq

# NOTE: ebseq installed several R packages in my local install of R
export PATH=$PATH:/home/dbickhart/RSEM
```
Here is some more documentation of [EBseq](https://www.biostat.wisc.edu/~kendzior/EBSEQ/) and the [EBseq-RSEM pipeline](http://deweylab.biostat.wisc.edu/rsem/rsem-run-ebseq.html). 

It looks like I have to run [rsem-run-ebseq](http://deweylab.biostat.wisc.edu/rsem/rsem-run-ebseq.html) then run [rsem-control-fdr](http://deweylab.biostat.wisc.edu/rsem/rsem-control-fdr.html) in order to get good differentially expressed transcripts. My only problem: I don't know how to designate samples vs replicates in the rsem-run-ebseq program. I guess I just need to start running trials of the program in order to intimate this.

<a name="rsem"></a>
## Identifying differentially expressed transcripts

OK, now that I've install the software, I think that I understand which files need to be loaded in order to identify differentially expressed transcripts. I've gotta just bite the bullet and try to run the program and get the error messages. Hopefully the sample context discovery is sensitive to the sample names. If it is sensitive to the order in the matrix file, I can reformat it using Perl.

> /mnt/iscsi/vnx_gliu_7/rna_seq/RBaldwin

```bash
# First, I need to find out how many samples and conditions there are
head -n 1 /mnt/cifs/bickhart-qnap/RBaldwin/EA12095_20130418/Batch.2013-04-18/rsem_gnorm_matrix.txt | perl -lane 'print scalar(@F);'
	39

# OK, now to divide them up. Let's do this by eye by looking at the organization of the sample names in the matrix header
head -n 1 /mnt/cifs/bickhart-qnap/RBaldwin/EA12095_20130418/Batch.2013-04-18/rsem_gnorm_matrix.txt
        H22821/2686_10_13_10_L7.LB27    H22821/2686_11_17_10_L8.LB1     H22821/2819_10_20_10_L7.LB8
     H22821/2866_1_21_11_L8.LB9     H22821/3006_12_29_10_L6.LB4     H22821/3006_1_31_11_L6.LB5  
    H22821/3008_11_17_10_L4.LB2     H22821/3013_1_26_11_L8.LB3     H22821/3025_2_9_11_L4.LB21
      H22821/3025_4_18_11_L6.LB22     H22821/3032_5_4_11_L4.LB6       H22821/3032_6_8_11_L6.LB7
      H22821/3032_7_5_11_L4.LB8       H22821/3036_5_25_11_L7.LB5      H22821/3041_12_30_10_L5.LB3
     H22821/3065_11_24_10_L4.LB18   H22821/3065_12_29_10_L5.LB19    H22821/3073_12_8_10_L7.LB6
      H22821/3073_2_14_11_L7.LB7      H22821/3076_1_12_11_L8.LB14    H22821/3076_2_18_11_L7.LB15
     H22821/3079_2_23_11_L8.LB10     H22821/3079_4_4_11_L8.LB11      H22821/3085_12_29_10_L5.LB23
   H22821/3087_12_8_10_L5.LB9      H22821/3087_1_12_11_L6.LB10     H22821/3087_2_16_11_L5.LB11
     H22821/3088_1_5_11_L5.LB12     H22821/3088_2_9_11_L4.LB13      H22821/3088_3_11_11_L6.LB14
     H22821/3094_2_16_11_L5.LB15     H22821/3094_4_22_11_L4.LB16    H22821/3097_3_16_11_L6.LB20
     W22942/3013_12_22_10_L2.LB2     W22942/3036_4_20_11_L2.LB4      W22942/3077_3_7_11_L2.LB16
     W22942/3078_12_22_10_L2.LB13    W22942/3083_6_13_11_L2.LB12     W22942/3085_2_2_11_L2.LB25
```

*5/15/2015*

--

Actually, I'm curious to see how the program will run without this division. 

> Blade14: /mnt/iscsi/vnx_gliu_7/rna_seq/RBaldwin/

```bash
rsem-run-ebseq /mnt/cifs/bickhart-qnap/RBaldwin/EA12095_20130418/Batch.2013-04-18/rsem_gnorm_matrix.txt 3,3,3,3,3,3,3,3,3,3,3,3,3 first_try_disordered
	Loading required package: blockmodeling
	Loading required package: gplots

	Attaching package: âgplotsâ

	The following object is masked from âpackage:statsâ:

	    lowess

	Removing transcripts with 75 th quantile < = 10
	14833transcripts will be tested"rsem-for-ebseq-find-DE /home/dbickhart/RSEM/EBSeq # /mnt/cifs/bickhart-qnap/RBaldwin/EA12095_20130418/Batch.2013-04-18/rsem_gnorm_matrix.txt first_try_disordered 3 3 3 3 3 3 3 3 3 3 3 3 3" failed! Plase check if you provide correct parameters/options for the pipeline!

```

OK that didn't work, time to contact the author. 

*5/18/2015*

--

While I do that, let's try a different setting that does not have as many divisions:

> Blade14:  /mnt/iscsi/vnx_gliu_7/rna_seq/RBaldwin/

```bash
rsem-run-ebseq /mnt/cifs/bickhart-qnap/RBaldwin/EA12095_20130418/Batch.2013-04-18/rsem_gnorm_matrix.txt 13,13,13 second_try_disordered
	rsem-for-ebseq-find-DE /home/dbickhart/RSEM/EBSeq # /mnt/cifs/bickhart-qnap/RBaldwin/EA12095_20130418/Batch.2013-04-18/rsem_gnorm_matrix.txt second_try_disordered 13 13 13
	Loading required package: blockmodeling
	Loading required package: gplots

	Attaching package: âgplotsâ

	The following object is masked from âpackage:statsâ:

    	lowess

	Removing transcripts with 75 th quantile < = 10
	14833transcripts will be testediteration 1 done

	time 712.94

	time 712.94

	iteration 2 done

	time 713.7

	iteration 3 done

	time 710.84

	iteration 4 done

	time 707.56

	iteration 5 done

	time 286.05

	Error in file(file, ifelse(append, "a", "w")) :
  	cannot open the connection
	Calls: write.table -> file
	In addition: Warning message:
	In file(file, ifelse(append, "a", "w")) :
  	cannot open file 'norm_/mnt/cifs/bickhart-qnap/RBaldwin/EA12095_20130418/Batch.2013-04-18/rsem_gnorm_matrix.txt': No such file or directory
	Execution halted
	"rsem-for-ebseq-find-DE /home/dbickhart/RSEM/EBSeq # /mnt/cifs/bickhart-qnap/RBaldwin/EA12095_20130418/Batch.2013-04-18/rsem_gnorm_matrix.txt second_try_disordered 13 13 13" failed! Plase check if you provide correct parameters/options for the pipeline!

```

Damn! Something in the script screwed up the write.table call. Let's dig into the guts of RSEM to try to run it from R instead. AH, ok, I see why, from the Rscript: [**rsem-for-ebseq-find-DE**](https://github.com/deweylab/RSEM/blob/master/EBSeq/rsem-for-ebseq-find-DE) the authors made a silly mistake here:

```R
path <- argv[1]
ngvector_file <- argv[2]
data_matrix_file <- argv[3]
output_file <- argv[4]
norm_out_file <- paste0("norm_", data_matrix_file)
# Above is the mistake ... all file paths are prepended with "norm_" regardless of file path
```

Let us dissect how the program handles the replicates (are they in the order of the columns of the data matrix?) in the meantime. Again, from the source code on **rsem-for-ebseq-find-DE**.

```R
conditions <- as.factor(rep(paste("C", 1:nc, sep=""), times = num_reps))
# OK, so this is according to the documentation -- each comma delimited 
Sizes <- MedianNorm(DataMat)
NormMat <- GetNormalizedMat(DataMat, Sizes)
ngvector <- NULL
if (ngvector_file != "#") {
  ngvector <- as.vector(data.matrix(read.table(ngvector_file)))
  stopifnot(!is.null(ngvector))
}

if (nc == 2) {
  	EBOut <- NULL
  	EBOut <- EBTest(Data = DataMat, NgVector = ngvector, Conditions = conditions, sizeFactors = Sizes, maxround = 5)
 	# OK, EBTest is part of the EBSeq package. I'll refer to that documentation for the remainder of the test
	 
	stopifnot(!is.null(EBOut))

	#...
}
```

Here is the [EBseq](http://www.bioconductor.org/packages/devel/bioc/vignettes/EBSeq/inst/doc/EBSeq_Vignette.pdf) documentation. It looks like my guess was correct -- basically the conditions must be ordered with respect to their column position in the matrix. So I would need to reformat the matrix so that each column is in order with the condition that the sample represents. 


OK, so I also need to copy the files to my working directory. Such poor design. Trying it again.

```bash
cp /mnt/cifs/bickhart-qnap/RBaldwin/EA12095_20130418/Batch.2013-04-18/rsem_gnorm_matrix.txt ./

# Running on the local file
rsem-run-ebseq rsem_gnorm_matrix.txt 13,13,13 second_try_disordered
```

I received three files. Here's a brief description of their contents:

* **second_try_disordered.pattern** : A condition similarity matrix that explains the pattern of significant differential expression
* **second_try_disordered** : the posteriori probabilities that the conditional gene expression profiles fit the patterns
* **second_try_disordered.condmeans** : the median expression across all conditions for the samples

*05/19/2015*

-- 
<a name="results"></a>
## Generating the results

OK, my strategy is to divide up the data into severate categories based on the information that Tony Capuco gave me. There are 9 different conditions, but rather than confuse the pipeline, I prefer to separate the gene expression matrix into three "bitesize" chunks. Here are the conditions:

* B
* C
* I

And it appears that there were three different biopsy tissue types:

* D
* L1
* L2

Rather than ask for the explanation of these tissues, I'm going to treat this like a "blind" study and treat the samples independently. One slight cause for concern is that cow 3094 in the **C - L2** group had a [RNA Integrity Number (RIN)](http://www.genomics.agilent.com/article.jsp?pageId=2181&_requestid=383577) lower than 4.8. I will check the pipeline's "condmeans" to see if the control group misbehaves significantly if I remove 3094 entirely.

##### Data manipulation

I'm going to subsection the matrix using Perl, just because I really don't feel like using R. Let's generate the sample table from the rsem_gnorm_matrix.txt:

> Blade14: /mnt/iscsi/vnx_gliu_7/rna_seq/RBaldwin

```bash
head -n 1 rsem_gnorm_matrix.txt | perl -lane 'print join("\n", @F);'
```

I'm taking care of the association in the excel file that Tony provided. I then cut and pasted the text over to a new text file in vim:

```bash
vim sample_key_file.tab
dos2unix sample_key_file.tab
# Just to print a table for my records
perl -e 'print "Cow | Cond | Biopsy | File\n"; print ":--- | :--- | :--- | :---\n"; while(<>){chomp; @s = split(/\t/); ($c) = $s[2] =~ m/.+\/(\d{4})_.+/; print "$c | $s[0] | $s[1] | $s[2]\n";}' < sample_key_file.tab
```

Cow | Cond | Biopsy | File
:--- | :--- | :--- | :---
2686 | I | L1 | H22821/2686_10_13_10_L7.LB27
2686 | I | D | H22821/2686_11_17_10_L8.LB1
2819 | B | L1 | H22821/2819_10_20_10_L7.LB8
2866 | B | L2 | H22821/2866_1_21_11_L8.LB9
3006 | B | L2 | H22821/3006_12_29_10_L6.LB4
3006 | B | D | H22821/3006_1_31_11_L6.LB5
3008 | C | D | H22821/3008_11_17_10_L4.LB2
3013 | B | D | H22821/3013_1_26_11_L8.LB3
3013 | B | L1 | W22942/3013_12_22_10_L2.LB2
3025 | I | L1 | H22821/3025_2_9_11_L4.LB21
3025 | I | L2 | H22821/3025_4_18_11_L6.LB22
3032 | C | L1 | H22821/3032_5_4_11_L4.LB6
3032 | C | D | H22821/3032_6_8_11_L6.LB7
3032 | C | L2 | H22821/3032_7_5_11_L4.LB8
3036 | B | L1 | H22821/3036_5_25_11_L7.LB5
3036 | B | D | W22942/3036_4_20_11_L2.LB4
3041 | I | L2 | H22821/3041_12_30_10_L5.LB3
3065 | C | L1 | H22821/3065_11_24_10_L4.LB18
3065 | C | D | H22821/3065_12_29_10_L5.LB19
3073 | B | L1 | H22821/3073_12_8_10_L7.LB6
3073 | B | L2 | H22821/3073_2_14_11_L7.LB7
3076 | I | D | H22821/3076_1_12_11_L8.LB14
3076 | I | L2 | H22821/3076_2_18_11_L7.LB15
3077 | I | L2 | W22942/3077_3_7_11_L2.LB16
3078 | I | D | W22942/3078_12_22_10_L2.LB13
3079 | B | D | H22821/3079_2_23_11_L8.LB10
3079 | B | L2 | H22821/3079_4_4_11_L8.LB11
3083 | C | L2 | W22942/3083_6_13_11_L2.LB12
3085 | I | L1 | H22821/3085_12_29_10_L5.LB23
3085 | I | D | W22942/3085_2_2_11_L2.LB25
3087 | C | D | H22821/3087_12_8_10_L5.LB9
3087 | C | L1 | H22821/3087_1_12_11_L6.LB10
3087 | C | L2 | H22821/3087_2_16_11_L5.LB11
3088 | C | L1 | H22821/3088_1_5_11_L5.LB12
3088 | C | D | H22821/3088_2_9_11_L4.LB13
3088 | C | L2 | H22821/3088_3_11_11_L6.LB14
3094 | C | L1 | H22821/3094_2_16_11_L5.LB15
3094 | C | L2 | H22821/3094_4_22_11_L4.LB16
3097 | I | L1 | H22821/3097_3_16_11_L6.LB20

OK, the next task is to create condition-biopsy specific gene matrices for the differential expression analysis. I can do this by writing a script to parse the key file and selecting only specific columns to write out. I then need to reorder the samples based on the following scheme: **C -> B -> I**.

I wrote a one-off script to process the matrix file: [separateGeneTableMatrixForRSEM.pl](https://github.com/njdbickhart/perl_toolchain/blob/master/bed_cnv_fig_table_pipeline/separateGeneTableMatrixForRSEM.pl)

> Blade14: /mnt/iscsi/vnx_gliu_7/rna_seq/RBaldwin

```bash
# Going to start with the "D" biopsy condition
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/separateGeneTableMatrixForRSEM.pl -f rsem_gnorm_matrix.txt -t C,B,I -s D -k sample_key_file.tab -o rsem_D_biopsy_ordered_matrix.txt

# Now let's check the fidelity of the file -- I don't want to find out later that it's mangled!
perl -lane 'print(scalar(@F));' < rsem_D_biopsy_ordered_matrix.txt | head -n 5
	13	<- I expected this, because the first entry is a blank line
	14
	14
	14
	14
# I'm still checking it visually
# OK, at first glance it looks good. Let's check the column number fidelity for the data rows

cat rsem_gnorm_matrix.txt | cut -f3 | head -n 5
	H22821/2686_11_17_10_L8.LB1
	50346.5691
	275.2660
	2.6596
	2558.4441

cat rsem_D_biopsy_ordered_matrix.txt | cut -f11 | head -n 5
	H22821/2686_11_17_10_L8.LB1
	50346.5691
	275.2660
	2.6596
	2558.4441

# Everything matches up -- looks good
# Let's replicate it for the remaining types
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/separateGeneTableMatrixForRSEM.pl -f rsem_gnorm_matrix.txt -t C,B,I -s L1 -k sample_key_file.tab -o rsem_L1_biopsy_ordered_matrix.txt

perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/separateGeneTableMatrixForRSEM.pl -f rsem_gnorm_matrix.txt -t C,B,I -s L2 -k sample_key_file.tab -o rsem_L2_biopsy_ordered_matrix.txt

# Just to be tidy, I'm going to move the files to new directories for processing
mkdir biopsy_D
mkdir biopsy_L1
mkdir biopsy_L2
```
<a name="running"></a>
## Running the analysis

The RSEM pipeline uses 5 iterations of EBSeq, and the authors of EBSeq say that this might not be enough. I'm going to run it with the straight pipeline anyways, as the RSEM pipeline also generates a normalized data matrix that can be plugged into R later for plotting/analysis.

> /mnt/iscsi/vnx_gliu_7/rna_seq/RBaldwin/biopsy_D

```bash
# Determining condition/replicate counts
head -n 1 rsem_D_biopsy_ordered_matrix.txt
	# 5,4,4
export PATH=/home/dbickhart/RSEM/:$PATH
rsem-run-ebseq rsem_D_biopsy_ordered_matrix.txt 5,4,4 rsem_D_biopsy.results
```

> /mnt/iscsi/vnx_gliu_7/rna_seq/RBaldwin/biopsy_L1

```bash
head -n 1 rsem_L1_biopsy_ordered_matrix.txt
	# 5,4,4
export PATH=/home/dbickhart/RSEM/:$PATH
rsem-run-ebseq rsem_L1_biopsy_ordered_matrix.txt 5,4,4 rsem_L1_biopsy.results
```

> /mnt/iscsi/vnx_gliu_7/rna_seq/RBaldwin/biopsy_L2

```bash
head -n 1 rsem_L2_biopsy_ordered_matrix.txt
	# 5,4,4
export PATH=/home/dbickhart/RSEM/:$PATH
rsem-run-ebseq rsem_L2_biopsy_ordered_matrix.txt 5,4,4 rsem_L2_biopsy.results
```
<a name="interpret"></a>
## Interpretting the results

OK, the program generates a posterior probability of the likelihood of differential expression (**PPDE**) in the flat "results" file. I'm going conservative and select a **PPDE** of > 95%, which would give me a false discovery rate of ~ 5% (not counting multiple testing). Let's check how many there are in each file.

> Blade 14: /mnt/iscsi/vnx_gliu_7/rna_seq/RBaldwin

```bash
for i in ./*/*.results; do echo -ne "$i\t"; perl -lane 'if($F[-1] > 0.95 && $F[-2] ne "\"Pattern1\""){print $_;}' < $i | wc -l; done
./biopsy_D/rsem_D_biopsy.results        231
./biopsy_L1/rsem_L1_biopsy.results      288
./biopsy_L2/rsem_L2_biopsy.results      284

# Interesting. So let's get the overlap
for i in ./*/*.results; do echo -ne "$i\t"; perl -lane 'if($F[-1] > 0.95 && $F[-2] ne "\"Pattern1\""){print $F[0];}' < $i; done | perl -e '%h; while(<>){chomp; $h{$_} += 1;} %v; foreach my $k (sort{$h{$a} <=> $h{$b}} keys(%h)){print "$k\t$h{$k}\n";}'  | wc -l
	672

# Not too much, sadly! Let's get a count of the number of multiple entries
for i in ./*/*.results; do echo -ne "$i\t"; perl -lane 'if($F[-1] > 0.95 && $F[-2] ne "\"Pattern1\""){print $F[0];}' < $i; done | perl -e '%h; while(<>){chomp; $h{$_} += 1;} %v; foreach my $k (sort{$h{$a} <=> $h{$b}} keys(%h)){$v{$h{$k}} += 1;} foreach my $b (sort{$a <=> $b}keys(%v)){print "$b\t$v{$b}\n";}'
	1       556
	2       101
	3       15

# So, about 116 are shared by two or more biopsies
# Let's print out the list of shared genes and run a David analysis on them.
# I will then do individual biopsy datasets in David to see if any trends emerge
mkdir analysis
for i in ./*/*.results; do echo -ne "$i\t"; perl -lane 'if($F[-1] > 0.95 && $F[-2] ne "\"Pattern1\""){print $F[0];}' < $i; done | perl -e '%h; while(<>){chomp; $h{$_} += 1;} %v; foreach my $k (sort{$h{$a} <=> $h{$b}} keys(%h)){if($h{$k} > 1){ ($g) = $k =~ m/\".*(ENSBTAG.+)\"/; print "$g\n";}}' > analysis/ensgene_two_or_more_biopsies.txt
```

The clustering analysis had an enrichment score of 6.35 for Milk protein functional classes and an enrichment score of 5.16 for secreted proteins/peptides. I believe Randy mentioned something about milk yield, so this is probably a very good (and expected) result! Benjamini p values for milk and milk protein categories was 5.9 x 10^-8. 

I also saved the functional annotation table for parsing. Since most of the information is tab-delimited, I can just send the file as part of an excel spreadsheet.

Let's prepare the ENSGene lists for the conditions, but we'll do it in a pattern specific manner.

```bash
# Pattern specific command
perl -lane 'if($F[-1] > 0.95 && $F[-2] ne "\"Pattern1\""){$F[-2] =~ s/\"//g; open(OUT, ">> analysis/rsem_D_biopsy.$F[-2].ensgene"); ($g) = $F[0] =~ m/\".*(ENSBTAG.+)\"/; print OUT "$g"; close OUT;}' < biopsy_D/rsem_D_biopsy.results
# Everything above 95% prior
perl -lane 'if($F[-1] > 0.95 && $F[-2] ne "\"Pattern1\""){ ($g) = $F[0] =~ m/\".*(ENSBTAG.+)\"/; print "$g";}' < biopsy_D/rsem_D_biopsy.results > analysis/rsem_D_biopsy.95.ensgene

# Now to do the other biopsies 
perl -lane 'if($F[-1] > 0.95 && $F[-2] ne "\"Pattern1\""){$F[-2] =~ s/\"//g; open(OUT, ">> analysis/rsem_L1_biopsy.$F[-2].ensgene"); ($g) = $F[0] =~ m/\".*(ENSBTAG.+)\"/; print OUT "$g"; close OUT;}' < biopsy_L1/rsem_L1_biopsy.results
perl -lane 'if($F[-1] > 0.95 && $F[-2] ne "\"Pattern1\""){$F[-2] =~ s/\"//g; open(OUT, ">> analysis/rsem_L2_biopsy.$F[-2].ensgene"); ($g) = $F[0] =~ m/\".*(ENSBTAG.+)\"/; print OUT "$g"; close OUT;}' < biopsy_L2/rsem_L2_biopsy.results

perl -lane 'if($F[-1] > 0.95 && $F[-2] ne "\"Pattern1\""){ ($g) = $F[0] =~ m/\".*(ENSBTAG.+)\"/; print "$g";}' < biopsy_L1/rsem_L1_biopsy.results > analysis/rsem_L1_biopsy.95.ensgene
perl -lane 'if($F[-1] > 0.95 && $F[-2] ne "\"Pattern1\""){ ($g) = $F[0] =~ m/\".*(ENSBTAG.+)\"/; print "$g";}' < biopsy_L2/rsem_L2_biopsy.results > analysis/rsem_L2_biopsy.95.ensgene

wc -l analysis/*
  116 analysis/ensgene_two_or_more_biopsies.txt
  231 analysis/rsem_D_biopsy.95.ensgene
  177 analysis/rsem_D_biopsy.Pattern2.ensgene
    4 analysis/rsem_D_biopsy.Pattern3.ensgene
   50 analysis/rsem_D_biopsy.Pattern4.ensgene
  288 analysis/rsem_L1_biopsy.95.ensgene
   18 analysis/rsem_L1_biopsy.Pattern2.ensgene
    6 analysis/rsem_L1_biopsy.Pattern3.ensgene
  261 analysis/rsem_L1_biopsy.Pattern4.ensgene
    3 analysis/rsem_L1_biopsy.Pattern5.ensgene
  284 analysis/rsem_L2_biopsy.95.ensgene
   18 analysis/rsem_L2_biopsy.Pattern2.ensgene
  113 analysis/rsem_L2_biopsy.Pattern3.ensgene
  135 analysis/rsem_L2_biopsy.Pattern4.ensgene
   18 analysis/rsem_L2_biopsy.Pattern5.ensgene
```

Here is the general "pattern" notation for reference:

Name | **C** |    **B**  |   **I**
:--- | :--- | :---| :---
"Pattern1"  |    1  |     1  |     1
"Pattern2"   |   1   |    1   |    2
"Pattern3"   |   1   |    2   |    1
"Pattern4"   |   1   |    2   |    2
"Pattern5"   |   1   |    2   |    3

So I need to talk to Randy to get a better sense of the type of significant expression that he's looking for. Let's send the excel spreadsheet and get his input. First, let's generate some quick QQ plots using the [EBSeq manual](http://www.bioconductor.org/packages/devel/bioc/vignettes/EBSeq/inst/doc/EBSeq_Vignette.pdf) and the source code for [rsem-for-ebseq-find-DE](https://github.com/deweylab/RSEM/blob/master/EBSeq/rsem-for-ebseq-find-DE). 

Since this takes forever, let's just work in the "D" directory for now.

> Blade14: /mnt/iscsi/vnx_gliu_7/rna_seq/RBaldwin/biopsy_D

```R
path <- "/home/dbickhart/RSEM/EBSeq"
ngvector_file <- "#"
data_matrix_file <- "rsem_D_biopsy_ordered_matrix.txt"
output_file <- "temp_test"
norm_out_file <- paste0("norm_", data_matrix_file)
nc <- 3
num_reps <- as.numeric(c(5, 4, 4))

# Now, this is just a paste of the pipeline
.libPaths(c(path, .libPaths()))
library(EBSeq)

DataMat <- data.matrix(read.table(data_matrix_file))
n <- dim(DataMat)[2]
if (sum(num_reps) != n) stop("Total number of replicates given does not match the number of columns from the data matrix!")

conditions <- as.factor(rep(paste("C", 1:nc, sep=""), times = num_reps))
Sizes <- MedianNorm(DataMat)
NormMat <- GetNormalizedMat(DataMat, Sizes)
ngvector <- NULL
if (ngvector_file != "#") {
  ngvector <- as.vector(data.matrix(read.table(ngvector_file)))
  stopifnot(!is.null(ngvector))
}

patterns <- GetPatterns(conditions)
eename <- rownames(patterns)[which(rowSums(patterns) == nc)]
stopifnot(length(eename) == 1)

MultiOut <- NULL
MultiOut <- EBMultiTest(Data = DataMat, NgVector = ngvector, Conditions = conditions, AllParti = patterns, sizeFactors = Sizes, maxround = 5)
stopifnot(!is.null(MultiOut))

MultiPP <- GetMultiPP(MultiOut)

PP <- as.data.frame(MultiPP$PP)
pos <- which(names(PP) == eename)
probs <- rowSums(PP[,-pos])

results <- cbind(PP, MultiPP$MAP[rownames(PP)], probs)
colnames(results) <- c(colnames(PP), "MAP", "PPDE")  
ord <- order(results[,"PPDE"], decreasing = TRUE)
results <- results[ord,]

pdf(file="biopsy_D_qqplot.pdf", useDingbats=FALSE)
par(mfrow=c(2,2))
QQP(MultiOut)
dev.off()
# I saved the workspace image -- just incase!
```

The QQ plots don't look too bad. A little right-skew for **C** and **B**, but other than that, not too shabby.

*5/20/2015*

--

OK, Let's start segregating the datat to look for group-specific gene expression results.

<a name="specific"></a>
## Group specific results

The Goal here is to treat the **D** group differently from the **L1** and **L2** groups. So I will attempt to identify DE genes that are common to the **L** types and specific to the **D** group.



## TODO List

* Generate large supplementary table
* Generate smaller table of DE genes with fold change
* Generate FC figures

*5/22/2015*

--

My goals here are to try to produce the tables for the publication. The larger supplementary table, the smaller summary table and the fold change figures as well. 

*5/26/2015*

OK, after some testing, I was able to write a script to format the EBseq results to a tab-delimited format. I will then take the formatted results and modify them in excel to give Randy and Tony the condensed results.

> Blade 14: /mnt/iscsi/vnx_gliu_7/rna_seq/RBaldwin

```bash
mkdir tab_output
perl ~/perl_toolchain/sequence_data_scripts/combineEBSeqResultsToTab.pl -c biopsy_D/rsem_D_biopsy.results.condmeans -f biopsy_D/rsem_D_biopsy.results -o tab_output/rsem_biopsy_D.tab

perl ~/perl_toolchain/sequence_data_scripts/combineEBSeqResultsToTab.pl -c biopsy_L1/rsem_L1_biopsy.results.condmeans -f biopsy_L1/rsem_L1_biopsy.results -o tab_output/rsem_biopsy_L1.tab

perl ~/perl_toolchain/sequence_data_scripts/combineEBSeqResultsToTab.pl -c biopsy_L2/rsem_L2_biopsy.results.condmeans -f biopsy_L2/rsem_L2_biopsy.results -o tab_output/rsem_biopsy_L2.tab
```

That was the easy part. I can create a new excel file with multiple spreadsheets to handle this larger supplementary table. Now I just need to manually pick out the portions that are likely to be the most interesting and generate the smaller table that Randy requested. I can do this by selecting the common genes that are present using my VENN comparison script, and then I should be able to select a few (10?) winners for further analysis.

```bash
cp -r tab_output /mnt/nfs/nfs2/dbickhart/rnaseq/

# Now to start the VENN comparison. I want the ENSGENE IDs that are >= 95% confidence, found in both the L1 and L2
# datasets, but are NOT found in the D dataset. I am also going to focus on patterns 2, 3 and 4
cd analysis

mkdir venns
# Starting with pattern 2
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o rsem_L1_biopsy.Pattern2.ensgene rsem_L2_biopsy.Pattern2.ensgene
	File Number 1: rsem_L1_biopsy.Pattern2.ensgene
	File Number 2: rsem_L2_biopsy.Pattern2.ensgene
	Set     Count
	1       18
	2       18

# Whoops! They're more different than I thought! Let's try the other patterns first
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl rsem_L1_biopsy.Pattern3.ensgene rsem_L2_biopsy.Pattern3.ensgene
	File Number 1: rsem_L1_biopsy.Pattern3.ensgene
	File Number 2: rsem_L2_biopsy.Pattern3.ensgene
	Set     Count
	1       5
	1;2     1	# Only one!
	2       112

perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl rsem_L1_biopsy.Pattern4.ensgene rsem_L2_biopsy.Pattern4.ensgene
	File Number 1: rsem_L1_biopsy.Pattern4.ensgene
	File Number 2: rsem_L2_biopsy.Pattern4.ensgene
	Set     Count
	1       210
	1;2     51	# OK, thats far more, and it is the pattern of expression expected
	2       84

# OK, so pattern 4 has the most number of shared genes. Let's save that and work with it. I will also prepare 
# a comparison with the entire >= 95% confidence list. Two comparisons, one small table.
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o rsem_L1_biopsy.Pattern4.ensgene rsem_L2_biopsy.Pattern4.ensgene
mv group_1_2.txt rsem_L1_L2_biopsy.Pattern4.ensgene
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o rsem_L1_L2_biopsy.Pattern4.ensgene rsem_D_biopsy.Pattern4.ensgene
	File Number 1: rsem_L1_L2_biopsy.Pattern4.ensgene
	File Number 2: rsem_D_biopsy.Pattern4.ensgene
	Set     Count
	1       49
	1;2     2
	2       48

	Group: 1 in output file: group_1.txt
	Group: 1;2 in output file: group_1_2.txt
	Group: 2 in output file: group_2.txt

mv group_1.txt venns/rsem_L1_L2_specific.Pattern4.ensgene

# OK, now to repeat the process with the whole gene list
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o rsem_L1_biopsy.95.ensgene rsem_L2_biopsy.95.ensgene
	File Number 1: rsem_L1_biopsy.95.ensgene
	File Number 2: rsem_L2_biopsy.95.ensgene
	Set     Count
	1       176
	1;2     112
	2       172

	Group: 1;2 in output file: group_1_2.txt
	Group: 2 in output file: group_2.txt
	Group: 1 in output file: group_1.txt

mv group_1_2.txt rsem_L1_L2_biopsy.95.ensgene
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o rsem_L1_L2_biopsy.95.ensgene rsem_D_biopsy.95.ensgene
	File Number 1: rsem_L1_L2_biopsy.95.ensgene
	File Number 2: rsem_D_biopsy.95.ensgene
	Set     Count
	1       96
	1;2     16
	2       215

	Group: 1;2 in output file: group_1_2.txt
	Group: 2 in output file: group_2.txt
	Group: 1 in output file: group_1.txt

mv group_1.txt venns/rsem_L1_L2_specific.95.ensgene
cp -r venns /mnt/nfs/nfs2/dbickhart/rnaseq/

```

So, if I focus on pattern specific expression just for Pattern4, I get 49 genes present in L1 and L2 but NOT in D. If I cast a wide net and select all genes, regardless of pattern (L1 U L2 && !D), then I get 96. Let's download them and select the specific information from the table to create a smaller entry.

I'm going to run DAVID on the pattern4 set first, just because I'm curious.

The DAVID results were not strong, but likely because there were so few genes and many fell within "ion binding" categories (there are so many in the genome). I saved the results files, regardless. Now to make the tables.

In order to make the smaller table, I'm going to select the L1/L2 Pattern4 specific transcripts for the table.

```bash
for i in `cat venns/rsem_L1_L2_specific.Pattern4.ensgene`; do grep $i ../tab_output/rsem_biopsy_L*.tab; done > ../tab_output/rsem_L1_L2_specific_pattern4.tab

cp ../tab_output/rsem_L1_L2_specific_pattern4.tab /mnt/nfs/nfs2/dbickhart/rnaseq/tab_output/
```

<a name="dissect"></a>
## Dissecting more useful information

*8/19/2015*

I need to prepare several distinct gene lists so that genes that are barely approaching significance can be identified and tabulated for each dataset. Here are the goals:

* Identify genes above 90% PPDE rather than 95% PPDE
* Isolate Patterns 2, 3 and 4 instead of just Pattern4
* Prepare individual timepoint gene lists instead of just L1_L2
* Prepare pairwise timepoint comparison lists

Much of this will mirror the commands I listed above, albeit on a more extensive scale.

> Blade 14: /mnt/iscsi/vnx_gliu_7/rna_seq/RBaldwin

```bash
# going to create a new folder for zipping
mkdir gene_lists

# Within the "analysis" folder in this directory (/mnt/iscsi/vnx_gliu_7/rna_seq/RBaldwin/analysis) I already
# have the different expression patterns separated, but I think that they're >= 95% PPDE
# I'll have to start from scratch
mkdir gene_lists/raw_gene_lists

# Dry time
perl -lane 'if($F[-1] > 0.90 && $F[-2] ne "\"Pattern1\"" && $F[-2] ne "\"Pattern5\""){$F[-2] =~ s/\"//g; open(OUT, ">> gene_lists/raw_gene_lists/rsem_D_biopsy.$F[-2].ensgene"); ($g) = $F[0] =~ m/\".*(ENSBTAG.+)\"/; print OUT "$g"; close OUT;}' < biopsy_D/rsem_D_biopsy.results
# Lactation 1 time
perl -lane 'if($F[-1] > 0.90 && $F[-2] ne "\"Pattern1\"" && $F[-2] ne "\"Pattern5\""){$F[-2] =~ s/\"//g; open(OUT, ">> gene_lists/raw_gene_lists/rsem_L1_biopsy.$F[-2].ensgene"); ($g) = $F[0] =~ m/\".*(ENSBTAG.+)\"/; print OUT "$g"; close OUT;}' < biopsy_L1/rsem_L1_biopsy.results
# Lactation 2 time
perl -lane 'if($F[-1] > 0.90 && $F[-2] ne "\"Pattern1\"" && $F[-2] ne "\"Pattern5\""){$F[-2] =~ s/\"//g; open(OUT, ">> gene_lists/raw_gene_lists/rsem_L2_biopsy.$F[-2].ensgene"); ($g) = $F[0] =~ m/\".*(ENSBTAG.+)\"/; print OUT "$g"; close OUT;}' < biopsy_L2/rsem_L2_biopsy.results

wc -l gene_lists/raw_gene_lists/*
  233 gene_lists/raw_gene_lists/rsem_D_biopsy.Pattern2.ensgene
    4 gene_lists/raw_gene_lists/rsem_D_biopsy.Pattern3.ensgene
   59 gene_lists/raw_gene_lists/rsem_D_biopsy.Pattern4.ensgene
   21 gene_lists/raw_gene_lists/rsem_L1_biopsy.Pattern2.ensgene
   11 gene_lists/raw_gene_lists/rsem_L1_biopsy.Pattern3.ensgene
  340 gene_lists/raw_gene_lists/rsem_L1_biopsy.Pattern4.ensgene
   23 gene_lists/raw_gene_lists/rsem_L2_biopsy.Pattern2.ensgene
  125 gene_lists/raw_gene_lists/rsem_L2_biopsy.Pattern3.ensgene
  176 gene_lists/raw_gene_lists/rsem_L2_biopsy.Pattern4.ensgene
  992 total

# Not too different from the previous lists! Just a bit larger
# Now let's do the pairwise comparisons

```
Pairwise comparisons will follow the following schema:

* L1 vs D
* L2 vs D
* L1 vs L2

```bash
# Making the directory
mkdir gene_lists/comparisons

# I've written a nested shell loop to process them all
for d in 'L1_vs_D' 'L2_vs_D' 'L1_vs_L2'; 
	do 
	echo $d; 
	for i in 2 3 4; 
		do 
		echo $i; 
		first=`echo $d | cut -d'_' -f1`; 
		second=`echo $d | cut -d'_' -f3`; 
		comp=`echo gene_lists/comparisons/rsem_${d}_pattern${i}_shared_gene.list`; 
		perl /home/dbickhart/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o gene_lists/raw_gene_lists/rsem_${first}_biopsy.Pattern${i}.ensgene gene_lists/raw_gene_lists/rsem_${second}_biopsy.Pattern${i}.ensgene; 
		mv group_1_2.txt $comp; 
	done; 
done

# NOTE: there were allot of comparisons that did not result in a usable file
wc -l gene_lists/comparisons/*
   3 gene_lists/comparisons/rsem_L1_vs_D_pattern4_shared_gene.list
   1 gene_lists/comparisons/rsem_L1_vs_L2_pattern3_shared_gene.list
  59 gene_lists/comparisons/rsem_L1_vs_L2_pattern4_shared_gene.list
   2 gene_lists/comparisons/rsem_L2_vs_D_pattern3_shared_gene.list
   2 gene_lists/comparisons/rsem_L2_vs_D_pattern4_shared_gene.list

# Pattern4, L1_vs_L2 is still the winner here. Pattern2 had no shared gene lists
```

OK, there could also be a list of genes that were above the significance threshold but were classified under different patterns. Let's create a list of intersecting genes regardless of pattern profile. By "regardless of pattern profile" I mean any pattern that is not Pattern1 or Pattern5 still! 

```bash
# combining the previously separated pattern files
cat gene_lists/raw_gene_lists/rsem_D_biopsy.Pattern*.ensgene | sort | uniq > gene_lists/raw_gene_lists/rsem_D_biopsy.allpat.ensgene
cat gene_lists/raw_gene_lists/rsem_L1_biopsy.Pattern*.ensgene | sort | uniq > gene_lists/raw_gene_lists/rsem_L1_biopsy.allpat.ensgene
cat gene_lists/raw_gene_lists/rsem_L2_biopsy.Pattern*.ensgene | sort | uniq > gene_lists/raw_gene_lists/rsem_L2_biopsy.allpat.ensgene

# OK, now to do the comparisons
for d in 'L1_vs_D' 'L2_vs_D' 'L1_vs_L2'; 
	do 
	echo $d; 
	first=`echo $d | cut -d'_' -f1`; 
	second=`echo $d | cut -d'_' -f3`; 
	comp=`echo gene_lists/comparisons/rsem_${d}_allpat_shared_gene.list`; 
	perl /home/dbickhart/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o gene_lists/raw_gene_lists/rsem_${first}_biopsy.allpat.ensgene gene_lists/raw_gene_lists/rsem_${second}_biopsy.allpat.ensgene; 
	mv group_1_2.txt $comp; 
done

wc -l gene_lists/comparisons/rsem_*allpat*.list
  18 gene_lists/comparisons/rsem_L1_vs_D_allpat_shared_gene.list
 104 gene_lists/comparisons/rsem_L1_vs_L2_allpat_shared_gene.list
  19 gene_lists/comparisons/rsem_L2_vs_D_allpat_shared_gene.list
 141 total

# L1_vs_L2 is still the winner!
```

I think that I just need to reformat the excel spreadsheet to just show the significant PPDE genes and then I'm done here.