# Processing Fescue RNA-seq results
---
*5/14/2015*

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