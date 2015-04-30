# Processing National DB NextSeq data
---
*04/28/2015*

### Gordon's notes and useful unix commands

This is the name of the editor that I'm using: [MarkdownPad2.](http://markdownpad.com/)
 
> Here are a few useful commands that I used to get basic information about the filesystem and datafiles.

```bash
# This command prints your current working directory (really helpful!)
pwd

# This command shows you all of the attached mounts
df

# How to calculate the number of reads from a fastq without unzipping
gunzip -c Undetermined_S0_L001_R1_001.fastq.gz | wc -l
# 11919260 / 4 = 2,978,915

# this does a word count on all fastqs that begin with 16 and prints the number of reads out
for i in 16*.gz; do gunzip -c $i | wc -l; done | perl -e '$c = 0; while(<>){chomp; $c += $_;} print ($c / 4); print "\n";'
# 188,433,212 

# Check the length of the reads in a fastq
perl -e '$s = "AGGAGCCTGGTAGGCTGCAGTCCATGGGGTCGCTAAGAGTCAGACACGACTGAGAGAAATAACATTCAAATTAAACTTTCATGAATTTGAGAAAGAAATGGAGACAGACTCAAAGAGTCTTGAATGGAGATTCAAAAGGACGGGGGAGACG"; print length($s); print "\n";'
151

# List the home directory
ls ~

# Reading files using Unix commands
less my_text_file.txt

# Using commands in tandem with unix "pipes"
grep 'my_pattern' ./*files | less
	# Searches for "my_pattern" in the files and allows you to read the output using "less"

```

#### Using Unix "less" to read files

Here is a little primer on useful commands in "less" that will help you read text files. All of these commands can be entered when you're in the interactive "less" screen.

| Command | Description|
| --- | :---: |
| q | Exits less |
| arrow up| Go up one line |
| arrow down| Go down one line |
| page up | Go up one page |
| page down | Go down one page |
| g | Go to the top of the file |
| G | Go to the bottom of the file | 
| **10**G | Go down to line number **10** from the top of the page |
| /(pattern) | searches for the first occurrence of the pattern|
| n | if used after "pattern searching, returns next pattern |
| N | Goes to the previous pattern |

#### Running longer commands

Sometimes commands take longer than you expect, or you might accidentally run the wrong command on your files! If you want to prematurely terminate a command, just type:

> Control + c

This is the universal Unix termination command.

Another way to terminate your commands is to exit your console session -- once you leave, the command terminates as well. This can also be a bad thing, as you might run programs that will take weeks to run, and you've gotta go home at some point! I use a Unix program called "screens" in order to manage long running jobs, as you can run a job on a "screen" and it will continue to run until one of the following things happens:

  1. The computer crashes
  2. You terminate the program

Here is how you run a screen session:

```bash
# this creates a new "screen." You can think of it like a separate console window
screen -a
```

Within the screen, you will still be in the same working directory, but now you can execute your program without worrying about it terminating when you leave for home. Another fringe benefit is when the internet goes down -- your job will still be running on the screen so long as it did not rely on internet access.

To leave the screen, type:

> Control + a (and then) d

When you want to reconnect to the screen, type:

```bash
screen -r
```

If you have more than one screen, you will be given a list of screens, and then you can select which screen to return to by typing it after the "-r" argument, like so:

```bash
# If given a choice of screens
screen -r 28738-screen.tty
```

To recap the use of screens:

| Command | Description |
| :---: | :---: |
| screen -a | Create new screen in the current directory |
| Control + a + d | leave the current screen |
| screen -r | Return to a previously created screen |

To kill a screen, type "exit" within that screen to terminate it.

--
### Data processing

#### Preliminary statistics

It is important to check data input before alignment.
One of the reasons for this is that you can actually "pre-treat" the fastq files (trimming, excluding reads) that will decrease your processing time and improve your data quality.

I use FastQC for fastq file data assessment. You can download [Fastqc here.](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

> pwd: **/mnt/cifs/bickhart-qnap/NextSeq/150403_NS500432_0006_AH1547BGXX/H1547BGXX/DBickhart-NatDB-2015-04-03**

```bash
# you can add as many fastq files at the end of this command as you want
# Runs fastqc analysis on the file
~/FastQC/fastqc ./16_S1_L001_R1_001.fastq.gz

```

Obviously, it takes a long time to process all of the fastq files, so you can actually use multiple cpu cores (termed "threading") in order to speed up the process.

```bash
# This command runs fastqc, but with 3 cpu "threads" or "cores" being used by the program
# The "-t" option controls the number of threads being used.
~/FastQC/fastqc -t 3 ./16_S1_L001_R1_001.fastq.gz
```

Also, you can run multiple files through fastqc side-by-side.

```bash
# This will process 10 files in parallel in the directory (threads = 10) and it will process all files in the directory.
# The '*' serves as a wildcard, indicating that all of the files must end in a ".gz" extension
~/FastQC/fastqc -t 10 *.gz
```

*04/29/2015*

--
### Testing a pipeline script

I created a larger pipeline that automates the quality control, alignment and SNP/INDEL calling of fastq files by queuing up a large list of them to run on several CPU processor cores at once. Here are my testing notes while I debug (numerous) issues with the prototype of this pipeline. 

#### Pipeline program prerequisites

1. Perl modules
  * Mouse
  * Mouse::NativeTraits
  * Forks::Super
  * namespace::autoclean
  * My [perl modules](https://github.com/njdbickhart/perl_toolchain)
2. [BWA](http://bio-bwa.sourceforge.net/)
3. [Samtools](http://samtools.sourceforge.net/)
4. [BCFtools](http://samtools.github.io/bcftools/)
5. [Fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
6. [Picard tools](http://broadinstitute.github.io/picard/)
7. [SNPeff](http://snpeff.sourceforge.net/)

#### Design philosophy

This pipeline was created to automate the processing of fastq files from beginning to (near) end. The user just enters a "spreadsheet" tab-delimited file containing the fastq files paired in rows, and the pipeline will align them to a reference genome, create bam files and call SNPs/INDELs on them. It will also eventually annotate the variants using SnpEFF, and "grep out" the predicted variants with high functional impact for the end user. 

#### Setting up a simulation test run to debug the code

I need to create a very small dataset to process the code so that I can identify any errors that I need to correct from the start. I will do this using the NextSeq fastq files in the qnap NAS folder for our latest run.

> Blade14: /mnt/cifs/bickhart-qnap/NextSeq/150403_NS500432_0006_AH1547BGXX/H1547BGXX/DBickhart-NatDB-2015-04-03

```bash
# First, to get a subset of reads from the fastq files that I want to work with
gunzip -c 16_S1_L001_R1_001.fastq.gz | head -n 4000 | gzip > sample_fastq.1.1.fq.gz
gunzip -c 16_S1_L001_R2_001.fastq.gz | head -n 4000 | gzip > sample_fastq.1.2.fq.gz

# That's the first pair, I'm going to make one more pair to test how the pipeline handles multiple samples
gunzip -c 25_S4_L001_R1_001.fastq.gz | head -n 4000 | gzip > sample_fastq.2.1.fq.gz
gunzip -c 25_S4_L001_R2_001.fastq.gz | head -n 4000 | gzip > sample_fastq.2.2.fq.gz

# Now to create the "sample spreadsheet" that I will use in the pipeline
# I am creating a variable "cwd" to store the current working directory to make this more readable
cwd=`pwd`; echo -e "$cwd/sample_fastq.1.1.fq.gz\t$cwd/sample_fastq.1.2.fq.gz\tlib1\tsamp1\n" > simulation_spreadsheet.tab

# The '>>' redirect command appends to the end of the file
cwd=`pwd`; echo -e "$cwd/sample_fastq.2.1.fq.gz\t$cwd/sample_fastq.2.2.fq.gz\tlib2\tsamp2\n" >> simulation_spreadsheet.tab

head simulation_spreadsheet.tab
	/mnt/cifs/bickhart-qnap/NextSeq/150403_NS500432_0006_AH1547BGXX/H1547BGXX/DBickhart-NatDB-2015-04-03/sample_fastq.1.1.fq.gz     /mnt/cifs/bickhart-qnap/NextSeq/150403_NS500432_0006_AH1547BGXX/H1547BGXX/DBickhart-NatDB-2015-04-03/sample_fastq.1.2.fq.gz        lib1       samp1
	/mnt/cifs/bickhart-qnap/NextSeq/150403_NS500432_0006_AH1547BGXX/H1547BGXX/DBickhart-NatDB-2015-04-03/sample_fastq.2.1.fq.gz     /mnt/cifs/bickhart-qnap/NextSeq/150403_NS500432_0006_AH1547BGXX/H1547BGXX/DBickhart-NatDB-2015-04-03/sample_fastq.2.2.fq.gz        lib2       samp2
```

Now I'm going to generate a configuration file that states where all of the program executables and other necessary files are located. My pipeline checks for a hidden configuration file within the home directory by default: ~/.mergedpipeline.cnfg

*04/30/2015*

--

I'm going to run the simulation pipeline on this test dataset on a different drive. Reasons: If I read from the NAS and write to the other drive, that spreads the work on the disks and increases the speed by which the pipeline can run. Also, it won't wear out the NAS prematurely!

> Blade14: /mnt/iscsi/vnx_gliu_7/test

```bash
# Getting the pipeline usage text
perl ~/perl_toolchain/sequence_data_pipeline/runMergedBamPipeline.pl
	#perl /home/dbickhart/perl_toolchain/sequence_data_pipeline/runMergedBamPipeline.pl --fastqs <spreadsheet to be processed> --output <base output folder> --reference <reference genome fasta>
	#Arguments:
    	    #--fastqs        A tab-delimited spreadsheet with fastq file locations and read group info
        	#--output        The base output folder for all bams and vcfs
        	#--reference     The reference genome fasta file

	#Optional:
    	    #--config        Override default configuration (searches for ~/.mergedpipeline.cnfg by default)
    	    #--coords        Reference genome fasta coordinates for threaded SNP calling
    	    #--threads       Number of threads to fork (default: 1)

# OK, now I'm going to run the pipeline on the simulation dataset
perl ~/perl_toolchain/sequence_data_pipeline/runMergedBamPipeline.pl \ 
	--fastqs /mnt/cifs/bickhart-qnap/NextSeq/150403_NS500432_0006_AH1547BGXX/H1547BGXX/DBickhart-NatDB-2015-04-03/simulation_spreadsheet.tab \
	--output sim --reference /mnt/iscsi/vnx_gliu_7/reference/umd3_kary_unmask_ngap.fa \
	--coords /mnt/iscsi/vnx_gliu_7/reference/samtools_chr_segs.txt \
	--threads 2

```

Let the bug-fest begin! I fixed several bugs related to command line options, and now I'm getting bugs related to fastqc:

> Unknown option: coords
> Can't locate object method "INFO" via package "simpleLogger" at /home/dbickhart/perl_toolchain/sequence_data_pipeline/runMergedBamPipeline.pl line 288, <IN> line 1.

> Can't locate object method "INFO" via package "simpleLogger" at /home/dbickhart/perl_toolchain/sequence_data_pipeline/runMergedBamPipeline.pl line 288, <IN> line 2.

> Use of uninitialized value in concatenation (.) or string at /home/dbickhart/perl_toolchain/personal_modules/fastqcParser.pm line 73.
> Use of uninitialized value in concatenation (.) or string at /home/dbickhart/perl_toolchain/personal_modules/fastqcParser.pm line 73.
Error with pipeline! [FQCPARSE] - Error accessing fastqc_data.txt in folder:  for file: /mnt/cifs/bickhart-qnap/NextSeq/150403_NS500432_0006_AH1547BGXX/H1547BGXX/DBickhart-NatDB-2015-04-03/sample_fastq.1.1.fq.gz!

I'm guessing that fastqc is not creating the folder using the naming conventions that I designed in the module.

```bash
ls /mnt/cifs/bickhart-qnap/NextSeq/150403_NS500432_0006_AH1547BGXX/H1547BGXX/DBickhart-NatDB-2015-04-03/
```

> sample_fastq.1.1.fq_fastqc
> sample_fastq.1.1.fq_fastqc.zip
> sample_fastq.1.1.fq.gz

Yup! Looks like I need to rewrite the name parsing for that module.

Damn, I can't use forks with objects, so I'll have to use threads instead. I'll rewrite the code to do this later.