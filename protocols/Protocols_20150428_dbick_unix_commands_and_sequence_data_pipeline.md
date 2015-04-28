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