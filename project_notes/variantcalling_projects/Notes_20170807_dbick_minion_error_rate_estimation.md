# Generating MinION read quality stats
---
*8/7/2017*

## MinION error rate estimation

These are my notes for rapidly generating an assessment of MinION read error rates using a Lambda control run. I will be calculating a per-read error rate using alignment of each read to the lambda reference genome, calculating reference edit distance and then dividing it by the read length.

> linux box: /home/dbickhart/share/nanopore/lambda_pass_test

```bash
# Preparing reference fasta downloaded from NCBI
bwa index lambda_phage.fasta

for i in workspace/*.fastq; do name=`basename $i | cut -d'.' -f1`; bwa mem lambda_phage.fasta $i | samtools sort -T $name.temp -o $name.sorted.bam -; done
for i in *.bam; do echo $i; samtools index $i; done
samtools merge lambda_run_merged.sorted.bam fastq_runid_25e3bb480d066d56689dd027a4fc35bb2a14eba5_0.sorted.bam fastq_runid_25e3bb480d066d56689dd027a4fc35bb2a14eba5_1.sorted.bam fastq_runid_25e3bb480d066d56689dd027a4fc35bb2a14eba5_2.sorted.bam fastq_runid_25e3bb480d066d56689dd027a4fc35bb2a14eba5_3.sorted.bam fastq_runid_25e3bb480d066d56689dd027a4fc35bb2a14eba5_4.sorted.bam fastq_runid_25e3bb480d066d56689dd027a4fc35bb2a14eba5_5.sorted.bam fastq_runid_25e3bb480d066d56689dd027a4fc35bb2a14eba5_6.sorted.bam fastq_runid_25e3bb480d066d56689dd027a4fc35bb2a14eba5_7.sorted.bam fastq_runid_25e3bb480d066d56689dd027a4fc35bb2a14eba5_8.sorted.bam fastq_runid_25e3bb480d066d56689dd027a4fc35bb2a14eba5_9.sorted.bam fastq_runid_38f77a1dc2ede5ed24748b180d768345392284d4_0.sorted.bam

samtools index lambda_run_merged.sorted.bam

# Now to start estimating error. I'll use the NM:I tag of the bam to estimate error rates as that is flat error rate
python3 ../../programs_source/python_toolchain/sequenceData/calcErrorRateBam.py lambda_run_merged.sorted.bam
#	BAM file:lambda_run_merged.sorted.bam	Num reads: 32815	Error High: 0.22826086956521738	Error Low: 0.0	Avg Error: 0.09429249393467176	Stdev Error: 0.06351189018395148
```

Then to take the output and plot a histogram of error rates in R:

```R
colnames(values) <- c("value")
library(ggplot2)
ggplot(data=values, aes(value)) + geom_histogram() + scale_x_continuous(limits = c(0, 0.22)) + xlab("Error ratio") + ylab("Count") + ggtitle("MinION Lambda Alignment Error Ratios") + theme(plot.title = element_text(hjust = 0.5))
dev.copy2pdf(file="minion_error_rate_simulation.pdf", useDingbats=FALSE)
```


Here is the python script I wrote to calculate the read stats for the alignment file:

#### calcErrorRateBam.py
```python
#!/usr/bin/python3
"""
This is a script designed to process a large bam file and generate statistics on read mapping error rates
The goal is to estimate error rates of long reads mapped against a reference genome (useful for Nanopore data)
"""

import sys
import re
import math
import subprocess
#from multiprocessing import Pool
USAGE = "Calculate edit distance error rate for aligned reads.\n python3 <input bam file>"

# Compiled regular expression to ensure capture of MDZ tag info
NMITAG = re.compile("NM:i")
NMIVAL = re.compile("NM:i:(\d+)")

"""
For each read, count the number of NMI tag differences and return a 
ratio against length of read
"""
def MapDiff(L): 
    segs = L.decode("utf-8").split("\t")
    if int(segs[1]) & 2048 == 2048:
        return -1
        
    if len(segs) < 12:
        return -1
    mo = re.match(NMIVAL, segs[11]);
    if mo:
        val = mo.group(1)
    
        seqlen = len(segs[9])
        return int(val) / seqlen
    else:
        return -1
    
"""
Reduce all ratios by taking the average and calculating the stdev
Returns the [avg, stdev]
"""
def ReduceAvgStd(L):
    num = len(L)
    
    c = 0
    high = 0
    low = 1
    for i in L:
        if i > high:
            high = i
        if i < low:
            low = i
        c += i
        
    avg= 0
    if num > 0:
        avg = c / num
    else:
        avg = 0
        
    # Now for the stdev calculation
    ss = 0
    for i in L:
        ss += (i - avg) * (i - avg)
    
    stdev = 0
    if num > 0: 
        stdev = math.sqrt(ss / num)
    else:
        stdev = 0
    
    return [num, high, low, avg, stdev]

if __name__ == '__main__':
    if(len(sys.argv) < 2):
        print(USAGE)
        sys.exit(1)
    
    f = subprocess.Popen(['samtools', 'view', sys.argv[1]], \
                        stdout=subprocess.PIPE)
    [num, high, low, avg, stdev] = ReduceAvgStd(list(filter(lambda x: x != -1, map(MapDiff, f.stdout.readlines()))))
    
    # Uncomment the following line if you want to get the list of ratios
    """
    with open("ratio_count.list", "w") as o:
        for l in list(filter(lambda x: x != -1, map(MapDiff, f.stdout.readlines()))):
            o.write("{}\n".format(l))
    """
        
    print("BAM file:{}\tNum reads: {}\tError High: {}\tError Low: {}\tAvg Error: {}\tStdev Error: {}"\
         .format(sys.argv[1], num, high, low, avg, stdev))
```


## MinION read length distribution

I will calculate this histogram from the first successful Yu and Morrison run just to save some time.

> Linux_box: /home/dbickhart/share/metagenomics/pilot_project

```bash
# Generating a quick list of read lengths
for i in nanopore/yu_and_morrison_3/*.fastq; do perl -e 'while(<>){$s = <>; $p = <>; $q = <>; print length($s); print "\n";}' < $i; done > yu_and_morrison_3_readlengths.list

wc -l yu_and_morrison_3_readlengths.list 
#	1268557 yu_and_morrison_3_readlengths.list
```

Now to load the values into R and to generate a histogram.

```R 
values <- read.delim("yu_and_morrison_3_readlengths.list", header=FALSE)
colnames(values) <- c("ReadLength")

ggplot(data=values, aes(ReadLength)) + geom_histogram() + scale_y_log10() + scale_x_continuous(labels = scales::comma) + ggtitle("MinION Rumen Sample Read Lengths") + ylab("Log10 Read Counts") + xlab("Read Length (bp)") + theme(plot.title = element_text(hjust = 0.5))
dev.copy2pdf(file="minion_readlength_rumen_count.pdf", useDingbats=FALSE)
```
