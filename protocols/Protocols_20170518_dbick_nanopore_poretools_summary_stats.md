# Nanopore poretools summary stats
---
*5/18/2017*

These are some brief notes on how to generate quick summary stats on nanopore data using the poretools package.

## Table of Contents
* [Summary plots of reads](#summary)
* [Winner, winner, chicken dinner](#winner)
* [Read squiggle plots](#squiggle)
* [Read length histograms](#hist)
* [Nanopore occupancy](#occupancy)
* [Quality score data](#quality)
* [Extracting fastqs](#fastqs)


<a name="summary"></a>
## Summary plots of reads

Generating a plot of basepairs over time:

```bash
# Saving as a PDF is an option as well
poretools yield_plot --plot-type basepairs --saveas tony_bp_reads.png tonys_reads/
```

<a name="winner"></a>
## Winner, winner, chicken dinner

Selecting the largest read from the group.

```bash
# Note: it dumps all output to stdout
poretools winner --type all tonys_reads/
>channel_22_c0686f97-e9a5-435d-93c1-f67ef7e6629b_template tonys_reads/pass/5/ARSWIMSN4DB1368_20170512_FNFAE32444_MN19849_sequencing_run_Fibro_UWT1_13838_ch22_read10813_strand.fast5

tail -n 1 tonys_reads/longest_read.fa | perl -lane 'print length($F[0]);'
	1,113,569   <- the longest read is 1 megabase, but it looks like it's a repetitive junk read
```
<a name="squiggle"></a>
## Read squiggle plots

Generating a list of "events" for each read file. Each "event" is a transition of electric signal that corresponds to a read/set of reads.

```bash
# This must not be working properly right now
poretools squiggle --saveas png tonys_reads/pass/5/ARSWIMSN4DB1368_20170512_FNFAE32444_MN19849_sequencing_run_Fibro_UWT1_13838_ch98_read6334_strand.fast5
	WARNING:poretools:Could not extract template events for read: tonys_reads/pass/5/	ARSWIMSN4DB1368_20170512_FNFAE32444_MN19849_sequencing_run_Fibro_UWT1_13838_ch98_read6334_strand.fast5.
```

<a name="hist"></a>
## Read length histograms

Generating a histogram of read lengths from a file.

```bash
poretools hist --saveas tony_bp_hist.png tonys_reads/

# better control on the output plot
poretools hist --max-length 65000 --num-bins 30 --saveas tony_bp_hist.png tonys_reads/
```

<a name="occupancy"></a>
## Nanopore occupancy

Generating a site map for pores showing how much data passed through each.

```bash
# This looks like it is currently broken
poretools occupancy --saveas tony_pore_occupancy.png tonys_reads/
```
<a name="quality"></a>
## Quality score data

Generate a plot of read position dependent quality scores.

```bash
# NOTE: png format doesn't work and PDF format is apparently super-verbose, leading to huge files!
poretools qualpos --max-length 60000 --bin-width 1000 --saveas tony_qual_score_whisker.pdf tonys_reads/

# Apparently neither is a jpg supported!
poretools qualpos --max-length 60000 --bin-width 6000 --saveas tony_qual_score_whisker.jpg tonys_reads

# Neither plot is sufficient, so let's generate the quality score tab delimited data for review
poretools qualdist tonys_reads/ > tony_qualdist_output.tab
```

<a name="fastqs"></a>
## Extracting fastqs

Finally, the most useful tool in "poretools."

```bash
# NOTE: data is all dumped to STDOUT
poretools fastq tonys_reads/ > tonys_fastq_minknow_data.fq
```
