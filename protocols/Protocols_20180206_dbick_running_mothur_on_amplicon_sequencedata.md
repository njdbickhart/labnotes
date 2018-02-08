# Mothur analysis of 16S amplicons
---
*2/6/2018*

These are my notes on how to perform a Mothur analysis on 16S amplicons using Madison's best practices workbook.



## Preparing for the analysis

```bash
mkdir buccal_sop
cp FASTQ_Generation_2017_12_17_10_37_25Z_67404143/*.fastq ./buccal_sop/

wget https://www.mothur.org/w/images/3/32/Silva.nr_v132.tgz
wget http://www.mothur.org/w/images/1/19/Gg_13_8_99.refalign.tgz
wget http://www.mothur.org/w/images/6/68/Gg_13_8_99.taxonomy.tgz

tar -xvf Silva.nr_v132.tgz
tar -xvf Gg_13_8_99.refalign.tgz; tar -xvf Gg_13_8_99.taxonomy.tgz
```

## Filtering and preparing sequence files for assembly/redundancy

> /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/buccal_samples/pilot_project

```bash
/mnt/nfs/nfs2/bickhart-users/binaries/mothur/mothur

make.file(inputdir=buccal_sop, type=fastq)
make.contigs(file=buccal_sop/stability.files, processors=60)

summary.seqs(fasta=buccal_sop/stability.trim.contigs.fasta)
# Because Mothur didn't parse the fastq files correctly, the group names were missing. I'm using the groups file from a previous run
screen.seqs(fasta=buccal_sop/stability.trim.contigs.fasta, group=FASTQ_Generation_2017_12_17_10_37_25Z_67404143/pilot.true.contigs.groups, maxambig=0, maxlength=300, maxhomop=8, processors=60)

system(mv buccal_sop/pilot.true.contigs.good.groups buccal_sop/stability.contigs.good.groups)

summary.seqs(fasta=buccal_sop/stability.trim.contigs.good.fasta)

Using 60 processors.

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       250     250     0       3       1
2.5%-tile:      1       252     252     0       3       20592
25%-tile:       1       253     253     0       4       205918
Median:         1       253     253     0       4       411835
75%-tile:       1       253     253     0       5       617752
97.5%-tile:     1       254     254     0       6       803078
Maximum:        1       298     298     0       8       823669
Mean:   1       252.929 252.929 0       4.42404
# of Seqs:      823669

unique.seqs(fasta=buccal_sop/stability.trim.contigs.good.fasta)

count.seqs(name=buccal_sop/stability.trim.contigs.good.names, group=buccal_sop/stability.contigs.good.groups)
```

## Running the alignment

```bash
align.seqs(fasta=buccal_sop/stability.trim.contigs.good.unique.fasta, reference=silva.nr_v132.align, flip=T)
summary.seqs(fasta=buccal_sop/stability.trim.contigs.good.unique.align, count=buccal_sop/stability.trim.contigs.good.count_table)

# Since I'm using the full db, I need to reduce the size to just the V4 region
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1046    1046    1       0       1       1
2.5%-tile:      13862   23444   252     0       3       20592
25%-tile:       13862   23444   253     0       4       205918
Median:         13862   23444   253     0       4       411835
75%-tile:       13862   23444   253     0       5       617752
97.5%-tile:     13862   23444   254     0       6       803078
Maximum:        43115   43116   297     0       8       823669
Mean:   13862.5 23444.5 252.898 0       4.42365
# of unique seqs:       98772
total # of seqs:        823669

system(mv silva.nr_v132.pcr.align silva.v132.v4.align)
# Realigning
align.seqs(fasta=buccal_sop/stability.trim.contigs.good.unique.fasta, reference=silva.v132.v4.align, flip=T)

summary.seqs(fasta=buccal_sop/stability.trim.contigs.good.unique.align, count=buccal_sop/stability.trim.contigs.good.count_table)

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        5       13      3       0       1       1
2.5%-tile:      8       9582    251     0       3       20592
25%-tile:       8       9582    252     0       4       205918
Median:         8       9582    252     0       4       411835
75%-tile:       8       9582    252     0       5       617752
97.5%-tile:     8       9582    253     0       6       803078
Maximum:        9578    9582    267     0       8       823669
Mean:   10.296  9581.71 251.876 0       4.42353
# of unique seqs:       98772
total # of seqs:        823669

screen.seqs(fasta=buccal_sop/stability.trim.contigs.good.unique.align, count=buccal_sop/stability.trim.contigs.good.count_table, summary=buccal_sop/stability.trim.contigs.good.unique.summary, start=8, end=9582)
summary.seqs(fasta=buccal_sop/stability.trim.contigs.good.unique.good.align, count=buccal_sop/stability.trim.contigs.good.good.count_table)

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        5       9582    217     0       3       1
2.5%-tile:      8       9582    251     0       3       20292
25%-tile:       8       9582    252     0       4       202920
Median:         8       9582    252     0       4       405840
75%-tile:       8       9582    252     0       5       608759
97.5%-tile:     8       9582    253     0       6       791387
Maximum:        8       9582    267     0       8       811678
Mean:   7.9999  9582    251.921 0       4.42846
# of unique seqs:       97400
total # of seqs:        811678

filter.seqs(fasta=buccal_sop/stability.trim.contigs.good.unique.good.align, vertical=T, trump=.)

unique.seqs(fasta=buccal_sop/stability.trim.contigs.good.unique.good.filter.fasta, count=buccal_sop/stability.trim.contigs.good.good.count_table)
```

## Preclustering

```bash
pre.cluster(fasta=buccal_sop/stability.trim.contigs.good.unique.good.filter.unique.fasta, count=buccal_sop/stability.trim.contigs.good.unique.good.filter.count_table, diffs=2, processors=60)

summary.seqs(fasta=buccal_sop/stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=buccal_sop/stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table)

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       550     217     0       3       1
2.5%-tile:      1       550     251     0       3       20292
25%-tile:       1       550     252     0       4       202920
Median:         1       550     252     0       4       405840
75%-tile:       1       550     252     0       5       608759
97.5%-tile:     1       550     253     0       6       791387
Maximum:        1       550     267     0       8       811678
Mean:   1       550     251.926 0       4.41688
# of unique seqs:       35763
total # of seqs:        811678


chimera.uchime(fasta=buccal_sop/stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=buccal_sop/stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t)

remove.seqs(fasta=buccal_sop/stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=buccal_sop/stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table, accnos=buccal_sop/stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.accnos)

summary.seqs(fasta=buccal_sop/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=buccal_sop/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.count_table)

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       550     217     0       3       1
2.5%-tile:      1       550     251     0       3       19978
25%-tile:       1       550     252     0       4       199774
Median:         1       550     252     0       4       399548
75%-tile:       1       550     252     0       5       599322
97.5%-tile:     1       550     253     0       6       779118
Maximum:        1       550     267     0       8       799095
Mean:   1       550     251.927 0       4.41876
# of unique seqs:       25680
total # of seqs:        799095


classify.seqs(fasta=buccal_sop/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=buccal_sop/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.count_table, reference=silva.v132.v4.align, taxonomy=silva.nr_v132.tax, cutoff=80)

remove.lineage(fasta=buccal_sop/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=buccal_sop/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.count_table, taxonomy=buccal_sop/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.taxonomy, taxon=unknown;-Archaea;-Eukaryota;)

summary.seqs(fasta=buccal_sop/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=buccal_sop/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.count_table)

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       550     220     0       3       1
2.5%-tile:      1       550     251     0       3       19667
25%-tile:       1       550     252     0       4       196669
Median:         1       550     252     0       4       393338
75%-tile:       1       550     252     0       5       590006
97.5%-tile:     1       550     252     0       6       767008
Maximum:        1       550     266     0       8       786674
Mean:   1       550     251.913 0       4.40849
# of unique seqs:       25330
total # of seqs:        786674


system(cp buccal_sop/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta buccal_sop/stability.final.fasta)

system(cp buccal_sop/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.count_table buccal_sop/stability.final.count_table)
```

## Defining OTUs

```bash
dist.seqs(fasta=buccal_sop/stability.final.fasta)
cluster.split(column=buccal_sop/stability.final.dist, count=buccal_sop/stability.final.count_table, method=average)

make.shared(list=buccal_sop/stability.final.an.unique_list.list, count=buccal_sop/stability.final.count_table, label=0.03)

classify.seqs(fasta=buccal_sop/stability.final.fasta, count=buccal_sop/stability.final.count_table, template=gg_13_8_99.fasta, taxonomy=gg_13_8_99.gg.tax, cutoff=80)

classify.otu(list=buccal_sop/stability.final.an.unique_list.list, taxonomy=buccal_sop/stability.final.gg.wang.taxonomy, count=buccal_sop/stability.final.count_table, label=0.03, cutoff=80, basis=otu, probs=F)
```

## Generating coverage estimates and general QC

```bash
rarefaction.single(shared=buccal_sop/stability.final.an.unique_list.shared, processors=60)

summary.single(shared=buccal_sop/stability.final.an.unique_list.shared, label=0.03, calc=nseqs-sobs-coverage)
```

