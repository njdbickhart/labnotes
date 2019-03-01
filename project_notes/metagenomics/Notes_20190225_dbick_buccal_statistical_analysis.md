# Buccal data analysis
---
*2/25/2019*

These are my notes on the statistical analysis of the Buccal swab OTU data with the overall goal of identifying which OTU are most likely part of the "oral-only" microbial community. 

## Table of Contents

## Preparing the data for analysis

I ran into some trouble getting Mothur to consistently analyze the samples I wanted to process. The problem was that the distance matrix I tried to calculate was absolutely huge! I am going to automate this via a script I've written and then use a few metrics to reduce the number of comparisons. 

I also want to include all of the samples taken thus far. The more the better. I want to include rumen contents and buccal samples to better train my models.

First, let's make sure all of the data is in the right place and is properly accounted for:

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/buccal

```bash
mkdir combined_fastqs
mv buccal_samples/evergreen_view/combined_fastq/*.fastq ./combined_fastqs/
mv buccal_samples/amelie_data/4_26_18/*.fastq ./combined_fastqs/
mv buccal_samples/pilot_project/buccal_sop/*.fastq ./combined_fastqs/

ls combined_fastqs/ | wc -l
248

mkdir tax_files
cp buccal_samples/pilot_project/silva.v132.v4.align ./tax_files/
cp buccal_samples/pilot_project/silva.nr_v132.tax ./tax_files/
cp buccal_samples/pilot_project/gg_13_8_99.fasta ./tax_files/
cp buccal_samples/pilot_project/gg_13_8_99.gg.tax ./tax_files/

# Now to make sure all of the fastq files are represented in the input file before queueing this up
perl -ne '$_ =~ s/4_26_18\///g; print $_;' < buccal_samples/amelie_data/true_listings.files > combined_fastqs/combined_fq_total.files
cat buccal_samples/evergreen_view/combined_fastq/db1.files >> combined_fastqs/combined_fq_total.files
perl -ne '$_ =~ s/FASTQ_Generation_2017-12-17_10_37_25Z-67404143\///g; print $_;' < buccal_samples/pilot_project/FASTQ_Generation_2017_12_17_10_37_25Z_67404143/pilot.true.files >> combined_fastqs/combined_fq_total.files

# Now to queue it up. The mothur batch file is: mothur_script.batch; The mothur shell wrapper is mothur_combined_run.sh

# Argh! The "input" directory setting appears to be useless! I need to put paths for all of the fastq files in the "files" file
perl -lane 'print "$F[0]\tcombined_fastqs/$F[1]\tcombined_fastqs/$F[2]";' < combined_fq_total.files > temp
mv temp combined_fq_total.files

sbatch mother_combined_run.sh

# Good and bad news
# The good news: the script ran with only one error!
# The bad news: the output files are in the working directory folder(really???), and I mispelled "label" 
# Cleaning up for a new run

rm combined_fq_total.contigs.* combined_fq_total.scrap.contigs.* combined_fq_total.trim.contigs.* combined_fq_total.filter
rm tax_files/*.8mer tax_files/*.numNonZero tax_files/*.prob tax_files/*.tree*

sbatch mother_combined_run.sh

# OK, that worked! Lots of warnings though from sequences that couldn't be classified.
mv *.rabund ./combined_output/
mv *.map ./combined_output/

# Cleaning up a few things and reducing file name lengths
mv combined_fq_total.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist combined_fq_total..final.dist
mv combined_fq_total.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta combined_fq_total..final.fasta
mv combined_fq_total.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.gg.wang.taxonomy combined_fq_total.gg.wang.taxonomy
mv combined_fq_total.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.gg.wang.tax.summary combined_fq_total.gg.wang.tax.summary
mv combined_fq_total.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.groups.summary combined_fq_total.goods.coverage.summary
mv combined_fq_total.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.groups.rarefaction combined_fq_total.rarefaction
mv combined_fq_total.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared combined_fq_total.final.prenorm.shared
mv combined_fq_total.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.10.subsample.shared combined_fq_total.final.subsample.shared
mv combined_fq_total.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.10.cons.taxonomy combined_fq_total.final.subsample.taxonomy
mv combined_fq_total.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.10.cons.tax.summary combined_fq_total.final.subsample.tax.summary

# moving the rest of the files
mv combined_fq_total.trim.contigs.good.* ./combined_output/
```

Taking the rarefaction data to R first.

```R
library(tidyr)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)

rare <- read_tsv(file="combined_fq_total.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.groups.rarefaction") %>% seleParsed with column specification:hci-")) %>% gather(-numsampled, key=sample, value=sobs) %>% mutate(sample=str_replace_all(sample, pattern="0.10-"cols(lacement="")) %>% drop_na()
pdf(file="combined_fq_rarefaction_plot.pdf", useDingbats=FALSE)
ggplot(rare, aes(x=numsampled, y=sobs, group=sample)) + geom_line() + theme_bw()
dev.off()

# The plot showed some samples were far more dense than expected and that we had not reached saturation
```


Now I'm going to dabble in Python for a random forest classifier experiment. I expect this to be very basic and not very informative. Consider this an exploration of the data.

Preparing files for input:

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/buccal

```bash
perl -lane 'print "$F[0]\t$F[3]";' < buccal_samples/evergreen_view/sample_list_categories.tab > combined_fq_list_categories_simple.tab
```

And reading data into Python3 for random forest analysis.


```python
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics

data = pd.read_table('combined_fq_total.final.prenorm.shared', delim_whitespace=True, header=0)
info = pd.read_table('combined_fq_list_categories_simple.tab', delim_whitespace=True, header=0)

# Map the categories to the data groups
data['class'] = data['Group'].map(info.set_index('Entry')['Type'])

# drop unnecessary columns
data = data.iloc[:,3:]

X = data.iloc[:,:-1]
Y = data.iloc[:,-1]

# Test set of 30% of the original input
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.3)

# I am going to need to look at the API for this more closely. I'm going to just test this out first using the 
# default parameters that I saw in the tutorial online.
clf = RandomForestClassifier(n_estimators=100)
clf.fit(X_train, Y_train)

Y_pred = clf.predict(X_test)
print(f'Accuracy:\t{metrics.accuracy_score(Y_test, Y_pred)}')

# OK, not bad, but still I'm skeptical of the results
f_importance = pd.Series(clf.feature_importances_, index=X.columns).sort_values(ascending=False)
f_importance.head(n=25)
Otu00009    0.021439
Otu00004    0.019115
Otu00066    0.018132
Otu00071    0.014613
Otu00069    0.013296
Otu00239    0.012755
Otu00001    0.012001
Otu00024    0.011650
Otu00242    0.010593
Otu00032    0.010217
Otu00140    0.009774
Otu00030    0.009761
Otu00010    0.009684
Otu00002    0.009593
Otu00012    0.009104
Otu00044    0.008123
Otu00060    0.007742
Otu00035    0.007671
Otu00011    0.007460
Otu00034    0.007114
Otu00075    0.006972
Otu00087    0.006930
Otu00064    0.006913
Otu00021    0.006794
Otu00105    0.006343
...

# It looks like there are no huge influencers, but OTU00009 looks like the top candidate
```
