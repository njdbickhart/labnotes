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

## Preliminary Random Forest Analysis

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

OK, new strategy. I want to remove the "controls" as these are likely contaminants that we don't want to keep. I will also drop the OTUs associated with the controls so that the data isn't used for training and prediction.

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

# Checking sum of OTU counts per swab/rumen type
groupsum = data.groupby('class').sum()

# This is very simple, but I'm trying to remove OTUs that are present ONLY in the control swabs
drop = [ (groupsum.iloc[0,k] > 0 and groupsum.iloc[1,k] + groupsum.iloc[2,k] <= 0) for k in range(len(groupsum.columns))]
# There were 575 OTUs that fit the bill here
groupsum.loc[:,drop]
         Otu00287  Otu00408  Otu00439  Otu00576  Otu00668  Otu00780    ...     Otu11853  Otu11860  Otu11877  Otu11927  Otu11954  Otu11975
class                                                                  ...
control       169        61        53        31        22        15    ...            1         1         1         1         1         1
rumen           0         0         0         0         0         0    ...            0         0         0         0         0         0
swab            0         0         0         0         0         0    ...            0         0         0         0         0         0

# OK, that required some wrangling but it's done. Now to reformat the data
drop.append(False) # to make sure that we don't drop the category!
filtered = data.drop(data.columns[drop], axis=1)
filtered = filtered[filtered['class'] != 'control'] # removing the controls so that they don't interfere in the classification

X = filtered.iloc[:,:-1]
Y = filtered.iloc[:,-1]

X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.3)
clf = RandomForestClassifier(n_estimators=100)
clf.fit(X_train, Y_train)

Y_pred = clf.predict(X_test)
print(f'Accuracy:\t{metrics.accuracy_score(Y_test, Y_pred)}')
# Accuracy:       0.9090909090909091

f_importance = pd.Series(clf.feature_importances_, index=X.columns).sort_values(ascending=False)
f_importance.head(n=25)
Otu00289    0.028751
Otu00352    0.025927
Otu00107    0.018546
Otu00056    0.017873
Otu00264    0.017768
Otu00042    0.016323
Otu00905    0.015268
Otu00063    0.015009
Otu00123    0.014952
Otu00004    0.014396
Otu00189    0.013485
Otu00001    0.013006
Otu00086    0.011654
Otu00381    0.011640
Otu00073    0.011589
Otu00219    0.011582
Otu00432    0.011373
Otu00454    0.011039
Otu00116    0.009815
Otu00033    0.009608
Otu00170    0.009255
Otu00140    0.009141
Otu00164    0.009141
Otu00192    0.008477
Otu00110    0.008451

# Again, many features contributed towards the classifier. Let's see how the data looks for these OTU.
top25 = list(f_importance.head(n=25).index)
top25.append('class')

# Checking to see how the features differentiate the data
filtered[top25].groupby('class').mean()
       Otu00289  Otu00352   Otu00107    Otu00056  Otu00264    ...       Otu00170   Otu00140   Otu00164   Otu00192   Otu00110
class                                                         ...
rumen  7.222222  3.833333  54.611111  134.611111  9.388889    ...      20.722222  29.944444  23.111111  15.111111  45.333333
swab   0.391304  0.304348   1.989130    7.934783  0.326087    ...       0.663043   1.597826   0.510870   0.880435   3.673913

# Looking to see how many OTU are likely oral related
filtered[top25].groupby('class').mean().transpose()
class          rumen         swab
Otu00004    0.666667  1231.739130 # Pasteurellaceae 
Otu00001    6.555556  3860.934783 # Pseudomonas veronii, which doesn't make much sense

# Let's save this data for future plotting
f_importance.to_csv("otu_feature_importance.tab", sep='\t', index_label = 'OTUName')
filtered.to_csv("filtered_otu_table_nocontrol.tab", sep='\t', index=False)
```

I think that I have the basic idea down. There are some hyperparameters I need to tweak (i.e. random forest permutation counts), but the pipeline is likely to remain similar. I do think that more samples would be beneficial, even if they are left out of the initial training and used as a test set for prediction accuracy measurement. Still, let's see if I can formulate the elements that I need for a manuscript.

## Manuscript planning

#### Necessary figures

* At least one decision tree to show the model and decisions
* Rarefaction curves
* A heat map showing relative abundance of taxa pre-normalization
* Feature contribution and Feature importance plots
* Post normalization heatmaps
* Comparative PCA plots showing differences in data clustering pre- and post-normalization
* k-means clustering to identify potential sources of substructure (ie. pre and post rumination swabs)

#### Necessary tables

* OTU counts pre-norm.
* Taxonomy of OTUs
* Top 25 important taxa
* Sample composition and read counts