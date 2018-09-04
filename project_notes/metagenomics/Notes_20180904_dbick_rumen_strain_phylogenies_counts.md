# Strain haplotyping
---
*9/4/2018*



> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/desman

```bash
for i in illumina illumina_rest illumina_sole pacbio pacbio_rest pacbio_sole; do echo $i; echo -e "experiment\tbin\tstrains" >> ${i}.combined.strain.counts; for j in $i/*/*.strain.count; do perl -e 'chomp(@ARGV); open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\s+/); print "$ARGV[0]\t$s[0]\t$s[1]\n";} close IN;' $i $j >> ${i}.combined.strain.counts; done; done
```


#### Testing strain counts


```R
library(ggplot2)
library(dplyr)

combined <- data.frame(experiment=character(), bin=character(), strains=numeric(), stringsAsFactors = FALSE); for(i in c("illumina_rest.combined.strain.counts", "illumina.combined.strain.counts", "pacbio_sole.combined.strain.counts", "illumina_sole.combined.strain.counts", "pacbio_rest.combined.strain.counts", "pacbio.combined.strain.counts")){print(i); temp <- read.delim(i, header=TRUE); temp$strains <- as.numeric(temp$strains); combined <- bind_rows(combined, temp);}
ggplot(combined, aes(x=experiment, y=strains, fill=experiment)) + geom_violin()
dev.copy2pdf(file="first_round_strain_counts.pdf", useDingbats=FALSE)
```


I'm a little concerned that the "illumina" and "pacbio" runs were initiated with a different version of the scripts and that may be a cause of discrepancy in the strain counts I see. Let me rerun them to be sure.

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/desman/illumina

```bash
mkdir old
mv illumina_* ./old
cp ../illumina_rest/run_desman_pipeline_illumina_rest.sh ./

for i in bin_lists/*.list; do name=`echo $i | perl -e '$h = <STDIN>; chomp($h); @bsegs = split(/[\._]/, $h); print "illumina_$bsegs[4]\_$bsegs[5]";'`; echo $name; bed=`basename $i | cut -d'.' -f1,2`; echo $bed; sbatch run_desman_pipeline_illumina.sh $name $bed; done
```

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/desman/pacbio

```bash
mkdir old
mv pacbio_* ./old/

cp ../pacbio_rest/run_desman_pipeline_pacbio_rest.sh ./
for i in bin_lists/*.list; do name=`echo $i | perl -e '$h = <STDIN>; chomp($h); @bsegs = split(/[\._]/, $h); print "pacbio_$bsegs[4]\_$bsegs[5]";'`; echo $name; bed=`basename $i | cut -d'.' -f1,2`; echo $bed; sbatch run_desman_pipeline_pacbio.sh $name $bed; done
```