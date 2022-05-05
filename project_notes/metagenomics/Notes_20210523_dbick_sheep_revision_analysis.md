# Sheep manuscript reviewer comments
---
*05/21/2021*

## Table of Contents
* [Kaiju analysis](#kaiju)

## SCAPP on CLR datasets

I predict that the reviewers will ask for summary statistics on the plasmid datasets similar to the mobile elements.

Let's try to run that preemptively to save time on reviewer revisions.

> Ceres: /lustre/project/rumen_longread_metagenome_assembly/assemblies/sheep/

```bash
module load python_3/3.6.6 miniconda/3.6

for i in clr1 clr2 clr3; do echo $i; sbatch -N 1 -n 2 --mem=85000 -p priority -q msn --wrap="python3 convertMetaFlyeGFA.py sheep_${i}/assembly_graph.gfa sheep_${i}/sheep_${i}.gfasta"; done

conda activate /KEEP/rumen_longread_metagenome_assembly/scapp/
module load samtools

for i in clr1 clr2 clr3; do sbatch -N 1 -n 70 --mem=320000 -p priority -q msn --wrap="scapp -g sheep_${i}/sheep_${i}.gfasta -o scapp_${i} -p 70 -r1 /project/rumen_longread_metagenome_assembly/sheep_poop/illumina_wgs/LIB101996_R1.fastq.gz -r2 /project/rumen_longread_metagenome_assembly/sheep_poop/illumina_wgs/LIB101996_R2.fastq.gz -sc false"; done
```

## Downsampling

One of the reviewers asked if we could downsample the HiFi reads to estimate the loss of fidelity through lower coverage. I will run the assemblies in pairs on the 8 Sequel II SMRTcells to generate these statistics. I am trying to equalize file sizes to reduce disparities in assembly size. Here's the pairing:

| Fastq1 | Fastq2 | Size 1 | Size 2| Total |
| :---   | :----  | :----- | :---- | :---- |
| m54337U_200227_201603.Q20.fasta | m54337U_200214_212611.Q20.fasta | 14G | 30G | 44G |
| m54337U_200222_075656.Q20.fasta | m54337U_200213_024013.Q20.fasta | 21G | 30G | 51G |
| m54337U_200203_184822.Q20.fasta | m54337U_200211_145859.Q20.fasta | 23G | 29G | 52G |
| m54337U_200220_195233.Q20.fasta | m54337U_200223_200916.Q20.fasta | 25G | 28G | 53G |

> Ceres: /lustre/project/rumen_longread_metagenome_assembly/assemblies/sheep

```bash
module load python_3/3.6.6 miniconda/3.6
source activate /KEEP/rumen_longread_metagenome_assembly/flye

sbatch --nodes=1 --mem=1000000 --ntasks-per-node=70 -p priority-mem -q msn-mem flye -g 1.0g --pacbio-hifi ../../sheep_poop/m54337U_200227_201603.Q20.fasta ../../sheep_poop/m54337U_200214_212611.Q20.fasta -o hifi_down_1 --meta -t 70 

sbatch --nodes=1 --mem=1000000 --ntasks-per-node=70 -p priority-mem -q msn-mem flye -g 1.0g --pacbio-hifi ../../sheep_poop/m54337U_200222_075656.Q20.fasta ../../sheep_poop/m54337U_200213_024013.Q20.fasta -o hifi_down_2 --meta -t 70

sbatch --nodes=1 --mem=1000000 --ntasks-per-node=70 -p priority-mem -q msn-mem flye -g 1.0g --pacbio-hifi ../../sheep_poop/m54337U_200203_184822.Q20.fasta ../../sheep_poop/m54337U_200211_145859.Q20.fasta -o hifi_down_3 --meta -t 70

sbatch --nodes=1 --mem=1000000 --ntasks-per-node=70 -p priority-mem -q msn-mem flye -g 1.0g --pacbio-hifi ../../sheep_poop/m54337U_200220_195233.Q20.fasta ../../sheep_poop/m54337U_200223_200916.Q20.fasta -o hifi_down_4 --meta -t 70

## Packaging to save important files before running new analysis
for i in hifi_down_?; do for j in assembly.fasta assembly_graph.gfa assembly_info.txt; do echo $i; echo $j; sbatch -N 1 -n 1 --mem=5000 -p priority -q msn --wrap="gzip $i/$j"; done; done
```

#### New strategy

At the co-authors' suggestions, I will be downsampling reads in 10% intervals to try to create about 9 subsampled assemblies for comparison.

```bash
# This should create the 9 fasta files
sheep=/lustre/project/rumen_longread_metagenome_assembly/sheep_poop/; sbatch -N 1 -n 3 --mem=25000 -p priority -q msn --wrap="perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/metagenomics_scripts/downsampleFastaForAssembly.pl $sheep/m54337U_200227_201603.Q20.fasta $sheep/m54337U_200214_212611.Q20.fasta $sheep/m54337U_200222_075656.Q20.fasta $sheep/m54337U_200213_024013.Q20.fasta $sheep/m54337U_200203_184822.Q20.fasta $sheep/m54337U_200211_145859.Q20.fasta $sheep/m54337U_200220_195233.Q20.fasta $sheep/m54337U_200223_200916.Q20.fasta"

# Getting sequence stats
for i in *.fasta; do sbatch -N 1 -n 1 --mem=5000 -p priority -q msn --wrap="samtools faidx $i" ;done

module load python_3/3.6.6 miniconda/3.6
source activate /KEEP/rumen_longread_metagenome_assembly/flye

# Queuing them up
for i in 10 20 30 40 50 60 70 80 90; do name="hifi_downsample_"$i; fasta="downsample."$i".fasta"; echo $name $fasta; sbatch --nodes=1 --mem=1000000 --ntasks-per-node=70 -p priority-mem -q msn-mem flye -g 1.0g --pacbio-hifi $fasta -o $name --meta -t 70;  done

mkdir links; for i in 10 20 30 40 50 60 70 80 90; do echo $i; ln -s /lustre/project/rumen_longread_metagenome_assembly/assemblies/sheep/hifi_downsample_${i}/assembly.fasta links/asm${i}.fa; done

Downsample  Read count  Actual percentage of reads  Gigabases   Ratio of total bases
0.9         16609388    0.899900065                 185.7920111 0.900155093
0.8         14764870    0.799963699                 165.1642142 0.800214216
0.7         12918151    0.699908083                 144.5019795 0.70010649
0.6         11072092    0.599888226                 123.8517703 0.600057027
0.5         9226397     0.499888091                 103.2008998 0.500004359
0.4         7380933     0.399900471                 82.55558086 0.399978589
0.3         5534992     0.299887007                 61.90968844 0.299950041
0.2         3688627     0.199850571                 41.26284906 0.199916904
0.1         1844693     0.099945847                 20.6337189  0.099969568


sbatch --nodes=1 --mem=1000000 --ntasks-per-node=70 -p priority-mem -q msn-mem flye -g 1.0g --pacbio-hifi $sheep/m54337U_200227_201603.Q20.fasta $sheep/m54337U_200214_212611.Q20.fasta $sheep/m54337U_200222_075656.Q20.fasta $sheep/m54337U_200213_024013.Q20.fasta $sheep/m54337U_200203_184822.Q20.fasta $sheep/m54337U_200211_145859.Q20.fasta $sheep/m54337U_200220_195233.Q20.fasta $sheep/m54337U_200223_200916.Q20.fasta -o asm00 --meta -t 70
```

> Ceres: /lustre/project/rumen_longread_metagenome_assembly/analysis/sheep 

```bash
module load python_3/3.6.6 miniconda/3.6 samtools

sbatch -N 1 -n 2 --mem=10000 -p priority -q msn snakemake -s ~/python_toolchain/snakeMake/hifiMAGManuscript/Downsample --cluster-config ~/python_toolchain/snakeMake/hifiMAGManuscript/cluster.json --cluster "sbatch -N 1 --ntasks-per-node={cluster.ntasks-per-node} --mem={cluster.mem} -p priority -q msn -o {cluster.stdout}" -p --use-conda --jobs 250 --verbose --latency-wait 40


# Creating a list of HQ contigs to intersect with the list of circular contigs
do echo -ne "$i\t"; perl -e '<>; <>; <>; while(<>){chomp; @s = split(/\s+/); if($s[-3] > 90 && $s[-2] < 10){print "$s[1]\n";}}' < $i > $i.hq.list; done

# Circular contig list
for i in 10 20 30 40 50 60 70 80 90; do file=/lustre/project/rumen_longread_metagenome_assembly/assemblies/sheep/hifi_downsample_${i}/assembly_info.txt; output=downsample_${i}.circ.list; echo $file; perl -lane 'if($F[0] =~ /#/){next;}elsif($F[3] eq "Y"){print $F[0];}' < $file > $output; done

# HQ circular contigs
for i in 10 20 30 40 50 60 70 80 90; do file=tables/asm${i}.contigs.checkm.txt.hq.list; output=downsample_${i}.circ.list; echo $file; perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl $file $output; done

# overwrote the protein faa files!
module load prodigalorffinder diamond/2.0.6; for i in asm10 asm20 asm30 asm40 asm50 asm60 asm70 asm80 asm90; do echo $i; sbatch -N 1 -n 30 --mem=60000 -p priority -q msn --wrap="DAS_Tool --outputbasename binning/DASTool/${i}.full --bins binning/bin3c/$i/bin3c.full.clusters.tab --labels bin3c --contigs assembly/${i}.fa --search_engine diamond --write_bin_evals 1 --threads 30 --debug"; done

# I need to prep bed files for my MAGPhase pipeline
for i in asm00 asm10 asm20 asm30 asm40 asm50 asm60 asm70 asm80 asm90; do echo $i; perl -e 'chomp(@ARGV); %keep; open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); $keep{$s[0]} = 1;} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $keep{$s[0]} = 1;} close IN; open(IN, "< binning/DASTool/$ARGV[2].full_proteins.faa");  while(<IN>){ unless($_ =~ />/){next;} chomp; @s = split(/\s+/); $s[0] =~ s/>//; if(exists($keep{$s[0]})){ $name = $s[0]; @d = split(/_/, $s[0]); $ctg = "$d[0]\_$d[1]"; print "$ctg\t$s[2]\t$s[4]\t$name\n";}} close IN;' binning/DASTool/${i}.full_proteins.faa.archaea.scg binning/DASTool/${i}.full_proteins.faa.bacteria.scg $i > binning/DASTool/${i}.full_proteins.scg.bed; done 


# Printing out only MAG beds
for i in asm00 asm10 asm20 asm30 asm40 asm50 asm60 asm70 asm80 asm90; do mkdir $i; rm -r $i/beds; mkdir $i/beds; echo $i; perl -e 'chomp(@ARGV); %data; open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); $s[1] =~ s/\./_/g; $data{$s[0]} = $s[1];} close IN;  open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); if(exists($data{$s[0]})){open(OUT, ">> $ARGV[2]/$data{$s[0]}.bed"); print {OUT} join("\t", @s) . "\n"; close OUT;}} close IN;' binning/DASTool/${i}.full_cluster_attribution.tsv binning/DASTool/${i}.full_proteins.scg.bed $i/beds ; done

# Printing out only MAG fastas
for i in asm00 asm10 asm20 asm30 asm40 asm50 asm60 asm70 asm80 asm90; do mkdir $i/mags; echo $i; perl -e 'chomp(@ARGV); %data; open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); $s[1] =~ s/\./_/g; $data{$s[0]} = $s[1];} close IN; open(IN, "< $ARGV[1].fai"); while(<IN>){chomp; @s = split(/\t/); if(exists($data{$s[0]})){system("samtools faidx $ARGV[1] $s[0] >> $ARGV[2]/$data{$s[0]}.fa"); }} close IN;' binning/DASTool/${i}.full_cluster_attribution.tsv assembly/${i}.fa $i/mags; done

for i in 00 10 20 30 40 50 60 70 80 90; do echo $i | perl -lane 'print "{\n\"assembly\" : \"/lustre/project/rumen_longread_metagenome_assembly/analysis/sheep/assembly/asm$F[0].fa\",\n\"reads\" : \"/lustre/project/rumen_longread_metagenome_assembly/assemblies/sheep/downsample.$F[0].fasta\"\n}";' > asm${i}/default.json; done

conda activate /KEEP/rumen_longread_metagenome_assembly/desman/
module load minimap2 samtools


for i in asm00 asm10 asm20 asm30 asm40 asm50 asm60 asm70 asm80 asm90; do echo $i; cd $i; rm mags/*.fai; sbatch -N 1 -n 1 --mem=10000 -p priority -q msn snakemake -s ~/python_toolchain/snakeMake/hifiMAGManuscript/MAGPhase --cluster-config ../cluster.json --cluster "sbatch -N {cluster.nodes} -n {cluster.ntasks-per-node} -p priority -q msn --mem={cluster.mem} -J {cluster.jobname}" --jobs 20; cd ..; done

for i in asm10 asm20 asm30 asm40 asm50 asm60 asm70 asm80 asm90; do echo $i; cd $i; cat consolidate/*/*.short > consolidate/consolidated.short; cat consolidate/*/*.long > consolidate/consolidated.long; cd ..; done

for i in asm00 asm10 asm20 asm30 asm40 asm50 asm60 asm70 asm80 asm90; do echo $i; perl -e 'chomp(@ARGV); $asm = $ARGV[0]; open(IN, "< binning/DASTool/$asm.full_DASTool_summary.txt"); %keep; <IN>; while(<IN>){chomp; @s = split(/\t/); $s[0] =~ s/\./_/; if($s[-2] >= 90 && $s[-1] < 5){$keep{$s[0]} = 1;}} close IN; open(IN, "< $asm/consolidate/consolidated.short"); while(<IN>){chomp; @s = split(/\t/); $s[0] =~ s/\.fa//g; if(exists($keep{$s[0]}) && $s[1] == 0){print "$s[0]\n";}} close IN;' $i | wc -l ; done
asm10	5
asm20	10
asm30	17
asm40	29
asm50	45
asm60	44
asm70	60
asm80	66
asm90	74


for j in asm00 asm10 asm20 asm30 asm40 asm50 asm60 asm70 asm80 asm90; do echo $j; for i in `perl -lane 'if($F[-2] > 90 && $F[-1] < 10){$F[0] =~ s/\./_/; print "$F[0]";}' < binning/DASTool/$j.full_DASTool_summary.txt`; do if [ -f $j/phase/${i}/${i}.fa.NO_SNPS_FOUND ]; then echo $i; fi; done | wc -l; done

for j in asm00 asm10 asm20 asm30 asm40 asm50 asm60 asm70 asm80 asm90; do echo $j; perl -lane 'if($F[-2] >= 90 && $F[-1] < 5){print "$F[0]";}' < binning/DASTool/$j.full_DASTool_summary.txt |wc -l; done

# It's surprisingly low! Let's test this pipeline out on the hifi assembly to make sure that things are Kosher!
mkdir hifi
mkdir hifi/mags
mkdir hifi/beds

cp /lustre/project/forage_assemblies/sheep_project/complete_flye/gtdbtk_bins/flye4.contigs/*.fna hifi/mags/
for i in hifi/mags/*.fna; do base=`basename $i | cut -d'.' -f1 | cut -d'_' -f2`; echo $base; cp /lustre/project/forage_assemblies/sheep_project/complete_flye/desman/bed_lists/flye4/bin3c.${base}.scg.bed ./hifi/beds/bin3c_${base}.bed; done

cd hifi/
sbatch -N 1 -n 1 --mem=10000 -p priority -q msn snakemake -s ~/python_toolchain/snakeMake/hifiMAGManuscript/MAGPhase --cluster "sbatch -N 1 -n 4 -p priority -q msn --mem=35000" --jobs 75
cd ..

perl -e 'chomp(@ARGV); $asm = $ARGV[0]; open(IN, "< /lustre/project/forage_assemblies/sheep_project/complete_flye/b3c_flye4_dastool/flye4.das_DASTool_summary.txt"); %keep; <IN>; while(<IN>){chomp; @s = split(/\t/); $s[0] =~ s/\./_/; if($s[-2] >= 90 && $s[-1] <= 10){$keep{$s[0]} = 1;}} close IN; open(IN, "< $asm/consolidate/consolidated.short"); while(<IN>){chomp; @s = split(/\t/); $s[0] =~ s/\.fna//g; if(exists($keep{$s[0]}) && $s[1] == 0){print "$s[0]\n";}} close IN;' hifi | wc -l
254   <- I used the wrong SCG bed files! AGH!
```

#### Summary statistics

```bash
python3 ~/python_toolchain/sequenceData/bamStats.py -f /lustre/project/forage_assemblies/sheep_project/complete_flye/flye4.contigs.fasta.ccs.bam
#Genome size: {'3,426,390,363'} and avg read len: {'11679'}. Cropping similar columns
File    TotalReads      MappedReads     MapPerc RawXCov MapXCov AvgMapChrXCov   AvgChrMapPerc
flye4.contigs.fasta.ccs.bam     36,362,602      35,806,885      0.985   123.944 122.049 691.551 0.975

for i in asm??; do echo $i; python3 ~/python_toolchain/sequenceData/bamStats.py -f $i/mapping/assembly.filtered.bam; done
asm10
#Genome size: {'885,697,461'} and avg read len: {'10736'}. Cropping similar columns
File    TotalReads      MappedReads     MapPerc RawXCov MapXCov AvgMapChrXCov   AvgChrMapPerc
assembly.filtered.bam   14,777,238      13,516,912      0.915   179.123 163.846 74.157  0.865
asm20
#Genome size: {'1,396,849,901'} and avg read len: {'9733'}. Cropping similar columns
File    TotalReads      MappedReads     MapPerc RawXCov MapXCov AvgMapChrXCov   AvgChrMapPerc
assembly.filtered.bam   15,433,548      14,463,069      0.937   107.538 100.776 43.563  0.823
asm30
#Genome size: {'1,778,316,973'} and avg read len: {'10681'}. Cropping similar columns
File    TotalReads      MappedReads     MapPerc RawXCov MapXCov AvgMapChrXCov   AvgChrMapPerc
assembly.filtered.bam   15,882,816      14,852,627      0.935   95.396  89.208  37.522  0.796
asm40
#Genome size: {'2,105,540,837'} and avg read len: {'10066'}. Cropping similar columns
File    TotalReads      MappedReads     MapPerc RawXCov MapXCov AvgMapChrXCov   AvgChrMapPerc
assembly.filtered.bam   16,094,347      15,299,866      0.951   76.943  73.144  30.073  0.778
asm50
#Genome size: {'2,395,203,741'} and avg read len: {'10401'}. Cropping similar columns
File    TotalReads      MappedReads     MapPerc RawXCov MapXCov AvgMapChrXCov   AvgChrMapPerc
assembly.filtered.bam   16,271,516      15,489,458      0.952   70.658  67.262  27.309  0.760
asm60
#Genome size: {'2,652,034,387'} and avg read len: {'11093'}. Cropping similar columns
File    TotalReads      MappedReads     MapPerc RawXCov MapXCov AvgMapChrXCov   AvgChrMapPerc
assembly.filtered.bam   16,322,995      15,715,282      0.963   68.276  65.734  26.942  0.748
asm70
#Genome size: {'2,728,802,791'} and avg read len: {'10602'}. Cropping similar columns
File    TotalReads      MappedReads     MapPerc RawXCov MapXCov AvgMapChrXCov   AvgChrMapPerc
assembly.filtered.bam   15,679,267      15,097,021      0.963   60.917  58.655  23.304  0.729
asm80
#Genome size: {'2,863,821,278'} and avg read len: {'10762'}. Cropping similar columns
File    TotalReads      MappedReads     MapPerc RawXCov MapXCov AvgMapChrXCov   AvgChrMapPerc
assembly.filtered.bam   15,440,668      14,748,970      0.955   58.025  55.425  23.095  0.705
asm90
#Genome size: {'3,018,971,575'} and avg read len: {'10778'}. Cropping similar columns
File    TotalReads      MappedReads     MapPerc RawXCov MapXCov AvgMapChrXCov   AvgChrMapPerc
assembly.filtered.bam   15,552,491      14,841,793      0.954   55.524  52.987  21.849  0.700
```


## Effects of illumina polishing

The question is whether or not Illumina polishing actually hurts HiFi assembly by mashing variant sites together. My thoughts are that we can address this by polishing the HiFi assembly and then assessing INDELs in ORFs using IDEEL.

> Ceres: /lustre/project/forage_assemblies/sheep_project/complete_flye

```bash
module load python_3/3.6.6

sbatch -N 1 -n 2 --mem=15000 -p priority -q msn --wrap="python3 ~/python_toolchain/sequenceData/slurmPilonCorrectionPipeline.py -f mapping/flye4.contigs/merged.bam -g flye4.contigs.fasta -o flye4_pilon1 -p priority -q msn -e 30000 -m"

for i in flye4_pilon1/*.fasta; do echo $i; perl -ne 'if($_ =~ /^>/){$_ =~ s/_pilon//;} print $_;' < $i >> flye4.pilon.fasta; done

## Initial comparisons
python3 ~/python_toolchain/utils/tabFileLeftJoinTable.py -f flye4.contigs.fasta.fai -f flye4.pilon.fasta.fai -c 0 -m 1 -o flye4.comparison.contiglengths.tab

perl -lane 'print $F[1] - $F[2];' < flye4.comparison.contiglengths.tab | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   60050
Sum:    145845
Minimum -59
Maximum 1236
Average 2.428726
Median  0
Standard Deviation      12.918394
Mode(Highest Distributed Value) 0

perl -lane 'if($F[1] - $F[2] == 0){print $F[1] - $F[2];}' < flye4.comparison.contiglengths.tab | wc -l
39914 (out of 60,050; or 66.5%).

git clone https://github.com/mw55309/ideel.git
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz

conda activate /KEEP/rumen_longread_metagenome_assembly/ideel

sbatch -N 1 -n 70 -p priority -q msn --mem=89000 --wrap="diamond makedb -p 70 --db uniprot_trembl.fasta.db --in uniprot_trembl.fasta.gz --taxonnodes /lustre/project/bostauruscnv/novel_seq/cddr_pipeline/taxonomy/nodes.dmp --taxonnames /lustre/project/bostauruscnv/novel_seq/cddr_pipeline/taxonomy/names.dmp"

cp flye4.contigs.fasta genomes/unpolished.fa
cp flye4.pilon.fasta genomes/pilon.fa

sbatch -N 1 -n 2 -p priority -q msn --mem=15000 --wrap='snakemake -s ideel/Snakefile --cluster "sbatch -N 1 -n 70 --mem=300000 -p priority -q msn" --jobs 10'

for i in pilon unpolished; do echo $i; sbatch -N 1 -n 2 --mem=36000 -p priority -q msn --wrap="Rscript ideel/scripts/hist.R lengths/$i.data hists/$i.png"; done


# Getting the count of pilon detections and corrections
grep 'Found' flye4_pilon1/outLog/*.out | perl -e 'print "SNPs\tAMB\tINS\tINSBP\tDEL\tDELBP\n"; while(<>){chomp; ($snp, $amb, $ins, $insbp, $del, $delbp) = $_ =~ m/Found (\d+) snps\; (\d+) ambiguous bases\; corrected (\d+) small insertions totaling (\d+) bases, (\d+) small deletions totaling (\d+) bases/; if($snp eq ""){next;}else{print "$snp\t$amb\t$ins\t$insbp\t$del\t$delbp\n";}}' > flye4.pilon.corrections.tab
```

#### Quick Q-Q plot to determine deviations in quartiles

I want to compare the distributions to see if there are significant deviations

```R
setwd("C:/SharedFolders/metagenomics/tim_sheep/ideel_check/")

pilon <- read.delim("pilon.data", header = FALSE)
unpolished <- read.delim("unpolished.data", header=FALSE)

pilon.vals <- pilon$V1 / pilon$V2
unpolished.vals <- unpolished$V1 / unpolished$V2

t.test(pilon.vals, unpolished.vals)

	Welch Two Sample t-test

data:  pilon.vals and unpolished.vals
t = 35.878, df = 5975261, p-value < 2.2e-16
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 0.008777689 0.009792129
sample estimates:
mean of x mean of y 
0.9436847 0.9343998

qq <- qqplot(x = pilon.vals, y=unpolished.vals, plot.it=FALSE)

library(ggplot2)
qq <- as.data.frame(qq)
ggplot(qq, aes(x=x, y=y)) + geom_point() + geom_abline(intercept=0, slope=1) + coord_fixed(ratio =1, xlim=c(0, 1.5), ylim=c(0,1.5)) + xlab("Pilon Corrected ORF Ratios") + ylab("Unpolished ORF Ratios") + theme_bw()

dev.copy2pdf(file="pilon_unpolished_qqplot.pdf", useDingbats=FALSE)

```

### IDEEL plot

> C:/SharedFolders/metagenomics/tim_sheep/ideel_check

```R
library(dplyr)
library(ggplot2)

pilon <- read.delim("pilon.data", header=FALSE)
unpolished <- read.delim("unpolished.data", header=FALSE)
pilon <- pilon %>% mutate(Ratio = V1 / V2, Name = c("Polished"))
unpolished <- unpolished %>% mutate(Ratio = V1 / V2, Name = c("Unpolished"))

pilon <- pilon[,c("Ratio", "Name")]
unpolished <- unpolished[,c("Ratio", "Name")]

full <- bind_rows(pilon, unpolished)
full$Name <- as.factor(full$Name)

ggplot(full, aes(x=Ratio, fill=Name, color=Name)) + geom_histogram(position="dodge", alpha=0.5) + scale_color_brewer(palette="Dark2") + scale_fill_brewer(palette="Dark2") + theme_bw() + xlim(0, 1.5)
```

### Quast analysis

> Ceres: /lustre/project/forage_assemblies/sheep_project/complete_flye/

```bash
module load python_3/3.6.6 miniconda/3.6
conda activate /KEEP/rumen_longread_metagenome_assembly/quast

sbatch -N 1 -n 70 --mem=300000 -p priority -q msn --wrap="metaquast flye4.contigs.fasta clr1.contigs.fasta clr2.contigs.fasta clr3_rerun/assembly/clr3.contigs.fa --max-ref-number 150 -o metaquast --no-check -t 70"
```

<a name="kaiju"></a>
### Centrifuge analysis

> Ceres: /lustre/project/rumen_longread_metagenome_assembly/analysis/sheep

```bash
sbatch -N 1 -n 70 --mem=900000 -p priority-mem -q msn-mem --wrap='/lustre/project/rumen_longread_metagenome_assembly/binaries/centrifuge/centrifuge -x refseq_abv_hum_mus_bos -1 /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/illumina_wgs/LIB101996_R1.fastq.gz -2 /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/illumina_wgs/LIB101996_R2.fastq.gz -p 70 --report-file sheep_illumina.centrifuge.report > sheep_illumina.centrifuge.out'

sbatch -N 1 -n 70 --mem=900000 -p priority-mem -q msn-mem --wrap='/lustre/project/rumen_longread_metagenome_assembly/binaries/centrifuge/centrifuge -x refseq_abv_hum_mus_bos -U /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/sheep_poop_sequel_combined_CCS.fasta.gz -f  -p 70 --report-file sheep_hifi.centrifuge.report > sheep_hifi.centrifuge.out'

for i in m54337U_200203_184822 m54337U_200211_145859 m54337U_200213_024013 m54337U_200214_212611 m54337U_200220_195233 m54337U_200222_075656 m54337U_200223_200916 m54337U_200227_201603; do file=/lustre/project/rumen_longread_metagenome_assembly/sheep_poop/${i}.fastq.gz; echo $file; sbatch -N 1 -n 70 --mem=900000 -p priority-mem -q msn-mem --wrap="/lustre/project/rumen_longread_metagenome_assembly/binaries/centrifuge/centrifuge -x refseq_abv_hum_mus_bos -U $file -p 70 --report-file sheep_${i}.centrifuge.report > sheep_${i}.centrifuge.out"; done

for i in *.report; do echo $i; perl -e '<>; while(<>){chomp; @s = split(/\t/); print "$s[0]\n";}' < $i > $i.list; done
cat sheep_m54337U*.list | sort | uniq > sheep_subreads.centrifuge.report.list

perl concatenate_reports.pl > sheep_subreads.centrifuge.report

for i in sheep_illumina.centrifuge.report sheep_hifi.centrifuge.report sheep_subreads.centrifuge.report; do echo $i; perl -e ' <>; while(<>){chomp; @s = split(/\t/); print "$s[0]\t$s[1]\t$s[5]\n";}' < $i > $i.prekrona; done
for i in sheep_illumina.centrifuge.report sheep_hifi.centrifuge.report sheep_subreads.centrifuge.report; do echo $i; /lustre/project/rumen_longread_metagenome_assembly/binaries/Krona/bin/ktImportTaxonomy -o $i.krona.html $i.prekrona; done


## Let's try this with Kaiju
module load kaiju/1.7.2

sbatch -N 1 -n 1 --mem=5000 -p priority -q msn --wrap="wget http://kaiju.binf.ku.dk/database/kaiju_db_nr_euk_2021-02-24.tgz"
sbatch -N 1 -n 1 --mem=5000 -p priority -q msn --wrap="tar -xvf kaiju_db_nr_euk_2021-02-24.tgz"

sbatch -N 1 -n 70 --mem=900000 -p priority-mem -q msn-mem --wrap='kaiju -t nodes.dmp -f kaiju_db_nr_euk.fmi -i /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/illumina_wgs/LIB101996_R1.fastq.gz -j /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/illumina_wgs/LIB101996_R2.fastq.gz -z 70 -o sheep_illumina_kaiju.out'

sbatch -N 1 -n 70 --mem=900000 -p priority-mem -q msn-mem --wrap='kaiju -t nodes.dmp -f kaiju_db_nr_euk.fmi -i pacbio_subreads.fixed.fastq.gz -z 70 -o sheep_subread_kaiju.out'

sbatch -N 1 -n 70 --mem=900000 -p priority-mem -q msn-mem --wrap='kaiju -t nodes.dmp -f kaiju_db_nr_euk.fmi -i /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/sheep_poop_sequel_combined_CCS.fasta.gz -z 70 -o sheep_hifi_kaiju.out'

sbatch -N 1 -n 2 -p priority -q msn --mem=6000 --wrap='ls /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/m54337U*.fastq.gz | grep -v clr | xargs -I {} gunzip -c {} | gzip > pacbio_subreads.fastq.gz'

for i in sheep_hifi_kaiju sheep_illumina_kaiju; do echo $i; sbatch -N 1 -n 2 -p priority -q msn --wrap="kaiju2krona -t nodes.dmp -n names.dmp -i $i.out -o $i.krona"; done

sbatch -N 1 -n 2 --mem=6000 -p priority -q msn --wrap="kaiju2krona -t nodes.dmp -n names.dmp -i sheep_subread_kaiju.out -o sheep_subread_kaiju.krona"

/lustre/project/rumen_longread_metagenome_assembly/binaries/Krona/bin/ktImportText -o hifi_illumina_reads.krona.html sheep_illumina_kaiju.krona,illumina sheep_hifi_kaiju.krona,hifi

/lustre/project/rumen_longread_metagenome_assembly/binaries/Krona/bin/ktImportText -o hifi_illumina_sub_reads.krona.html sheep_illumina_kaiju.krona,illumina sheep_hifi_kaiju.krona,hifi sheep_subread_kaiju.krona,subreads


# Checking read lengths
grep 'Proteobacteria' names.dmp | perl -lane 'print "$F[0]\tProteobacteria";' > taxa_groups.list
grep 'Firmicutes' names.dmp | perl -lane 'print "$F[0]\tFirmicutes";' >> taxa_groups.list
grep 'Bacteroidetes' names.dmp | perl -lane 'print "$F[0]\tBacteroidetes";' >> taxa_groups.list

sbatch -N 1 -n 70 --mem=900000 -p priority-mem -q msn-mem --wrap='kaiju -t nodes.dmp -f kaiju_db_nr_euk.fmi -i /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/sheep_poop_sequel_combined_CCS.fasta.gz -v -z 70 -o sheep_hifi_kaiju.verbose.out'

sbatch -N 1 -n 70 --mem=900000 -p priority-mem -q msn-mem --wrap='kaiju -t nodes.dmp -f kaiju_db_nr_euk.fmi -i /pacbio_subreads.fixed.fastq.gz -v -z 70 -o sheep_subread_kaiju.verbose.out'


perl -e 'chomp(@ARGV); %taxa; open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); $taxa{$s[0]} = $s[1];} close IN; print "Length\tTaxa\n"; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); if(exists($taxa{$s[2]})){print "$s[3]\t$taxa{$s[2]}\n";}} close IN;' taxa_groups.list sheep_hifi_kaiju.verbose.out > sheep_hifi_kaiju.verbose.taxa.lens
```

#### Kmer venn analysis

```bash
PATH=$PATH:/software/7/apps/merqury/1.1/:/software/7/apps/meryl/1.0/Linux-amd64/bin/
MERQURY=/software/7/apps/merqury/1.1

sbatch -N 1 -n 30 --mem=50000 -p priority -q msn --wrap="meryl count k=21 ./pacbio_subreads.fixed.fastq.gz  threads=30 memory=48 output subreads.meryl"

sbatch -N 1 -n 30 --mem=50000 -p priority -q msn --wrap="meryl count k=21 lustre/project/rumen_longread_metagenome_assembly/sheep_poop/illumina_wgs/LIB101996_R1.fastq.gz  threads=30 memory=48 output illumina.meryl"

sbatch -N 1 -n 30 --mem=50000 -p priority -q msn --wrap="meryl count k=21 /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/sheep_poop_sequel_combined_CCS.fasta.gz  threads=30 memory=48 output hifi.meryl"

conda activate /KEEP/rumen_longread_metagenome_assembly/meryl
sbatch -N 1 -n 4 --mem=50000 -p priority -q msn -t 2-0 --wrap="python ~/python_toolchain/sequenceData/merylVennUpset.py -m /lustre/project/rumen_longread_metagenome_assembly/binaries/meryl/build/bin/meryl -o sheep_reads_venn_comp -d subreads.meryl -d hifi.meryl -d illumina.meryl"
```

### HiFi read mapping to lineage-resolved assemblies

> Ceres: /lustre/project/forage_assemblies/sheep_project/complete_flye/clr3_rerun

```bash
module load samtools

perl -lane 'if($F[-2] >= 90 && $F[-1] <= 10 && $F[1] == 0){print "$F[0]";}' < magphase/consolidated/clr3.consolidated.filt.short.tab > magphase/consolidated/clr3.consolidated.lineagemags.list
head -n 1 clr3_master_table_03_2021.tab | perl -lane 'for($x = 0; $x < scalar(@F); $x++){print "$x\t$F[$x]";}'

for i in clr1 clr2 clr3 flye4; do echo $i; perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %not; while(<IN>){chomp; @s = split(/\t/); $not{$s[1]} = 1;} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); if(exists($not{$s[1]})){next;}else{print "$s[0]\n";}} close IN;' binning/DASTool/$i.contigs.full_cluster_attribution.tsv binning/bin3c/$i.contigs/bin3c.full.clusters.tab > $i.lowqual.contigs.list; done

perl tabulate_read_mapping.pl 14 8 21 22 clr3_master_table_03_2021.tab assembly/clr3.contigs.fa.ccs.check.bam magphase/consolidated/clr3.consolidated.lineagemags.list clr3.lowqual.contigs.list clr3_hifi_read_categories.tab

for i in clr1 clr2 flye4; do echo $i; perl -lane 'if($F[-2] >= 90 && $F[-1] <= 10 && $F[1] == 0){print "$F[0]";}' < ../magphase/consolidated/$i.consolidated.short.tab > ../magphase/consolidated/$i.consolidated.lineagemags.list; done

for i in clr1 clr2; do echo $i; perl tabulate_read_mapping.pl 13 7 20 21 ../${i}_master_table_12_2020.tab ../$i.contigs.fasta.ccs.bam ../magphase/consolidated/$i.consolidated.lineagemags.list ${i}.lowqual.contigs.list ${i}_hifi_read_categories.tab; done

perl tabulate_read_mapping.pl 13 7 19 20 ../flye4_master_table_12_2020.tab ../flye4.contigs.fasta.ccs.bam flye4.lineage_resolved_HQ_bowers.list flye4.lowqual.contigs.list flye4.bowers.hqmags.list flye4_hifi_read_categories.tab

perl tabulate_read_mapping.revised.pl 13 7 19 20 ../flye4_master_table_12_2020.tab ../flye4.contigs.fasta.ccs.bam ../magphase/consolidated/flye4.consolidated.lineagemags.list flye4.lowqual.contigs.list flye4_hifi_read_categories.revised.tab

perl -e 'chomp(@ARGV); %perc; @rows = ("Kingdom", "MAG"); @files; foreach $f (@ARGV){ @name = split(/_/, $f); push(@files, $name[0]); open(IN, "< $f"); <IN>; while(<IN>){chomp; @s = split(/\t/); $perc{$s[0]}->{$s[1]}->{$name[0]} = $s[3]; } close IN;} print join("\t", @rows) . "\t" . join("\t", @files); print "\n"; foreach $t (keys(%perc)){foreach $m (keys(%{$perc{$t}})){print "$t\t$m"; foreach $d (@files){print "\t" . $perc{$t}->{$m}->{$d}} print "\n"; }}' flye4_hifi_read_categories.tab clr1_hifi_read_categories.tab clr2_hifi_read_categories.tab clr3_hifi_read_categories.tab > combined_hifi_read_categories.tab
```

```R
data <- read.delim("flye4_mapping_categories.txt", header=TRUE)
summary(data)
data$Kingdom <- as.factor(data$Kingdom)
data$MAG <- as.factor(data$MAG)

data <- data %>% mutate(ASM = "HiFi")
sum(data$ReadNum)

data <- data %>% mutate(Perc = ReadNum / 36362431 * 100)

data.trim <- data %>% group_by(MAG) %>% summarize(count=sum(ReadNum)) %>% summarise(MAG = MAG, Perc = count / 36362431 * 100)
data.trim <- data.trim %>% mutate(ASM = "HiFi")

ggplot(data.trim, aes(x = ASM, y=Perc, fill=MAG)) + geom_col() + geom_text(aes(label = paste0(Perc, "%")), position = position_stack(vjust = 0.5)) + scale_fill_brewer(palette = "Set2") + theme_bw() + xlab(NULL) + ylab("Percentage Reads Mapped")
dev.copy2pdf(file="hifi_mag_mapping_stats.pdf", useDingbats=FALSE)

ggplot(data, aes(x=Kingdom, y=Perc, fill=MAG)) + geom_bar(stat="identity") + scale_fill_brewer(palette = "Set2") + theme_bw()
dev.copy2pdf(file="hifi_mag_mapping_by_kingdom.pdf", useDingbats=FALSE)
```


### TRNA and rRNA identification for HQ MAG designation

> Ceres: /lustre/project/rumen_longread_metagenome_assembly/analysis/sheep

```bash
module load infernal/1.1.2
mkdir logs/trnascan
mkdir trnas

for i in asm10; do mkdir trnas/$i; for j in $i/mags/*.fa; do echo $j; name=`basename $j | cut -d'.' -f1`; echo $name; sbatch trna_scanner.sh $j trnas/$i/$name.out; done; done

for i in asm20 asm30 asm40 asm50 asm60 asm70 asm80 asm90; do mkdir trnas/$i; for j in $i/mags/*.fa; do echo $j; name=`basename $j | cut -d'.' -f1`; echo $name; sbatch trna_scanner.sh $j trnas/$i/$name.out; done; done

## Making the clr and hifi bins available
for i in clr1 clr2 clr3 flye4; do echo $i; mkdir $i; mkdir $i/mags; for j in /lustre/project/forage_assemblies/sheep_project/complete_flye/clr3_rerun/gtdbtk_bins/$i.contigs/*.fna; do name=`basename $j | cut -d'.' -f1`; ln -s $j $i/mags/$name.fa; done; done

for i in asm00 clr1 clr2 clr3 flye4; do mkdir trnas/$i; for j in $i/mags/*.fa; do echo $j; name=`basename $j | cut -d'.' -f1`; echo $name; sbatch trna_scanner.sh $j trnas/$i/$name.out; done; done


module load checkm/1.1.2
mkdir checkm
for i in asm10 asm20 asm30 asm40 asm50 asm60 asm70 asm80 asm90; do mkdir checkm/$i; sbatch -N 1 -n 30 --mem=100000 -p priority -q msn --wrap="checkm lineage_wf -f checkm/$i/$i.checkm.txt -t 30 --pplacer_threads 3 -x fa $i/mags checkm/$i"; done

for i in asm00 clr1 clr2 clr3 flye4; do mkdir checkm/$i; sbatch -N 1 -n 30 --mem=100000 -p priority -q msn --wrap="checkm lineage_wf -f checkm/$i/$i.checkm.txt -t 30 --pplacer_threads 3 -x fa $i/mags checkm/$i"; done


for i in asm10 asm20 asm30 asm40 asm50 asm60 asm70 asm80 asm90 asm00 clr1 clr2 clr3 flye4; do echo -ne "$i\t"; perl -ne '<>; <>; <>; $c = 0; while(<>){$_; $_ =~ s/^\s+//; @s = split(/\s+/); if($s[-3] >= 90 && $s[-2] <= 5){$c++;}} print "$c\n";' < checkm/$i/$i.checkm.txt; done
asm10   85
asm20   132
asm30   173
asm40   192
asm50   214
asm60   237
asm70   240
asm80   268
asm90   269
asm00   250
clr1    225
clr2    217
clr3    234
flye4   283

for i in clr1 clr2 clr3 flye4; do echo $i; ln -s /lustre/project/forage_assemblies/sheep_project/complete_flye/clr3_rerun/binning/DASTool/$i.contigs.full_DASTool_summary.txt binning/DASTool/$i.full_DASTool_summary.txt; done

for i in asm10 asm20 asm30 asm40 asm50 asm60 asm70 asm80 asm90 asm00 clr1 clr2 clr3 flye4; do echo -ne "$i\t"; perl -ne '<>; $c = 0; while(<>){$_; $_ =~ s/^\s+//; @s = split(/\s+/); if($s[-2] >= 90 && $s[-1] <= 5){$c++;}} print "$c\n";' < binning/DASTool/$i.full_DASTool_summary.txt; done


## TRNA discovery
for i in asm10 asm20 asm30 asm40 asm50 asm60 asm70 asm80 asm90 asm00 clr1 clr2 clr3 flye4; do echo $i; rm trnas/$i.trna.counts.tab; for j in trnas/$i/*.out; do name=`basename $j | cut -d'.' -f1`; echo -ne "$name\t" >> trnas/$i.trna.counts.tab; grep -v Undet $j | grep -v Pseudo | perl -ne '<>; <>; <>; while(<>){$t = $_; chomp $t; @s = split(/\t/, $t); if($s[8] > 25){print $_;}}' | python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f stdin -c 4 -d '\t' | grep -v Entry | wc -l >> trnas/$i.trna.counts.tab; done; done

## RRNA discovery
conda activate /KEEP/rumen_longread_metagenome_assembly/barrnap
mkdir rrna
for i in asm10 asm20 asm30 asm40 asm50 asm60 asm70 asm80 asm90 asm00 clr1 clr2 clr3 flye4; do echo $i; mkdir rrna/$i; for j in $i/mags/*.fa; do echo $j; name=`basename $j | cut -d'.' -f1`; echo $name; sbatch rrna_scanner.sh $j rrna/$i/$name.out; done; done

for i in asm10 asm20 asm30 asm40 asm50 asm60 asm70 asm80 asm90 asm00 clr1 clr2 clr3 flye4; do echo $i; rm rrna/$i.rrna.count; for j in rrna/$i/*.out; do name=`basename $j | cut -d'.' -f1`; echo -ne "$name\t" >> rrna/$i.rrna.count; perl -e '<>; while(<>){chomp; ($d) = $_ =~ m/Name=(.{6,10});/; print "$d\n";}' < $j | sort | uniq | wc -l >> rrna/$i.rrna.count; done; done

## get the full set of HQ MAGs
for i in asm10 asm20 asm30 asm40 asm50 asm60 asm70 asm80 asm90 asm00 clr1 clr2 clr3 flye4; do echo $i; perl -e '%data; chomp(@ARGV); open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); if($s[1] >= 3){$data{$s[0]} += 1;}} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); if($s[1] >= 18){$data{$s[0]} += 1;}} close IN; open(IN, "< $ARGV[2]"); <IN>; while(<IN>){$_ =~ s/^\s+//; @s = split(/\s+/); $s[0] =~ s/\./_/; if($s[-2] >= 90 && $s[-1] <= 5){$data{$s[0]} += 1;}} close IN; foreach $k (keys(%data)){if($data{$k} == 3){print "$k\n";}}' rrna/$i.rrna.count trnas/$i.trna.counts.tab binning/DASTool/$i.full_DASTool_summary.txt | wc -l; done
```

#### Making list of Bowers et al MAGs

> Ceres: /project/rumen_longread_metagenome_assembly/analysis/sheep

```bash
wd=/OLD/project/rumen_longread_metagenome_assembly/analysis/sheep
for i in asm10 asm20 asm30 asm40 asm50 asm60 asm70 asm80 asm90 asm00 clr1 clr2 clr3 flye4; do echo $i; perl -e '%data; chomp(@ARGV); open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); if($s[1] >= 3){$data{$s[0]} += 1;}} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); if($s[1] >= 18){$data{$s[0]} += 1;}} close IN; open(IN, "< $ARGV[2]"); <IN>; while(<IN>){$_ =~ s/^\s+//; @s = split(/\s+/); $s[0] =~ s/\./_/; if($s[-2] >= 90 && $s[-1] <= 5){$data{$s[0]} += 1;}} close IN; foreach $k (keys(%data)){if($data{$k} == 3){print "$k\n";}}' $wd/rrna/$i.rrna.count $wd/trnas/$i.trna.counts.tab $wd/binning/DASTool/$i.full_DASTool_summary.txt > $i.mimag.tab; done

# Just for flye4
for i in flye4; do echo $i; perl -e '%data; chomp(@ARGV); open(IN, "< $ARGV[0]"); %values; print "MAG\tRRNA\tTRNA\tCOMP\tCON\n"; while(<IN>){chomp; @s = split(/\t/); push(@{$values{$s[0]}}, $s[1]); if($s[1] >= 3){$data{$s[0]} += 1;}} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); push(@{$values{$s[0]}}, $s[1]); if($s[1] >= 18){$data{$s[0]} += 1;}} close IN; open(IN, "< $ARGV[2]"); <IN>; while(<IN>){$_ =~ s/^\s+//; @s = split(/\s+/); $s[0] =~ s/\./_/; push(@{$values{$s[0]}}, ($s[-2], $s[-1])); if($s[-2] >= 90 && $s[-1] <= 5){$data{$s[0]} += 1;}} close IN; foreach $k (keys(%data)){if($data{$k} == 3){print "$k\t" . join("\t", @{$values{$k}}) . "\n";}}' $wd/rrna/$i.rrna.count $wd/trnas/$i.trna.counts.tab $i.contigs.full_DASTool_summary.txt > $i.mimag.tab; done
flye4
```


```R
library(dplyr)
library(reshape)
library(ggplot2)

data <- read.delim("final_statistics.txt", header=TRUE)
data.melt <- melt(data, id=c("Downsample", "Gigabases"))
ggplot(data.melt, aes(x=Gigabases, y=value, fill=variable)) + geom_bar(stat="identity") + theme_bw() + scale_fill_brewer(palette="Dark2") + xlab("HiFi Read Subsets (Gbp)") + facet_wrap(~ variable, scales="free")
dev.copy2pdf(file="downsample_analysis_results.pdf", useDingbats=FALSE)
```

#### Misha's dedup bins

```bash
module load python_3/3.6.6 miniconda/3.6 usearch/11.0.667

perl -e 'use File::Basename; @files = `ls bins_dedup_2/*.fa`; chomp(@files); open(OUT, "> misha_dedup.list"); foreach $f (@files){@names = split(/\./, basename($f)); open(IN, "<$f"); while(<IN>){if($_ =~ /^>/){$_ =~ s/>//; chomp; print {OUT} "$_\t$names[0]\_$names[1]\n";}} close IN;} close OUT;'

conda activate /KEEP/rumen_longread_metagenome_assembly/das_tool
module load diamond/2.0.6; sbatch -N 1 -n 70 --mem=200000 -p priority -q msn --wrap="DAS_Tool --outputbasename dedup_dastool --bins misha_dedup.list --labels bin3cDedup --contigs /lustre/project/forage_assemblies/sheep_project/complete_flye/flye4.contigs.fasta --search_engine diamond --write_bin_evals 1 --create_plots 1 --threads 70"
```