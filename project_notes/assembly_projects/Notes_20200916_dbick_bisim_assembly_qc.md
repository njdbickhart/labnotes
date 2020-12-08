# Bison Simmental analysis 
---
*9/16/2020*

These are my notes on quality control and publication of the Bison x Simmental manuscript

## Table of Contents


## Setting up analysis pipeline

I am going to queue up these assemblies for analysis in separate sets. The bison and simmental in separate folders.

#### Bison

The Bison's name is Woody.

> Ceres: /lustre/project/cattle_genome_assemblies/bison_x_simmental/bison_x_simmental_qc/bison_qc

```bash
module load miniconda/3.6
cp ~/python_toolchain/snakeMake/assemblyValidation/default.json ./

ls /lustre/project/cattle_genome_assemblies/bison_x_simmental/illumina_data/Woody/*.fastq.gz | perl -lane 'print "    \"$F[0]\",";'

vim default.json

sbatch --nodes=1 --mem=5000 --ntasks-per-node=2 -p priority -q msn -t 8-0 snakemake --cluster-config ~/python_toolchain/snakeMake/assemblyValidation/cluster.json --cluster "sbatch --nodes={cluster.nodes} --ntasks-per-node={cluster.ntasks-per-node} --mem={cluster.mem} --partition={cluster.partition} -q {cluster.qos} -o {cluster.stdout} -t 8-0" -p --jobs 250 -s ~/python_toolchain/snakeMake/assemblyValidation/assemblyValidation --use-conda
```

Venn diagram analysis:

```bash
conda activate /KEEP/rumen_longread_metagenome_assembly/meryl

sbatch -N 1 -n 4 --mem=30000 -p priority -q msn -t 2-0 --wrap="python ~/python_toolchain/sequenceData/merylVennUpset.py -m /lustre/project/rumen_longread_metagenome_assembly/binaries/meryl/build/bin/meryl -o simmental_manuscript_comp --d merqury/arsucd/arsucd.meryl -d merqury/contigs/contigs.meryl -d merqury/polished/polished.meryl -d mapped/meryl_db.meryl"

sbatch -N 1 -n 4 --mem=30000 -p priority -q msn -t 2-0 --wrap="python ~/python_toolchain/sequenceData/merylVennUpset.py -m /lustre/project/rumen_longread_metagenome_assembly/binaries/meryl/build/bin/meryl -o simmental_manuscript_asmcomp --d merqury/arsucd/arsucd.meryl -d merqury/polished/polished.meryl"
```

#### Simmental

And the Simmental's name is Hollary (cow)

> Ceres: /lustre/project/cattle_genome_assemblies/bison_x_simmental/bison_x_simmental_qc/simmental_qc

```bash
module load miniconda/3.6
cp ~/python_toolchain/snakeMake/assemblyValidation/default.json ./

ls /lustre/project/cattle_genome_assemblies/bison_x_simmental/illumina_data/Hollary/*.fastq.gz | perl -lane 'print "    \"$F[0]\",";'

vim default.json
sbatch --nodes=1 --mem=5000 --ntasks-per-node=2 -p priority -q msn -t 5-0 snakemake --cluster-config ~/python_toolchain/snakeMake/assemblyValidation/cluster.json --cluster "sbatch --nodes={cluster.nodes} --ntasks-per-node={cluster.ntasks-per-node} --mem={cluster.mem} --partition={cluster.partition} -q {cluster.qos} -o {cluster.stdout} -t 5-0" -p --jobs 250 -s ~/python_toolchain/snakeMake/assemblyValidation/assemblyValidation --use-conda
```

#### Pulling missing buscos

```bash
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %data; while(<IN>){chomp; if($_ =~ /^#/){next;} $data{$_} = 1;} close IN; open(IN, "gunzip -c $ARGV[1] |"); while($h = <IN>){chomp $h; $s = <IN>; chomp $s; $h =~ s/>//; if(exists($data{$h})){$s =~ s/(.{1,60})/$1\n/g; print ">$h\n$s\n";}} close IN;' btemp_arsucd/run_mammalia_odb10/missing_busco_list.tsv busco_downloads/lineages/mammalia_odb10/refseq_db.faa.gz > btemp_arsucd/missing_genes.faa
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %data; while(<IN>){chomp; if($_ =~ /^#/){next;} $data{$_} = 1;} close IN; open(IN, "gunzip -c $ARGV[1] |"); while($h = <IN>){chomp $h; $s = <IN>; chomp $s; $h =~ s/>//; if(exists($data{$h})){$s =~ s/(.{1,60})/$1\n/g; print ">$h\n$s\n";}} close IN;' btemp_contigs/run_mammalia_odb10/missing_busco_list.tsv busco_downloads/lineages/mammalia_odb10/refseq_db.faa.gz > btemp_contigs/missing_genes.faa
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %data; while(<IN>){chomp; if($_ =~ /^#/){next;} $data{$_} = 1;} close IN; open(IN, "gunzip -c $ARGV[1] |"); while($h = <IN>){chomp $h; $s = <IN>; chomp $s; $h =~ s/>//; if(exists($data{$h})){$s =~ s/(.{1,60})/$1\n/g; print ">$h\n$s\n";}} close IN;' btemp_polished/run_mammalia_odb10/missing_busco_list.tsv busco_downloads/lineages/mammalia_odb10/refseq_db.faa.gz > btemp_polished/missing_genes.faa

for i in arsucd contigs polished; do echo $i; perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %data; while(<IN>){chomp; if($_ =~ /^#/){next;} $data{$_} = 1;} close IN; open(IN, "< $ARGV[1]"); while($h = <IN>){chomp $h; if($h =~ /^#/){next;} @s = split(/\t/, $h); if(exists($data{$s[0]})){print "$s[0]\t$s[1]\t$s[2]\n";}} close IN;' btemp_${i}/run_mammalia_odb10/missing_busco_list.tsv busco_downloads/lineages/mammalia_odb10/links_to_ODB10.txt > btemp_${i}/missing_genes_list.tab; done
```