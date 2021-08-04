# Rambouillet assembly analysis
---
*8/11/2020*

These are my notes on the QC and analysis of the Rambouillet genome.

## Table of Contents

## Repeatmasking

> Ceres: /lustre/project/gaur_genome_assembly/Rambouillet/PBJelly

```bash
sbatch --nodes=1 --mem=45000 --ntasks-per-node=50 -p priority -q msn --wrap="/lustre/project/rumen_longread_metagenome_assembly/binaries/RepeatMasker/RepeatMasker -pa 50 -q -species sheep -no_is -gff rambouillet.PBJelly.fa"

# Getting repeat families
perl -e '<>; <>; <>; while(<>){ $_ =~ s/^\s+//; @s = split(/\s+/); $orient = ($s[8] eq "+")? "+" : "-"; $qlen = $s[12] - $s[11]; print "$s[4]\t$s[5]\t$s[6]\t$orient\t$s[9]\t$s[10]\t$qlen\n";}' < rambouillet.PBJelly.fa.out > rambouillet.PBJelly.fa.famrepeat.bed
``` 

## Assembly QC

Assembly locations:

>/lustre/project/gaur_genome_assembly/Rambouillet/rambouillet_v1.9_contigs.fasta
>/lustre/project/gaur_genome_assembly/Rambouillet/NanoGrid_Rambouillet/rambouillet_v1.9_contigs.polished.fa
>/lustre/project/gaur_genome_assembly/Rambouillet/PurgeDups_config/rambouillet_v1.9_contigs.polished.purged.fa
>/lustre/project/gaur_genome_assembly/Rambouillet/Salsa/rambouillet_scaffolds_3/scaffolds_FINAL.fasta
>/lustre/project/gaur_genome_assembly/Rambouillet/scaffold_edits4/salsa_proposed_fixes4.fasta
>/lustre/project/gaur_genome_assembly/Rambouillet/rambouillet.contigs.polish.purge.salsa.pbjelly.polish.NUMT.MT.fasta
>/lustre/project/gaur_genome_assembly/Rambouillet/rambouillet.contigs.polish.purge.salsa.pbjelly.polish.NUMT.MT.haps.fasta

OK, setting up the assembly QC pipeline and letting it go to town!

> Ceres: /lustre/project/gaur_genome_assembly/Rambouillet/publication_qc/

```bash
module load miniconda/3.6

sbatch --nodes=1 --mem=5000 --ntasks-per-node=2 -p priority -q msn -t 8-0 snakemake --cluster-config ~/python_toolchain/snakeMake/assemblyValidation/cluster.json --cluster "sbatch --nodes={cluster.nodes} --ntasks-per-node={cluster.ntasks-per-node} --mem={cluster.mem} --partition={cluster.partition} -q {cluster.qos} -o {cluster.stdout} -t 8-0" -p --jobs 250 -s ~/python_toolchain/snakeMake/assemblyValidation/assemblyValidation --use-conda

bash /home/derek.bickharhth/python_toolchain/snakeMake/assemblyValidation/scripts/spectra-cn.revised.sh mapped/meryl_db.meryl /lustre/project/gaur_genome_assembly/Rambouillet/rambouillet.contigs.polish.purge.salsa.pbjelly.polish.NUMT.MT.fasta finalt finalt

mv finalt.* ./merqury/finalt/
Rscript /software/7/apps/merqury/1.0/plot/plot_spectra_cn.R -f merqury/finalt/finalt.spectra-asm.hist -o finalt.spectra-asm -z merqury/finalt/finalt.dist_only.hist --pdf

asm=finalt
name=finalt
asm_fa=fastas/final.fa
echo "# QV statistics"
ASM_ONLY=`meryl statistics merqury/finalt/$name.${asm}.0.meryl  | head -n4 | tail -n1 | awk '{print $2}'`
    TOTAL=`meryl statistics merqury/finalt/${asm}.meryl  | head -n4 | tail -n1 | awk '{print $2}'`
    ERROR=`echo "$ASM_ONLY $TOTAL" | awk -v k=$k '{print (1-(1-$1/$2)^(1/k))}'`
    QV=`echo "$ASM_ONLY $TOTAL" | awk -v k=$k '{print (-10*log(1-(1-$1/$2)^(1/k))/log(10))}'`
    echo -e "$asm\t$ASM_ONLY\t$TOTAL\t$QV\t$ERROR" >> merqury/$asm/$name.qv
echo

echo "# Per seq QV statistics"
meryl-lookup -existence -sequence $asm_fa -mers merqury/finalt/$name.$asm.0.meryl/ | \
awk -v k=21 '{print $1"\t"$NF"\t"$(NF-2)"\t"(-10*log(1-(1-$NF/$(NF-2))^(1/k))/log(10))"\t"(1-(1-$NF/$(NF-2))^(1/k))}' > merqury/finalt/$name.$asm.qv
echo
	read_solid=merqury/finalt/finalt.read.k21.finalt.gt4.meryl
    echo "# k-mer completeness (recoveray rate) with solid k-mers for $asm with > $filt counts"
meryl intersect output merqury/finalt/$asm.solid.meryl merqury/finalt/$asm.meryl $read_solid
    TOTAL=`meryl statistics $read_solid | head -n3 | tail -n1 | awk '{print $2}'`
    ASM=`meryl statistics merqury/finalt/$asm.solid.meryl | head -n3 | tail -n1 | awk '{print $2}'`
    echo -e "${asm}\tall\t${ASM}\t${TOTAL}" | awk '{print $0"\t"((100*$3)/$4)}' > merqury/$asm/$name.completeness.stats
rm -r $asm.solid.meryl
echo
```

#### Rerunning busco 3

The runtime and standards for busco 4 are all borked. I'm running on the previous version.

```bash
module load busco3/3.1.0

export BUSCO_CONFIG_FILE=/software/7/apps/busco3/3.1.0/config/config.ini.default
for i in fastas/*.fa; do name=`basename $i | cut -d'.' -f1`; echo $name; mkdir $name"_btemp"; cp -Rp /software/apps/augustus/gcc/64/3.2.3/config $name"_btemp/AUGUSTUS_CONFIG"; done

for i in fastas/*.fa; do name=`basename $i | cut -d'.' -f1`; echo $name; sbatch -N 1 -n 70 --mem=350000 -p priority -q msn --wrap="export AUGUSTUS_CONFIG_PATH=/lustre/project/gaur_genome_assembly/Rambouillet/rambouillet_qc/${name}_btemp/AUGUSTUS_CONFIG; run_BUSCO.py -i $i -o temp_${name} -c 70 -l /reference/data/BUSCO/v3/mammalia_odb9 -m genome"; done

for i in contigs ctgpolish ctgpolpurge final finhap oar4 ram1 scaffoldfixed scaffolds; do echo $i; cp run_temp_$i/short_summary*.txt busco/$i/busco_summary.txt; done

# NOTE: I had to run a snakemake instance with the --cleanup-metadata argument for two of the files

# Rerunning on freebayes_qc

export BUSCO_CONFIG_FILE=/software/7/apps/busco3/3.1.0/config/config.ini.default
for i in fastas/*.fa; do name=`basename $i | cut -d'.' -f1`; echo $name; mkdir $name"_btemp"; cp -Rp /software/apps/augustus/gcc/64/3.2.3/config $name"_btemp/AUGUSTUS_CONFIG"; done

for i in fastas/*.fa; do name=`basename $i | cut -d'.' -f1`; echo $name; sbatch -N 1 -n 70 --mem=350000 -p priority -q msn -t 6-0 --wrap="export AUGUSTUS_CONFIG_PATH=/lustre/project/gaur_genome_assembly/Rambouillet/freebayes_qc/${name}_btemp/AUGUSTUS_CONFIG; run_BUSCO.py -i $i -o temp_${name} -c 70 -l /reference/data/BUSCO/v3/mammalia_odb9 -m genome"; done
```


#### Testing kmer venn addition to Meryl

> Ceres: /lustre/project/gaur_genome_assembly/Rambouillet/freebayes_qc

```bash
conda activate /KEEP/rumen_longread_metagenome_assembly/meryl

sbatch -N 1 -n 4 --mem=30000 -p priority -q msn -t 2-0 --wrap="python3 ~/python_toolchain/sequenceData/merylVennUpset.py -m /lustre/project/rumen_longread_metagenome_assembly/binaries/meryl/build/bin/meryl -o ram_freebayes_comp -d merqury/finalt/finalt.meryl -d merqury/finhap/finhap.meryl -d merqury/freebayes/freebayes.meryl -d merqury/freehap/freehap.meryl -d merqury/freethap/freethap.meryl -d merqury/freetwo/freetwo.meryl -d merqury/ram1/ram1.meryl -d mapped/meryl_db.meryl"
```


#### Running individual Assembly QC on Rambouillet reads

> Ceres: /lustre/project/gaur_genome_assembly/Rambouillet/publication_qc/rambouillet_reads

```bash
sbatch --nodes=1 --mem=5000 --ntasks-per-node=2 -p priority -q msn -t 8-0 snakemake --cluster-config ~/python_toolchain/snakeMake/assemblyValidation/cluster.json --cluster "sbatch --nodes={cluster.nodes} --ntasks-per-node={cluster.ntasks-per-node} --mem={cluster.mem} --partition={cluster.partition} -q {cluster.qos} -o {cluster.stdout} -t 8-0" -p --jobs 250 -s ~/python_toolchain/snakeMake/assemblyValidation/assemblyValidation --use-conda

# The pipeline got screwed up
for i in oar4 ram1; do echo $i; sbatch depth.sh mapped/$i/merged.bam calls/$i/merged_depth.txt; done

module load python_3/3.6.6 miniconda/3.6
conda activate /KEEP/rumen_longread_metagenome_assembly/lumpy
for i in oar4 ram1 ram2; do echo $i; sbatch -N 1 -n 2 --mem=55000 -p priority -q msn --wrap="lumpyexpress -B mapped/$i/merged.bam -o calls/$i/merged_lumpy.vcf -v 2> logs/$i/lumpy_rerun.log"  ; done
```


#### Running individual Assembly QC on Texel reads

> Ceres: /lustre/project/gaur_genome_assembly/Rambouillet/publication_qc/texel_reads

```bash
sbatch --nodes=1 --mem=5000 --ntasks-per-node=2 -p priority -q msn -t 8-0 snakemake --cluster-config ~/python_toolchain/snakeMake/assemblyValidation/cluster.json --cluster "sbatch --nodes={cluster.nodes} --ntasks-per-node={cluster.ntasks-per-node} --mem={cluster.mem} --partition={cluster.partition} -q {cluster.qos} -o {cluster.stdout} -t 8-0" -p --jobs 250 -s ~/python_toolchain/snakeMake/assemblyValidation/assemblyValidation --use-conda
```