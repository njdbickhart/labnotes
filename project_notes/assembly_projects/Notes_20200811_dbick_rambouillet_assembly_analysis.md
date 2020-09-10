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

> Ceres:

```bash
module load miniconda/3.6

sbatch --nodes=1 --mem=5000 --ntasks-per-node=2 -p priority -q msn snakemake --cluster-config ~/python_toolchain/snakeMake/assemblyValidation/cluster.json --cluster "sbatch --nodes={cluster.nodes} --ntasks-per-node={cluster.ntasks-per-node} --mem={cluster.mem} --partition={cluster.partition} -q {cluster.qos} -o {cluster.stdout}" -p --jobs 250 -s ~/python_toolchain/snakeMake/assemblyValidation/assemblyValidation --use-conda

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