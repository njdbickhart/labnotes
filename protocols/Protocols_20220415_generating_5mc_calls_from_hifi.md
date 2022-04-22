# HiFi 5mc conversion

## Workflow

* Generate HiFi reads + kinetics from subreads
* Call bases from HiFi bam + kinetics to generate annotated bam file
* Use bam annotations in downstream applications

## Trial run on Ceres

> Ceres: /90daydata/project/sheep_genome_assemblies/Churro_x_Friesian/HiFi_data/subreads.bam

```bash
module load python_3/3.6.6 miniconda/3.6 samtools minimap2
for i in `seq  1 10`; do echo $i; sbatch -N 1 -n 30 --mem=100000 -p priority -q msn --wrap="ccs -j 30 --chunk ${i}/10 --hifi-kinetics m54337U_211106_060943.subreads.bam m54337U_211106_060943.hifikin.${i}.bam"; done

conda activate /project/rumen_longread_metagenome_assembly/environments/primrose/

for i in `seq  1 10`; do echo $i; sbatch -N 1 -n 30 --mem=100000 -p priority -q msn --wrap="primrose -j 30 m54337U_211106_060943.hifikin.${i}.bam m54337U_211106_060943.primrose.${i}.bam"; done

# Now to align to the reference genome
sbatch -N 1 -n 25 --mem=150000 -p priority -q msn --wrap="pbmm2 align -j 18 -J 7 --sort --log-level INFO --preset HIFI /90daydata/sheep_genome_assemblies/sergek/verkko_beta2/8-trio/gapped/renamed_gapped.fasta m54337U_211106_060943.primrose.bam sheept2t_test.primrose.bam"
```