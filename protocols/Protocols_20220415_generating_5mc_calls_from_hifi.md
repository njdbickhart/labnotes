# HiFi 5mc conversion

## Workflow

* Generate HiFi reads + kinetics from subreads
* Call bases from HiFi bam + kinetics to generate annotated bam file
* Use bam annotations in downstream applications

## Trial run on Ceres

> Ceres: /90daydata/project/sheep_genome_assemblies/Churro_x_Friesian/HiFi_data/subreads.bam

```bash
for i in `seq  1 10`; do echo $i; sbatch -N 1 -n 30 --mem=100000 -p priority -q msn --wrap="ccs -j 30 --chunk ${i}/10 --hifi-kinetics m54337U_211106_060943.subreads.bam m54337U_211106_060943.hifikin.${i}.bam"; done
```