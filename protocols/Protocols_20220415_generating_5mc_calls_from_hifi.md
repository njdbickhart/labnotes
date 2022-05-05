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

conda activate /project/rumen_longread_metagenome_assembly/environments/pb_cpg_tools
for i in count model; do echo $i; sbatch -N 1 -n 36 --mem=150000 -p priority -q msn --wrap="python3 ~/rumen_longread_metagenome_assembly/binaries/pb-CpG-tools/aligned_bam_to_cpg_scores.py -b sheept2t_test.primrose.bam -f /90daydata/sheep_genome_assemblies/sergek/verkko_beta2/8-trio/gapped/renamed_gapped.fasta -o m54337U_211106_060943.pbcpg.${i} -p $i -d /project/rumen_longread_metagenome_assembly/binaries/pb-CpG-tools/pileup_calling_model -t 36"; done
```

I have some options for visualization to compare methylation predicted near gene regulatory regions. First I will need to generate some annotation tracks for genes and then plot windows containing clusters of those tracks in a wider area.

```bash
module load python_3/3.6.6 miniconda/3.6 minimap2
wget https://hgdownload.soe.ucsc.edu/gbdb/oviAri4/ncbiRefSeq/seqNcbiRefSeq.rna.fa
sbatch -N 1 -n 3 --mem=35000 -p priority -q msn --wrap='minimap2 -x asm20 /90daydata/sheep_genome_assemblies/sergek/verkko_beta2/8-trio/gapped/renamed_gapped.fasta seqNcbiRefSeq.rna.fa > seqNcbiRefSeq.paf'

python3 ~/python_toolchain/sequenceData/pafToBed12.py -p seqNcbiRefSeq.paf -o seqNcbiRefSeq.sheep.bed

```

Now let's try [pyGenomeTracks](https://github.com/deeptools/pyGenomeTracks) to plot methylation calls as a window within the genome.

```bash
module load python_3/3.6.6 miniconda/3.6 minimap2

conda activate /project/rumen_longread_metagenome_assembly/environments/pygenometracks
make_tracks_file --trackFiles m54337U_211106_060943.pbcpg.count.combined.denovo.bw m54337U_211106_060943.pbcpg.model.combined.denovo.bw seqNcbiRefSeq.sheep.bed --out primrose_pygenometracks.ini

# Now to edit the INI

# Now to run the script
pyGenomeTracks --region path_from_utig4-1000:6370000-65700000 --tracks primrose_pygenometracks.ini --outFileName utig4_1000_tracks.png

# It looks like the "-" in the chromosome name confused the plot. Let's try with a bed file instead
echo -e "path_from_utig4-1000\t63700000\t65700000" > test_plot.bed
pyGenomeTracks --BED test_plot.bed --tracks primrose_pygenometracks.ini --outFileName utig4_1000_tracks.png

# Testing multiple plots from haplotigs
echo -e "path_from_utig4-1746\t67251409\t67415759\npath_from_utig4-1000\t64647559\t64811972" > test_plot.bed
pyGenomeTracks --BED test_plot.bed --tracks primrose_pygenometracks.ini --outFileName haplotig_test.png
```