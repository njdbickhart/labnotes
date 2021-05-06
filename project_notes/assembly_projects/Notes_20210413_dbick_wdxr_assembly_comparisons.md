# White Dorper x Romanov genome assembly comparison
---
*5/20/2019*

These are my notes on the comparative assembly assessment of our WDxR trio.

## Table of Contents

## Themis-ASM run and comparisons

> Ceres: /lustre/project/rumen_longread_metagenome_assembly/wdxrom

```bash
# Downloading the assemblies from NCBI first
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/298/735/GCA_000298735.2_Oar_v4.0/GCA_000298735.2_Oar_v4.0_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/170/295/GCA_011170295.1_ASM1117029v1/GCA_011170295.1_ASM1117029v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/765/115/GCA_000765115.1_Oori1/GCA_000765115.1_Oori1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/121/645/GCA_003121645.1_O_ammon_KGZ_v1.0/GCA_003121645.1_O_ammon_KGZ_v1.0_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/017/524/585/GCA_017524585.1_CAU_O.aries_1.0/GCA_017524585.1_CAU_O.aries_1.0_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/039/535/GCA_001039535.1_ASM103953v1/GCA_001039535.1_ASM103953v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/903/231/385/GCA_903231385.1_OvNiv1.0/GCA_903231385.1_OvNiv1.0_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/014/523/465/GCA_014523465.1_CAU_Oori_1.0/GCA_014523465.1_CAU_Oori_1.0_genomic.fna.gz


for i in *.gz; do echo $i; sbatch -N 1 -n 1 --mem=5000 -p priority -q msn --wrap="gunzip $i"; done

mkdir assemblies
cp *a ./assemblies/

cp -r /lustre/project/cattle_genome_assemblies/bison_x_simmental/bison_x_simmental_qc/bison_manuscript_qc/busco_downloads ./


module load miniconda/3.6 "java/11.0.2" bedtools

PATH=$PATH:/software/7/apps/merqury/1.1/:/software/7/apps/meryl/1.0/Linux-amd64/bin/
MERQURY=/software/7/apps/merqury/1.1

sbatch -N 1 -n 2 --mem=6000 -p priority -q msn --wrap='python3 ~/rumen_longread_metagenome_assembly/binaries/Themis-ASM/themisASM.py -a /lustre/project/rumen_longread_metagenome_assembly/wdxrom/assemblies/Dorper_sire_nano_purge_salsa_fix_arrow_MT_fb2.fasta -n Dorper -a /lustre/project/rumen_longread_metagenome_assembly/wdxrom/assemblies/GCA_000298735.2_Oar_v4.0_genomic.fna -n OAR4 -a /lustre/project/rumen_longread_metagenome_assembly/wdxrom/assemblies/GCA_000765115.1_Oori1_genomic.fna -n Mouflon -a /lustre/project/rumen_longread_metagenome_assembly/wdxrom/assemblies/GCA_001039535.1_ASM103953v1_genomic.fna -n BigHorn -a /lustre/project/rumen_longread_metagenome_assembly/wdxrom/assemblies/GCA_003121645.1_O_ammon_KGZ_v1.0_genomic.fna -n Argili -a /lustre/project/rumen_longread_metagenome_assembly/wdxrom/assemblies/GCA_011170295.1_ASM1117029v1_genomic.fna -n Hu -a /lustre/project/rumen_longread_metagenome_assembly/wdxrom/assemblies/GCA_014523465.1_CAU_Oori_1.0_genomic.fna -n AMouflon -a /lustre/project/rumen_longread_metagenome_assembly/wdxrom/assemblies/GCA_017524585.1_CAU_O.aries_1.0_genomic.fna -n Tibetan -a /lustre/project/rumen_longread_metagenome_assembly/wdxrom/assemblies/GCA_903231385.1_OvNiv1.0_genomic.fna -n SnowSheep -a /lustre/project/rumen_longread_metagenome_assembly/wdxrom/assemblies/Romanov_dam_nano_purge_salsa_fix_arrow_MT_fb2.fasta -n Romanov -a /lustre/project/sheep_genome_assemblies/Rambouillet/GCA_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna -n Rambouillet -b mammalia_odb10 -f /lustre/project/sheep_genome_assemblies/WhiteDorper_x_Romanov/illumina_data/Romanov_dam/LIB102054_S1_L001_R1_001.fastq.gz,/lustre/project/sheep_genome_assemblies/WhiteDorper_x_Romanov/illumina_data/Romanov_dam/LIB102054_S1_L001_R2_001.fastq.gz -s RomDam -f /lustre/project/sheep_genome_assemblies/WhiteDorper_x_Romanov/illumina_data/White_Dorper_sire/LIB102055_S1_L001_R1_001.fastq.gz,/lustre/project/sheep_genome_assemblies/WhiteDorper_x_Romanov/illumina_data/White_Dorper_sire/LIB102055_S1_L001_R2_001.fastq.gz -s WDSire -c "sbatch --nodes={cluster.nodes} --ntasks-per-node={cluster.ntasks-per-node} --mem={cluster.mem} --partition={cluster.partition} -q {cluster.qos} -o {cluster.stdout}" -j 100 --resume'
```

#### Meryl Venn analysis

```bash
conda activate /KEEP/rumen_longread_metagenome_assembly/meryl

for i in AMouflon Argili BigHorn Dorper Hu Mouflon OAR4 Rambouillet Romanov SnowSheep Tibetan; do echo -n "-d merqury/$i/$i.meryl "; done; echo

sbatch -N 1 -n 4 --mem=50000 -p priority -q msn -t 2-0 --wrap="python ~/python_toolchain/sequenceData/merylVennUpset.py -m /lustre/project/rumen_longread_metagenome_assembly/binaries/meryl/build/bin/meryl -o sheep_venn_kmer_comp -d merqury/AMouflon/AMouflon.meryl -d merqury/Argili/Argili.meryl -d merqury/BigHorn/BigHorn.meryl -d merqury/Dorper/Dorper.meryl -d merqury/Hu/Hu.meryl -d merqury/Mouflon/Mouflon.meryl -d merqury/OAR4/OAR4.meryl -d merqury/Rambouillet/Rambouillet.meryl -d merqury/Romanov/Romanov.meryl -d merqury/SnowSheep/SnowSheep.meryl -d merqury/Tibetan/Tibetan.meryl"
```