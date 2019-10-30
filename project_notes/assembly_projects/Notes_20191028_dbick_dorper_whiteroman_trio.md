# Dorper White Romanov cross assembly
---
*10/28/2019*

These are my notes on assisting with the Dorper x White Romanov assembly project.

## Table of contents


## Generating trio canu corrected reads

I first want to bin all of the reads from Tim's pacbio run.

> Ceres: /project/forage_assemblies/assemblies/romanov_whitedorp

```bash
module load canu/1.8 samtools
sbatch --nodes=1 --ntasks-per-node=10 --mem=10G -q msn -p msn --wrap="canu -p dorp_wroman -d dorp_wroman genomeSize=2800m -haplotypeSire /project/gaur_genome_assembly/WhiteDorper_x_Romanov/Dorper_sire.fastq  -haplotypeDam /project/gaur_genome_assembly/WhiteDorper_x_Romanov/Romanov_dam.fastq corMhapSensitivity=normal corOutCoverage=200 saveReadCorrections=true 'gridOptions=-q msn -p msn' -pacbio-raw /project/gaur_genome_assembly/WhiteDorper_x_Romanov/pacbio_data/all_reads.fq"

# I just cancelled the queueing jobs and now I'm going to queue up two Flye tasks
# The location of the separated files
ls dorp_wroman/haplotype/
0-kmers                 haplotype.log            haplotype-unknown.fasta.gz    splitHaplotype.jobSubmit-01.out  splitHaplotype.sh
haplotype-Dam.fasta.gz  haplotype-Sire.fasta.gz  splitHaplotype.1268593_1.out  splitHaplotype.jobSubmit-01.sh

# OK, now to queue up the two separate flye jobs
module load miniconda
source activate /KEEP/rumen_longread_metagenome_assembly/flye

sbatch --nodes=1 --mem=300000 --ntasks-per-node=70 -p msn -q msn flye --pacbio-raw /project/forage_assemblies/assemblies/romanov_whitedorp/dorp_wroman/haplotype/haplotype-Dam.fasta.gz -g 2800m -t 70 -i 2 -m 10000 --asm-coverage 40 -o rw_flye_dam

sbatch --nodes=1 --mem=300000 --ntasks-per-node=70 -p msn -q msn flye --pacbio-raw /project/forage_assemblies/assemblies/romanov_whitedorp/dorp_wroman/haplotype/haplotype-Sire.fasta.gz -g 2800m -t 70 -i 2 -m 10000 --asm-coverage 40 -o rw_flye_sire

# The haplotype resolved reads had duplicate IDs! Checking and removing duplicates
sbatch check_names.pl /project/forage_assemblies/assemblies/romanov_whitedorp/dorp_wroman/haplotype/haplotype-Sire.fasta.gz haplotype-Sire.nonuniq
sbatch check_names.pl /project/forage_assemblies/assemblies/romanov_whitedorp/dorp_wroman/haplotype/haplotype-Dam.fasta.gz haplotype-Dam.nonuniq

sbatch --nodes=1 --mem=10000 --ntasks-per-node=2 -p msn -q msn --wrap='python3 ~/python_toolchain/sequenceData/filterFastaqFromList.py -f dorp_wroman/haplotype/haplotype-Dam.fasta.gz -l haplotype-Dam.nonuniq -o haplotype_Dam.uniq.fasta -v'
sbatch --nodes=1 --mem=10000 --ntasks-per-node=2 -p msn -q msn --wrap='python3 ~/python_toolchain/sequenceData/filterFastaqFromList.py -f dorp_wroman/haplotype/haplotype-Dam.fasta.gz -l haplotype-Dam.nonuniq -o haplotype_Dam.nonuniq.fasta'
sbatch --nodes=1 --mem=10000 --ntasks-per-node=2 -p msn -q msn --wrap='python3 ~/python_toolchain/sequenceData/filterFastaqFromList.py -f dorp_wroman/haplotype/haplotype-Sire.fasta.gz -l haplotype-Sire.nonuniq -o haplotype_Sire.uniq.fasta -v'
sbatch --nodes=1 --mem=10000 --ntasks-per-node=2 -p msn -q msn --wrap='python3 ~/python_toolchain/sequenceData/filterFastaqFromList.py -f dorp_wroman/haplotype/haplotype-Sire.fasta.gz -l haplotype-Sire.nonuniq -o haplotype_Sire.nonuniq.fasta'

# The nonuniq files have the duplicated sequence still. Let's see if that is present in true duplicated (ordered) fashion
head -n 5055054 haplotype_Sire.nonuniq.fasta > haplotype_Sire.nowuniq.fasta
head -n 5894688 haplotype_Dam.nonuniq.fasta > haplotype_Dam.nowuniq.fasta

sbatch --nodes=1 --mem=1000 --ntasks-per-node=1 -p msn -q msn --wrap='cat haplotype_Dam.nowuniq.fasta haplotype_Dam.uniq.fasta > haplotype_Dam.total.fasta'
sbatch --nodes=1 --mem=1000 --ntasks-per-node=1 -p msn -q msn --wrap='cat haplotype_Sire.nowuniq.fasta haplotype_Sire.uniq.fasta > haplotype_Sire.total.fasta'

sbatch --dependency=afterok:1270592 --nodes=1 --mem=300000 --ntasks-per-node=70 -p msn -q msn flye --pacbio-raw haplotype_Dam.total.fasta -g 2800m -t 70 -i 2 -m 10000 --asm-coverage 40 -o rw_flye_dam

sbatch --dependency=afterok:1270595 --nodes=1 --mem=300000 --ntasks-per-node=70 -p msn -q msn flye --pacbio-raw haplotype_Sire.total.fasta -g 2800m -t 70 -i 2 -m 10000 --asm-coverage 40 -o rw_flye_sire
```