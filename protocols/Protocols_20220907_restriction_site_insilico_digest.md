# In silico RE digest for RRS assay development
---

## Table of Contents

## Setting up the environment

For this exercise, I will be using Biopython to generate a list of RE sites. The goal is to generate the list of sites and then export the cut sites on a per-chromosome basis so that I can calculate the fragment sizes locally. 

> Anunna: /lustre/backup/HG/Sequence/PROTIX/restriction

```bash
conda deactivate
conda activate /lustre/backup/HG/bickha01/conda_envs/biopython
```

```python
from Bio.Restriction import *
from Bio.Seq import Seq
from Bio import SeqIO

# Picking a bunch of small (6bp or less) REs with blunt end cutting
rb = RestrictionBatch([SspI, ScaI, PvuII, SmaI, EcoRV, DpnI, HaeIII, HpaI, SnaBI, FspI, DraI])

chrs = dict()
with open("/lustre/backup/HG/Sequence/PROTIX/genome/GCF_905115235.1_iHerIll2.2.curated.20191125_genomic.fna", "r") as seqio:
     for srec in SeqIO.parse(seqio, "fasta"):
             analysis = Analysis(rb, srec.seq)
             chrs[srec.id] = analysis.full()
             print(srec.id)

with open("restriction_sites.tab","w") as output:
     for c, rs in chrs.items():
             for rn, rl in rs.items():
                     output.write(f'{c}\t{rn}\t{rl}\n')
```

Now that should have printed all of the restriction sites for the selected enzymes in a tab delimited format with the columns being:

* Chr name
* Restriction enzyme name
* [list of all cut positions in the sequence]