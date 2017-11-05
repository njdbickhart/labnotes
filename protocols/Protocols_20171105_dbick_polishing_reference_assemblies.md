# Polishing reference assemblies using PB tools
---
*11/5/2017*

These are my notes on how to polish reference assemblies using PacBio tools on a commandline.


## Running arrow and PBJelly

I will be running these tools on Hank. First I need to generate a reference XML and use that for arrow. 

> Hank:

```bash
# Download the reference
python download_from_googledrive.py 0BxbRXPCzWa5-amh4RUo0WXFxT1U ./ARS-UCDv1.0.18.fasta.gz

# Generate reference indicies
qsub call_fasta_to_ref.sh ARS-UCDv1.0.18.fasta ARS-UCDv1.0.18.ref
