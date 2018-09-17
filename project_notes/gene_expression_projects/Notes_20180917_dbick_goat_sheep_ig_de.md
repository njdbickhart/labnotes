# Goat and Sheep IG differential expression analysis
---
*9/17/2018*

These are my notes on running the EBseq pipeline on John H and John S's differential expression data on Goat and Sheep.


## Preparing the data

So, I first need to format the data into a useable matrix 

```bash
for i in *.tab; do echo $i; dos2unix $i; done
```