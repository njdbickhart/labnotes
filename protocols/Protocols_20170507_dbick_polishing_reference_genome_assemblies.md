# Polishing reference genome assemblies
---
*5/7/2017*

## Table of Contents

<a name="pbjelly"></a>
## PbJelly run

First things first, we need to polish the reference assembly using the filtered subreads and pbJelly from PBSuite. I have installed PBsuite, with the following notes:

* It needs Blasr v1.3.1
* Blasr is a pain to install locally -- friends don't let friends install Blasr!
* You need to install networkx via pip
* NOTE: you need to have a reference "FASTA" extension for jelly to work!

Now to generate the protocol.xml. The working directory:

> hank: /Jake/home/derek.bickhart/ars_ucd_asm

```xml
<jellyProtocol>
    <reference>/Jake/home/derek.bickhart/ars_ucd_asm/ARS-UCD1.0.9.fasta</reference>
    <outputDir>/Jake/home/derek.bickhart/ars_ucd_asm/pbjelly/</outputDir>
    <blasr>-minMatch 8 -minPctIdentity 70 -bestn 1 -nCandidates 20 -maxScore -500 -nproc 15 -noSplitSubreads</blasr>
    <input baseDir="/Jake/PacBio/transfer/">
        <job>Dominette_genome_USMARC_179cells.fastq</job>
        <job>UC_Davis_filtered_subreads.fastq</job>
    </input>
</jellyProtocol>
```

And now to submit via SGE.

```bash
qsub pbjelly_run_sge.sh
```

