# Identifying IGC regions that are resolved from optical mapping
---
*10/14/2015*

These are my notes on doing a very rapid identification of immune gene cluster (IGC) exons that span different PacBio contigs, but are linked by optical mapping technologies -- namely BioNano Genomics Irys!


## Preparing materials for the analysis

I'm going to have to get exons of the genes that I'm interested in interrogating. The first step is getting Human and Cattle exon sequence. Human is likely better annotated, but Cattle will be the closest/highest quality ruminant reference genome sequence to Goat. I will then map both sets to the Goat reference assembly and try to identify exons that span different contigs. Then I will see if any of those exons lie within Super-scaffolds.

First, let's get the fasta files for the list of contigs from Ensembl.

OK, ran into a problem: to restrict bandwidth, Ensembl restricts gene accession to EnsGene and EnsTranscript IDs. Here is the XML of how I retrieved the Human genes with GO type "Immune response:"

```xml
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
			
	<Dataset name = "hsapiens_gene_ensembl" interface = "default" >
		<Filter name = "transcript_source" value = "ensembl"/>
		<Filter name = "go_parent_name" value = "immune response"/>
		<Filter name = "transcript_biotype" value = "protein_coding"/>
		<Attribute name = "ensembl_gene_id" />
		<Attribute name = "ensembl_transcript_id" />
	</Dataset>
</Query>
```

There were 400+ genes with that GO annotation. Now for Cattle.

```xml
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
			
	<Dataset name = "btaurus_gene_ensembl" interface = "default" >
		<Filter name = "transcript_source" value = "ensembl"/>
		<Filter name = "go_parent_name" value = "immune response"/>
		<Filter name = "transcript_biotype" value = "protein_coding"/>
		<Attribute name = "ensembl_gene_id" />
		<Attribute name = "ensembl_transcript_id" />
	</Dataset>
</Query>
```

Strange, cattle had 600 genes this time. Oh well, I'm going to try pulling gene lists from the UCSC.

The UCSC genome browser does have a nice list of exon coordinates for the assembly, but it does not contain the exon fastas themselves. I'll have to generate that myself with a custom script. I wrote a script [here](https://github.com/njdbickhart/perl_toolchain/blob/master/bed_cnv_fig_table_pipeline/generateExonFastaFromUCSCTable.pl) that should take care of the problem.

> pwd: /home/dbickhart/share/goat_assembly_paper/igc_annotation

```bash
perl -e '<>; while(<>){chomp; @s = split(/\t/); print "$s[1]\n";}' < cattle_ensid_immune_response.txt > cattle_enstranscript_immune_response.txt

perl ../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/generateExonFastaFromUCSCTable.pl -b umd3EnsGene.txt -f ../../umd3_data/umd3_kary_unmask_ngap.fa -e cattle_enstranscript_immune_response.txt -o cattle_enstranscript_immune_response.fa
```

I'm really curious to see how this turns out. Let's start the alignment and see how the data piles up. I'm going to use the frozen BNG version 5 assembly in order to test this out -- I should have the super scaffolds nicely broken up and tagged already.

> Blade14: /mnt/iscsi/vnx_gliu_7/goat_assembly/annotation

```bash
bwa mem ../papadum-v5bng-pilon-split2kbgap.fa cattle_enstranscript_immune_response.fa | samtools view -bS - | samtools sort -T cattle.temp -o cattle_igc_goat_align.bam -

samtools index cattle_igc_goat_align.bam
samtools flagstat cattle_igc_goat_align.bam
	5275 + 0 in total (QC-passed reads + QC-failed reads)
	0 + 0 secondary
	35 + 0 supplementary
	0 + 0 duplicates
	5210 + 0 mapped (98.77%:-nan%)
	0 + 0 paired in sequencing
	0 + 0 read1
	0 + 0 read2
	0 + 0 properly paired (-nan%:-nan%)
	0 + 0 with itself and mate mapped
	0 + 0 singletons (-nan%:-nan%)
	0 + 0 with mate mapped to a different chr
	0 + 0 with mate mapped to a different chr (mapQ>=5)

# So there are 35 that have multiple mappings
# Lots have aligned to the super-reads -- let's see if we can find alternative mappings

samtools view cattle_igc_goat_align.bam | perl -e '%h; while(<>){chomp; @s = split(/\t/); $s[0] =~ s/\_\d{1,3}$//; $h{$s[0]}->{$s[2]} = 1;} foreach my $contig (keys %h){foreach my $chr (keys %{$h{$contig}}){ print "$contig\t$chr\n";}}' | less
ENSBTAT00000007438      Scaffold_154.1
ENSBTAT00000007438      Scaffold_94.1

ENSBTAT00000007743      Scaffold_94.1
ENSBTAT00000007743      Scaffold_134.1

# Very good example Interleukin 1 alpha
ENSBTAT00000013665      Scaffold_22.3
ENSBTAT00000013665      Scaffold_22.1
```