# Polishing reference assemblies using PB tools
---
*11/5/2017*

These are my notes on how to polish reference assemblies using PacBio tools on a commandline.


## Running arrow and PBJelly

I will be running these tools on Hank. First I need to generate a reference XML and use that for arrow. 

> Hank: /Jake/home/derek.bickhart/ARS-UCDv1.0.18_polish

```bash

# Download the reference
python download_from_googledrive.py 0BxbRXPCzWa5-amh4RUo0WXFxT1U ./ARS-UCDv1.0.18.fasta.gz

# Generate reference indicies
qsub call_fasta_to_ref.sh ARS-UCDv1.0.18.fasta ARS-UCDv1.0.18.ref

# Now to print out the list of pipelines and determine if they have the right settings
pbsmrtpipe show-templates
pbsmrtpipe show-template-details pbsmrtpipe.pipelines.sa3_ds_resequencing_fat -o resequencing_template.xml

# I now need to edit the resequencing_template.xml file to change the settings of the pipeline run
# Editted settings:
# "genomic_consensus.task_options.algorithm" to "arrow"
# And that was it. Crossing fingers

# Now to run the pipeline
pbsmrtpipe pipeline-id pbsmrtpipe.pipelines.sa3_ds_resequencing_fat --preset-xml resequencing_template.xml -e eid_subread:/ext/smrtlink/userdata/jobs_root/000/000632/tasks/pbcoretools.tasks.gather_subreadset-1/file.subreadset.xml -e eid_ref_dataset:ARS_UCDv1_0_18_ref/referenceset.xml

# It exited with the following error:
WARNING:pbsmrtpipe.pb_io:UNKNOWN Task Option with id 'genomic_consensus.task_options.diploid'. Ignoring option
successfully initialized datastore.
...
Starting worker pbcoretools.tasks.filterdataset-0 (1 workers running, 1 total proc in use)
ls: cannot access /Jake/home/derek.bickhart/ARS-UCDv1.0.18_polish/tasks/pbcoretools.tasks.filterdataset-0/stderr: No such file or directory
Task pbcoretools.tasks.filterdataset-0 task pbcoretools.tasks.filterdataset failed (exit-code 1) after 0.52 sec

Unable to run job: Job was rejected because job requests unknown queue "default".
Exiting.

Workflow status 2/13 completed/total tasks
Shutting down.
Completed execution pbsmrtpipe v0.51.2. Workflow Failed in 4.84 sec (0.08 min) with exit code 1

##############
# Retrying   #
##############

# OK, so the problem is that the queue is wrong! I need to change the start and stop values for the SGE to proceed in a workflows xml file
mkdir sge_smrtq
cp /ext/smrtlink/install/smrtlink-release_5.0.1.9585/bundles/smrttools/install/smrttools-release_5.0.1.9578/private/pacbio/pythonpkgs/pbsmrtpipe/lib/python2.7/site-packages/pbsmrtpipe/cluster_templates/sge_pacbio_def66/start.tmpl ./sge_smrtq/
cp /ext/smrtlink/install/smrtlink-release_5.0.1.9585/bundles/smrttools/install/smrttools-release_5.0.1.9578/private/pacbio/pythonpkgs/pbsmrtpipe/lib/python2.7/site-packages/pbsmrtpipe/cluster_templates/sge_pacbio_def66/stop.tmpl ./sge_smrtq/

# I then modified the start.tmpl to use the correct queue: smrt.q
# now to generate a workflow xml
pbsmrtpipe show-workflow-options -o resequencing_workflow_template.xml

# I changed ""pbsmrtpipe.options.cluster_manager"" to /Jake/home/derek.bickhart/ARS-UCDv1.0.18_polish/sge_smrtq/
# OK, running it again with the change in workflows
pbsmrtpipe pipeline-id pbsmrtpipe.pipelines.sa3_ds_resequencing_fat --preset-xml resequencing_template.xml --preset-xml resequencing_workflow_template.xml -e eid_subread:/ext/smrtlink/userdata/jobs_root/000/000632/tasks/pbcoretools.tasks.gather_subreadset-1/file.subreadset.xml -e eid_ref_dataset:ARS_UCDv1_0_18_ref/referenceset.xml
```

#### PBjelly reattempt

I am going to reattempt PBjelly on the files using scripts that Serge provided. We noticed that the pacbio tools hate "." and "-" delimiters in file names (apart from the last extension) so I am preemptively changing those delimiters to avoid any problems down the road.

> Hank: /Jake/home/derek.bickhart/ARS_UCDv1_0_19_polish

```bash
mv ARS-UCD1.0.19.base.fasta ARS_UCD1_0_19_base.fasta

cp ../PBSuite_15.8.24/docs/TemplateProtocol.xml ./
mv TemplateProtocol.xml Protocol.xml
mkdir ARS_v19_jellied

vim Protocol.xml
```