# Sequenza Workflow
Workflow for allele-specific copy number analysis from tumour and matched normal sequencing data 

Takes a tumour BAM file and matched normal BAM file as inputs, and outputs allele-specific copy numbers and related R plots. Pipeline was developed by Center for Biological Sequence Analysis, Technical University of Denmark and can be found in the [Sequenza Bitbucket](https://sequenzatools.bitbucket.io/#/home).

The original workflow has been modified to not create and output an archive for seqz partial files (parts_seqz.tar.gz) as it was not needed for our purposes. 

---

**sequenza.inputs.json:** 

* Replace "sample_ID" with the hypenated sample IDs of the tumor and normal samples you want to analyze.
* Replace "normal_reads.bam" with the full filepath (or SAS token) to the input matched normal BAM.
* Replace "tumor_reads.bam" with the full filepath (or SAS token) to the input tumor BAM.

**sequenza.trigger.json:**

* Replace "WorkflowUrl" with the URL to either a local version of the WDL (in an Azure Storage Account), or the URL to the version available in this repository online.

* Replace "WorkflowInputsUrl" with the URL to a local version of the inputs.json file (in an Azure Storage Account).

* Optional: Replace "WorkflowOptionsUrl" and/or "WorkflowDependenciesUrl" with the URL to a local version of the options.json and/or dependencies.json files, respectively.

**sequenza.wdl:**

* No changes necessary.

---

### Expected Running Time
* To be determined. 
