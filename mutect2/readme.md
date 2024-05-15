# Mutect2 Workflow
Workflow for somatic variant calling; outputs both VCF and MAF format functionally annotated files.

Source GATK version of this workflow can be found at the [GATK GitHub repository](https://github.com/broadinstitute/gatk/tree/master/scripts/mutect2_wdl) for mutect2 workflows.

This workflow follows the general outline shown below:

<p align="center"><img src="https://user-images.githubusercontent.com/107152811/181548163-3fd1b990-e8dc-428a-be39-678e26b9ab6c.PNG" width="550"></p>

---

**mutect2.inputs.json:** 

* Replace "Mutect2.normal_reads" with the filepath to the BAM file containing matched normal reads (associated with the tumor reads).
* Replace "Mutect2.normal_reads_index" with the filepath to the index file (.bai) of the matched normal reads.
* Replace "Mutect2.tumor_reads" with the filepath to the BAM file containing tumor reads.
* Replace "Mutect2.tumor_reads_index" with the filepath to the index file (.bai) of the tumor reads.
* Replace lines 17-23 with the filepaths to the appropriate reference file listed.
* Replace "Mutect2.funco_data_sources_tar_gz" with the filepath to the tar and gzipped Funcotator references (see Confluence page on Mutect2 for more information).


**mutect2.trigger.json:**

* Replace "WorkflowUrl" with the URL to either a local version of the WDL (in an Azure Storage Account), or the URL to the version available in this repository online.

* Replace "WorkflowInputsUrl" with the URL to a local version of the inputs.json file (in an Azure Storage Account), updated with the above-mentioned normal and tumor reads/indices.

* Optional: Replace "WorkflowOptionsUrl" and/or "WorkflowDependenciesUrl" with the URL to a local version of the options.json and/or dependencies.json files, respectively.


**Mutect2.wdl Versions:**

*Primary WDLs:*
* **mutect2-cbio.wdl:** Default for cfDNA and other tumor samples.
* **mutect2-cbio-BigBam.wdl:** For combined BAM files exceeding 150GB.
* **mutect2-cbio-FFPE.wdl:** Specifically designed for FFPE samples, includes additional FFPE filtering.

*Other Variants:*
* **mutect2-cbio-more-funcotate.wdl:** Enhanced runtime metrics for Funcotate task. Includes an additional 50MB of memory and 10GB of disk space. (Code lines: 1037-1039, 1093, 1097)
* **mutect2-cbio-more-M2.wdl:** Enhanced runtime metrics for M2 task. Increased disk space allocation. (Code lines: 166, 545, 620, 623)
* **mutect2-cbio-more-maffuncotate.wdl:** Enhanced runtime metrics for MafFuncotate task. Increased memory, disk, and CPU allocations. (Code lines: 1158-1159, 1220-1221, 1224)
* **mutect2-cbio-odd.wdl:** Enhanced runtime metrics for M2, FilterAlignmentArtifacts, Funcotate, and MafFuncotate tasks. This version has extended runtime. (Code lines: 157, 168, 208, 547, 622, 625, 959, 1007, 1009, 1137-1138)

---

### Expected Running Time
* For a 85GB normal and a 143GB tumor input, mutect2 will take ~6.5 hours to finish.
