version 1.0

## Copyright Broad Institute, 2018 and Aparicio Lab (BC Cancer Research Centre), 2022
##
## This WDL pipeline implements data pre-processing according to the GATK Best Practices 
## (June 2016) for human whole-genome data.
##
## Requirements/expectations :
## - Human whole-genome pair-end sequencing data in unmapped BAM (uBAM) format
## - One or more read groups, one per uBAM file, all belonging to a single sample (SM)
## - Input uBAM files must additionally comply with the following requirements:
## - - filenames all have the same suffix (we use ".unmapped.bam")
## - - files must pass validation by ValidateSamFile
## - - reads are provided in query-sorted order
## - - all reads must have an RG tag
## - Reference genome must be Hg38 with ALT contigs
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.
## 
## UPDATE NOTES:
## Last Updated: September 20, 2023 by Matthew Cho (Aparicio Lab) - fixed bug in GetBamHeader task
## 
##
## This pipeline has been modified from its original, which can be found at 
## https://github.com/microsoft/gatk4-genome-processing-pipeline-azure. Major changes include
## removing all germline SNP/indel calling functionality; pipeline is now just used for
## converting unmapped BAM files (uBAMs) into analysis-ready BAM files, that can be
## used in later analysis (ex., somatic variant calling); and renaming the workflow to "PreProcessing".
## It also includes the call to BamToCram for Cram output files.

import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/custom-pre-pro/tasks/UnmappedBamToAlignedBam.wdl" as ToBam
import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/custom-pre-pro/tasks/BamToCram.wdl" as ToCram
import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/custom-pre-pro/tasks/Qc.wdl" as QC
import "https://raw.githubusercontent.com/microsoft/gatk4-genome-processing-pipeline-azure/az1.1.0/structs/GermlineStructs.wdl"

# WORKFLOW DEFINITION
workflow PreProcessing {

  String pipeline_version = "1.4"

  input {
    String study # metadata for clean-up automation 
    
    SampleAndUnmappedBams sample_and_unmapped_bams
    GermlineSingleSampleReferences references
    PapiSettings papi_settings
    File wgs_coverage_interval_list

    File? haplotype_database_file
    Boolean provide_bam_output = true
    Boolean use_gatk3_haplotype_caller = false
  }

  # Not overridable:
  Int read_length = 250
  Float lod_threshold = -20.0
  String cross_check_fingerprints_by = "READGROUP"
  String recalibrated_bam_basename = sample_and_unmapped_bams.base_file_name + ".aligned.duplicates_marked.recalibrated"

  call ToBam.UnmappedBamToAlignedBam {
    input:
      sample_and_unmapped_bams    = sample_and_unmapped_bams,
      references                  = references,
      papi_settings               = papi_settings,

      contamination_sites_ud = references.contamination_sites_ud,
      contamination_sites_bed = references.contamination_sites_bed,
      contamination_sites_mu = references.contamination_sites_mu,

      cross_check_fingerprints_by = cross_check_fingerprints_by,
      haplotype_database_file     = haplotype_database_file,
      lod_threshold               = lod_threshold,
      recalibrated_bam_basename   = recalibrated_bam_basename
  }

  if (provide_bam_output) {
    File provided_output_bam = UnmappedBamToAlignedBam.output_bam
    File provided_output_bam_index = UnmappedBamToAlignedBam.output_bam_index
  }

  call QC.CollectWgsMetrics as CollectWgsMetrics {
    input:
      input_bam = UnmappedBamToAlignedBam.output_bam,
      input_bam_index = UnmappedBamToAlignedBam.output_bam_index,
      metrics_filename = sample_and_unmapped_bams.base_file_name + ".wgs_metrics",
      ref_fasta = references.reference_fasta.ref_fasta,
      ref_fasta_index = references.reference_fasta.ref_fasta_index,
      wgs_coverage_interval_list = wgs_coverage_interval_list,
      read_length = read_length,
      preemptible_tries = papi_settings.agg_preemptible_tries
  }

  ## NOTE: uncomment if workflow requires .cram output
  # call ToCram.BamToCram {
  #   input:
  #     input_bam = UnmappedBamToAlignedBam.output_bam,
  #     ref_fasta = references.reference_fasta.ref_fasta,
  #     ref_fasta_index = references.reference_fasta.ref_fasta_index,
  #     ref_dict = references.reference_fasta.ref_dict,
  #     base_file_name = sample_and_unmapped_bams.base_file_name,
  #     agg_preemptible_tries = papi_settings.agg_preemptible_tries
  # }

  call Flagstat {
    input:
      input_bam = UnmappedBamToAlignedBam.output_bam,
      output_name = sample_and_unmapped_bams.base_file_name + ".flagstat.txt"
  }

  call WriteFlagstatOut {
    input: 
      input_flagstat = Flagstat.output_file,
      output_name = sample_and_unmapped_bams.base_file_name + ".writeout.flagstat.txt", 
      sample_ID = sample_and_unmapped_bams.base_file_name
  }

  call GetBamHeader {
    input: 
      input_bam = UnmappedBamToAlignedBam.output_bam, 
      output_name = sample_and_unmapped_bams.base_file_name + ".bam.header.txt"
  }

  # Outputs that will be retained when execution is complete
  output {
    File flagstat_output = WriteFlagstatOut.output_file
  }
} 


task Flagstat {
  input {
    File input_bam
    String output_name
  }

  Int disk = ceil(size(input_bam, "GB") * 2)
    
  command <<< 
    samtools flagstat ~{input_bam} > ~{output_name}
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330"
    disk: disk + " GB"
    cpu: 4
    preemptible: true
  }

  output {
    File output_file = "~{output_name}"
  }
}

task WriteFlagstatOut {
  input {
    File input_flagstat
    String output_name
    String sample_ID
  }

  command <<<
    head -1 ~{input_flagstat} > temp.txt
    VAR=$(awk -F ' ' '{print $1, $3}' temp.txt)
    echo -n ~{sample_ID} >> ~{output_name}
    echo -n " " >> ~{output_name}
    echo -n $VAR >> ~{output_name}
    rm temp.txt
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330"
    cpu: 4
    preemptible: true
  }

  output { 
    File output_file = "~{output_name}"
  }
}

task GetBamHeader {
  input {
    File input_bam
    String output_name
  }

  Int disk = ceil(size(input_bam, "GB") * 2)

  command <<<
    samtools view -H input_bam > ~{output_name}
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330"
    disk: disk + " GB"
    cpu: 4
    preemptible: true
  }

  output {
    File output_file = "~{output_name}"
  }
}