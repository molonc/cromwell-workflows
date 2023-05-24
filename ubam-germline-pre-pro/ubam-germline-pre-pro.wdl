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
## UPDATE NOTES :
## Updated by Jenna Liebe at the Aparicio Lab (BC Cancer Research Centre) May/June 2022.
## 
## This pipeline has been modified from its original, which can be found at 
## https://github.com/microsoft/gatk4-genome-processing-pipeline-azure. Major changes 
## include adding BAM to uBAM conversion, flagstat, and BAM to CRAM tasks; changing the default pipeline
## settings to output a VCF instead of a gVCF file; and renaming the workflow to 
## "UbamGermlinePrePro".


import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/ubam-germline-pre-pro/tasks/BamToUnmappedBam.wdl" as ToUbam
import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/ubam-germline-pre-pro/tasks/UnmappedBamToAlignedBam.wdl" as ToBam
import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/ubam-germline-pre-pro/tasks/BamToCram.wdl" as ToCram
import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/ubam-germline-pre-pro/tasks/AggregatedBamQC.wdl" as AggregatedQC
import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/ubam-germline-pre-pro/tasks/Qc.wdl" as QC
import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/ubam-germline-pre-pro/tasks/VariantCalling.wdl" as ToGvcf
import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/ubam-germline-pre-pro/tasks/GermlineStructs.wdl"


# WORKFLOW DEFINITION
struct Runtime {
    String gatk_docker
    File? gatk_override
    Int max_retries
    Int preemptible
    Int cpu
    Int machine_mem
    Int command_mem
    Int disk
    #Int boot_disk_size
}

workflow UbamGermlinePrePro {

  String pipeline_version = "1.4"

  input {
    SampleInfo sample_info
    File input_bam
    GermlineSingleSampleReferences references
    PapiSettings papi_settings
    File wgs_coverage_interval_list

    File? haplotype_database_file
    Boolean provide_bam_output = true
    Boolean use_gatk3_haplotype_caller = false

    # [START] Added inputs for funcotator
    Boolean? compress_vcfs
    Boolean? run_orientation_bias_mixture_model_filter
    Boolean? make_bamout
    Int? preemptible
    Int? max_retries
    Int preemptible_or_default = select_first([preemptible, 2])
    Int max_retries_or_default = select_first([max_retries, 2])
    Int small_task_cpu = 2 
    Int small_task_mem = 4 
    Int small_task_disk = 100 

    # Use as a last resort to increase the disk given to every task in case of ill behaving data
    Int? emergency_extra_disk

    # Funcotator inputs
    Boolean? run_funcotator
    String? sequencing_center
    String? sequence_source
    String? funco_reference_version
    String? funco_output_format
    Boolean? funco_compress
    Boolean? funco_use_gnomad_AF
    File? funco_data_sources_tar_gz
    String? funco_transcript_selection_mode
    File? funco_transcript_selection_list
    Array[String]? funco_annotation_defaults
    Array[String]? funco_annotation_overrides
    Array[String]? funcotator_excluded_fields
    Boolean? funco_filter_funcotations
    String? funcotator_extra_args

    String funco_default_output_format = "MAF"

    # runtime (from mutect2 funcotator)
    String gatk_docker
    File? gatk_override
    String basic_bash_docker = "ubuntu:16.04"
    Boolean? filter_funcotations

    # [END] Added inputs for funcotator
    # ----------------------------------------

    # These are multipliers to multipler inputs by to make sure we have enough disk to accommodate for possible output sizes
    # Large is for Bams/WGS vcfs
    # Small is for metrics/other vcfs
    Float large_input_to_output_multiplier = 2.25
    Float small_input_to_output_multiplier = 2.0
    Float cram_to_bam_multiplier = 6.0
  }

  Boolean compress = select_first([compress_vcfs, false])
  Boolean run_funcotator_or_default = select_first([run_funcotator, false])
  Boolean filter_funcotations_or_default = select_first([filter_funcotations, true])

  # If no tar is provided, the task downloads one from broads ftp server
  Int funco_tar_size = if defined(funco_data_sources_tar_gz) then ceil(size(funco_data_sources_tar_gz, "GB") * 3) else 100
  Int gatk_override_size = if defined(gatk_override) then ceil(size(gatk_override, "GB")) else 0

  # This is added to every task as padding, should increase if systematically you need more disk for every call
  Int disk_pad = 10 + gatk_override_size + select_first([emergency_extra_disk,0])

  Runtime standard_runtime = {"gatk_docker": gatk_docker, "gatk_override": gatk_override,
          "max_retries": max_retries_or_default, "preemptible": preemptible_or_default, "cpu": small_task_cpu,
          "machine_mem": small_task_mem * 1000, "command_mem": small_task_mem * 1000 - 500,
          "disk": small_task_disk + disk_pad}

  # Not overridable:
  Int read_length = 250
  Float lod_threshold = -20.0
  String cross_check_fingerprints_by = "READGROUP"
  String recalibrated_bam_basename = sample_info.base_file_name + ".aligned.duplicates_marked.recalibrated"

  call ToUbam.BamToUnmappedBams {
    input:
      input_bam             = input_bam
  }

  call ToBam.UnmappedBamToAlignedBam {
    input:
      unmapped_bams         = BamToUnmappedBams.unmapped_bams,
      sample_info           = sample_info,
      references            = references,
      papi_settings         = papi_settings,

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

   call AggregatedQC.AggregatedBamQC {
    input:
      base_recalibrated_bam = UnmappedBamToAlignedBam.output_bam,
      base_recalibrated_bam_index = UnmappedBamToAlignedBam.output_bam_index,
      base_name = sample_info.base_file_name,
      sample_name = sample_info.sample_name,
      recalibrated_bam_base_name = recalibrated_bam_basename,
      haplotype_database_file = haplotype_database_file,
      references = references,
      papi_settings = papi_settings
  }

  call ToCram.BamToCram as BamToCram {
    input:
      input_bam = UnmappedBamToAlignedBam.output_bam,
      ref_fasta = references.reference_fasta.ref_fasta,
      ref_fasta_index = references.reference_fasta.ref_fasta_index,
      ref_dict = references.reference_fasta.ref_dict,
      duplication_metrics = UnmappedBamToAlignedBam.duplicate_metrics,
      chimerism_metrics = AggregatedBamQC.agg_alignment_summary_metrics,
      base_file_name = sample_info.base_file_name,
      agg_preemptible_tries = papi_settings.agg_preemptible_tries
  }

  # QC the sample WGS metrics (stringent thresholds)
  call QC.CollectWgsMetrics as CollectWgsMetrics {
    input:
      input_bam = UnmappedBamToAlignedBam.output_bam,
      input_bam_index = UnmappedBamToAlignedBam.output_bam_index,
      metrics_filename = sample_info.base_file_name + ".wgs_metrics",
      ref_fasta = references.reference_fasta.ref_fasta,
      ref_fasta_index = references.reference_fasta.ref_fasta_index,
      wgs_coverage_interval_list = wgs_coverage_interval_list,
      read_length = read_length,
      preemptible_tries = papi_settings.agg_preemptible_tries
  }

  # QC the sample raw WGS metrics (common thresholds)
  call QC.CollectRawWgsMetrics as CollectRawWgsMetrics {
    input:
      input_bam = UnmappedBamToAlignedBam.output_bam,
      input_bam_index = UnmappedBamToAlignedBam.output_bam_index,
      metrics_filename = sample_info.base_file_name + ".raw_wgs_metrics",
      ref_fasta = references.reference_fasta.ref_fasta,
      ref_fasta_index = references.reference_fasta.ref_fasta_index,
      wgs_coverage_interval_list = wgs_coverage_interval_list,
      read_length = read_length,
      preemptible_tries = papi_settings.agg_preemptible_tries
  }

  call ToGvcf.VariantCalling as BamToGvcf {
    input:
      calling_interval_list = references.calling_interval_list,
      evaluation_interval_list = references.evaluation_interval_list,
      haplotype_scatter_count = references.haplotype_scatter_count,
      break_bands_at_multiples_of = references.break_bands_at_multiples_of,
      contamination = UnmappedBamToAlignedBam.contamination,
      input_bam = UnmappedBamToAlignedBam.output_bam,
      input_bam_index = UnmappedBamToAlignedBam.output_bam_index,
      ref_fasta = references.reference_fasta.ref_fasta,
      ref_fasta_index = references.reference_fasta.ref_fasta_index,
      ref_dict = references.reference_fasta.ref_dict,
      dbsnp_vcf = references.dbsnp_vcf,
      dbsnp_vcf_index = references.dbsnp_vcf_index,
      base_file_name = sample_info.base_file_name,
      final_vcf_base_name = sample_info.final_gvcf_base_name,
      agg_preemptible_tries = papi_settings.agg_preemptible_tries,
      use_gatk3_haplotype_caller = use_gatk3_haplotype_caller
  }

  if (run_funcotator_or_default) {
      File funcotate_vcf_input = BamToGvcf.output_vcf # unsure of this 
      File funcotate_vcf_input_index = BamToGvcf.output_vcf_index

      call Funcotate {
          input:
            ref_fasta = references.reference_fasta.ref_fasta,
            ref_fai = references.reference_fasta.ref_fasta_index,
            ref_dict = references.reference_fasta.ref_dict,
            input_vcf = funcotate_vcf_input,
            input_vcf_idx = funcotate_vcf_input_index,
            reference_version = select_first([funco_reference_version, "hg38"]),
            output_file_base_name = basename(funcotate_vcf_input, ".vcf") + ".annotated",
            output_format = if defined(funco_output_format) then "" + funco_output_format else funco_default_output_format,
            compress = if defined(funco_compress) then select_first([funco_compress]) else false,
            use_gnomad = if defined(funco_use_gnomad_AF) then select_first([funco_use_gnomad_AF]) else false,
            data_sources_tar_gz = funco_data_sources_tar_gz,
            sequencing_center = sequencing_center,
            sequence_source = sequence_source,
            transcript_selection_mode = funco_transcript_selection_mode,
            transcript_selection_list = funco_transcript_selection_list,
            annotation_defaults = funco_annotation_defaults,
            annotation_overrides = funco_annotation_overrides,
            funcotator_excluded_fields = funcotator_excluded_fields,
            filter_funcotations = filter_funcotations_or_default,
            extra_args = funcotator_extra_args,
            runtime_params = standard_runtime,
            disk_space = ceil(size(funcotate_vcf_input, "GB") * large_input_to_output_multiplier)  + funco_tar_size + disk_pad
      }

      call MafFuncotate {
          input:
              ref_fasta = references.reference_fasta.ref_fasta,
              ref_fai = references.reference_fasta.ref_fasta_index,
              ref_dict = references.reference_fasta.ref_dict,
              input_vcf = funcotate_vcf_input,
              input_vcf_idx = funcotate_vcf_input_index,
              reference_version = select_first([funco_reference_version, "hg38"]),
              output_file_base_name = basename(funcotate_vcf_input, ".vcf") + ".annotated",
              output_format = "MAF",
              compress = if defined(funco_compress) then select_first([funco_compress]) else false,
              use_gnomad = if defined(funco_use_gnomad_AF) then select_first([funco_use_gnomad_AF]) else false,
              data_sources_tar_gz = funco_data_sources_tar_gz,
              sequencing_center = sequencing_center,
              sequence_source = sequence_source,
              transcript_selection_mode = funco_transcript_selection_mode,
              transcript_selection_list = funco_transcript_selection_list,
              annotation_defaults = funco_annotation_defaults,
              annotation_overrides = funco_annotation_overrides,
              funcotator_excluded_fields = funcotator_excluded_fields,
              filter_funcotations = filter_funcotations_or_default,
              extra_args = funcotator_extra_args,
              runtime_params = standard_runtime,
              disk_space = ceil(size(funcotate_vcf_input, "GB") * large_input_to_output_multiplier)  + funco_tar_size + disk_pad
      }
  }

  call Flagstat {
    input:
      input_bam = UnmappedBamToAlignedBam.output_bam,
      output_name = sample_info.base_file_name + ".flagstat.txt"
  }

  if (provide_bam_output) {
    File provided_output_bam = UnmappedBamToAlignedBam.output_bam
    File provided_output_bam_index = UnmappedBamToAlignedBam.output_bam_index
  }


  # Outputs that will be retained when execution is complete
  output {
    Array[File] quality_yield_metrics = UnmappedBamToAlignedBam.quality_yield_metrics

    File? output_bam = provided_output_bam
    File? output_bam_index = provided_output_bam_index

    File read_group_alignment_summary_metrics = AggregatedBamQC.read_group_alignment_summary_metrics
    File read_group_gc_bias_detail_metrics = AggregatedBamQC.read_group_gc_bias_detail_metrics
    File agg_alignment_summary_metrics = AggregatedBamQC.agg_alignment_summary_metrics
    File agg_gc_bias_detail_metrics = AggregatedBamQC.agg_gc_bias_detail_metrics

    File wgs_metrics = CollectWgsMetrics.metrics
    File raw_wgs_metrics = CollectRawWgsMetrics.metrics
    File duplicate_metrics = UnmappedBamToAlignedBam.duplicate_metrics
    File output_bqsr_reports = UnmappedBamToAlignedBam.output_bqsr_reports

    File gvcf_summary_metrics = BamToGvcf.vcf_summary_metrics
    File gvcf_detail_metrics = BamToGvcf.vcf_detail_metrics
    File output_vcf = BamToGvcf.output_vcf
    File output_vcf_index = BamToGvcf.output_vcf_index

    File output_cram = BamToCram.output_cram
    File output_cram_index = BamToCram.output_cram_index
    File output_cram_md5 = BamToCram.output_cram_md5

    # Funcotator outputs 
    File? funcotated_file = Funcotate.funcotated_output_file
    File? funcotated_file_index = Funcotate.funcotated_output_file_index
    File? maf_funcotated_file = MafFuncotate.maf_funcotated_output_file
    File? maf_funcotated_file_index = MafFuncotate.maf_funcotated_output_file_index

    File output_flagstat = Flagstat.output_file
  }
}

task Funcotate {
     input {
       File ref_fasta
       File ref_fai
       File ref_dict
       File input_vcf
       File input_vcf_idx
       String reference_version
       String output_file_base_name
       String output_format
       Boolean compress
       Boolean use_gnomad
       # This should be updated when a new version of the data sources is released
       # TODO: Make this dynamically chosen in the command.
       #File? data_sources_tar_gz = "gs://broad-public-datasets/funcotator/funcotator_dataSources.v1.7.20200521s.tar.gz"
       File? data_sources_tar_gz
       String? control_id
       String? sequencing_center
       String? sequence_source
       String? transcript_selection_mode
       File? transcript_selection_list
       Array[String]? annotation_defaults
       Array[String]? annotation_overrides
       Array[String]? funcotator_excluded_fields
       Boolean? filter_funcotations
       File? interval_list

       String? extra_args
       String? gcs_project_for_requester_pays

       # ==============
       Runtime runtime_params
       Int? disk_space   #override to request more disk than default small task params

       # You may have to change the following two parameter values depending on the task requirements
       Int default_ram_mb = 3000
       # WARNING: In the workflow, you should calculate the disk space as an input to this task (disk_space_gb).  Please see [TODO: Link from Jose] for examples.
       Int default_disk_space_gb = 150
     }

     # ==============
     # Process input args:
     String output_maf = output_file_base_name + ".maf"
     String output_maf_index = output_maf + ".idx"
     String output_vcf = output_file_base_name + if compress then ".vcf.gz" else ".vcf"
     String output_vcf_idx = output_vcf +  if compress then ".tbi" else ".idx"
     String output_file = if output_format == "MAF" then output_maf else output_vcf
     String output_file_index = if output_format == "MAF" then output_maf_index else output_vcf_idx
     String transcript_selection_arg = if defined(transcript_selection_list) then " --transcript-list " else ""
     String annotation_def_arg = if defined(annotation_defaults) then " --annotation-default " else ""
     String annotation_over_arg = if defined(annotation_overrides) then " --annotation-override " else ""
     String filter_funcotations_args = if defined(filter_funcotations) && (filter_funcotations) then " --remove-filtered-variants " else ""
     String excluded_fields_args = if defined(funcotator_excluded_fields) then " --exclude-field " else ""
     String interval_list_arg = if defined(interval_list) then " -L " else ""
     String extra_args_arg = select_first([extra_args, ""])

     String dollar = "$"

     parameter_meta{
      ref_fasta: {localization_optional: true}
      ref_fai: {localization_optional: true}
      ref_dict: {localization_optional: true}
      input_vcf: {localization_optional: true}
      input_vcf_idx: {localization_optional: true}
     }

     command <<<
         set -e
         export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}
         # Extract our data sources:
         echo "Extracting data sources zip file..."
         mkdir datasources_dir
         tar zxvf ~{data_sources_tar_gz} -C datasources_dir --strip-components 1
         DATA_SOURCES_FOLDER="$PWD/datasources_dir"

         # Handle gnomAD:
         if ~{use_gnomad} ; then
             echo "Enabling gnomAD..."
             for potential_gnomad_gz in gnomAD_exome.tar.gz gnomAD_genome.tar.gz ; do
                 if [[ -f ~{dollar}{DATA_SOURCES_FOLDER}/~{dollar}{potential_gnomad_gz} ]] ; then
                     cd ~{dollar}{DATA_SOURCES_FOLDER}
                     tar -zvxf ~{dollar}{potential_gnomad_gz}
                     cd -
                 else
                     echo "ERROR: Cannot find gnomAD folder: ~{dollar}{potential_gnomad_gz}" 1>&2
                     false
                 fi
             done
         fi
         # Run Funcotator:
         gatk --java-options "-Xmx~{runtime_params.command_mem}m" Funcotator \
             --data-sources-path $DATA_SOURCES_FOLDER \
             --ref-version ~{reference_version} \
             --output-file-format ~{output_format} \
             -R ~{ref_fasta} \
             -V ~{input_vcf} \
             -O ~{output_file} \
             ~{interval_list_arg} ~{default="" interval_list} \
             --annotation-default normal_barcode:~{default="Unknown" control_id} \
             --annotation-default Center:~{default="Unknown" sequencing_center} \
             --annotation-default source:~{default="Unknown" sequence_source} \
             ~{"--transcript-selection-mode " + transcript_selection_mode} \
             ~{transcript_selection_arg}~{default="" sep=" --transcript-list " transcript_selection_list} \
             ~{annotation_def_arg}~{default="" sep=" --annotation-default " annotation_defaults} \
             ~{annotation_over_arg}~{default="" sep=" --annotation-override " annotation_overrides} \
             ~{excluded_fields_args}~{default="" sep=" --exclude-field " funcotator_excluded_fields} \
             ~{filter_funcotations_args} \
             ~{extra_args_arg} \
             ~{"--gcs-project-for-requester-pays " + gcs_project_for_requester_pays}
         # Make sure we have a placeholder index for MAF files so this workflow doesn't fail:
         if [[ "~{output_format}" == "MAF" ]] ; then
            touch ~{output_maf_index}
         fi
     >>>

    runtime {
        docker: runtime_params.gatk_docker
        #bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        disk: select_first([disk_space, runtime_params.disk]) + " GB"
        preemptible: true
        maxRetries: runtime_params.max_retries
        cpu: 6
    }

    output {
        File funcotated_output_file = "~{output_file}"
        File funcotated_output_file_index = "~{output_file_index}"
    }
}


task MafFuncotate {
     input {
       File ref_fasta
       File ref_fai
       File ref_dict
       File input_vcf
       File input_vcf_idx
       String reference_version
       String output_file_base_name
       String output_format
       Boolean compress
       Boolean use_gnomad
       File? data_sources_tar_gz
       String? control_id
       String? sequencing_center
       String? sequence_source
       String? transcript_selection_mode
       File? transcript_selection_list
       Array[String]? annotation_defaults
       Array[String]? annotation_overrides
       Array[String]? funcotator_excluded_fields
       Boolean? filter_funcotations
       File? interval_list

       String? extra_args
       String? gcs_project_for_requester_pays

       # ==============
       Runtime runtime_params
       Int? disk_space   #override to request more disk than default small task params

       # You may have to change the following two parameter values depending on the task requirements
       Int default_ram_mb = 3000
       Int default_disk_space_gb = 150
     }

     # ==============
     # Process input args:
     String output_maf = output_file_base_name + ".maf"
     String output_maf_index = output_maf + ".idx"
     String output_vcf = output_file_base_name + if compress then ".vcf.gz" else ".vcf"
     String output_vcf_idx = output_vcf +  if compress then ".tbi" else ".idx"
     String output_file = if output_format == "MAF" then output_maf else output_vcf
     String output_file_index = if output_format == "MAF" then output_maf_index else output_vcf_idx
     String transcript_selection_arg = if defined(transcript_selection_list) then " --transcript-list " else ""
     String annotation_def_arg = if defined(annotation_defaults) then " --annotation-default " else ""
     String annotation_over_arg = if defined(annotation_overrides) then " --annotation-override " else ""
     String filter_funcotations_args = if defined(filter_funcotations) && (filter_funcotations) then " --remove-filtered-variants " else ""
     String excluded_fields_args = if defined(funcotator_excluded_fields) then " --exclude-field " else ""
     String interval_list_arg = if defined(interval_list) then " -L " else ""
     String extra_args_arg = select_first([extra_args, ""])

     String dollar = "$"

     parameter_meta{
      ref_fasta: {localization_optional: true}
      ref_fai: {localization_optional: true}
      ref_dict: {localization_optional: true}
      input_vcf: {localization_optional: true}
      input_vcf_idx: {localization_optional: true}
     }

     command <<<
         set -e
         export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}
         # Extract our data sources:
         echo "Extracting data sources zip file..."
         mkdir datasources_dir
         tar zxvf ~{data_sources_tar_gz} -C datasources_dir --strip-components 1
         DATA_SOURCES_FOLDER="$PWD/datasources_dir"
         # Handle gnomAD:
         if ~{use_gnomad} ; then
             echo "Enabling gnomAD..."
             for potential_gnomad_gz in gnomAD_exome.tar.gz gnomAD_genome.tar.gz ; do
                 if [[ -f ~{dollar}{DATA_SOURCES_FOLDER}/~{dollar}{potential_gnomad_gz} ]] ; then
                     cd ~{dollar}{DATA_SOURCES_FOLDER}
                     tar -zvxf ~{dollar}{potential_gnomad_gz}
                     cd -
                 else
                     echo "ERROR: Cannot find gnomAD folder: ~{dollar}{potential_gnomad_gz}" 1>&2
                     false
                 fi
             done
         fi
         # Run Funcotator:
         gatk --java-options "-Xmx~{runtime_params.command_mem}m" Funcotator \
             --data-sources-path $DATA_SOURCES_FOLDER \
             --ref-version ~{reference_version} \
             --output-file-format ~{output_format} \
             -R ~{ref_fasta} \
             -V ~{input_vcf} \
             -O ~{output_file} \
             ~{interval_list_arg} ~{default="" interval_list} \
             --annotation-default normal_barcode:~{default="Unknown" control_id} \
             --annotation-default Center:~{default="Unknown" sequencing_center} \
             --annotation-default source:~{default="Unknown" sequence_source} \
             ~{"--transcript-selection-mode " + transcript_selection_mode} \
             ~{transcript_selection_arg}~{default="" sep=" --transcript-list " transcript_selection_list} \
             ~{annotation_def_arg}~{default="" sep=" --annotation-default " annotation_defaults} \
             ~{annotation_over_arg}~{default="" sep=" --annotation-override " annotation_overrides} \
             ~{excluded_fields_args}~{default="" sep=" --exclude-field " funcotator_excluded_fields} \
             ~{filter_funcotations_args} \
             ~{extra_args_arg} \
             ~{"--gcs-project-for-requester-pays " + gcs_project_for_requester_pays}         
         fi
     >>>

    runtime {
        docker: runtime_params.gatk_docker
        memory: runtime_params.machine_mem + " MB"
        disk: select_first([disk_space, runtime_params.disk]) + " GB"
        preemptible: true
        maxRetries: runtime_params.max_retries
        cpu: 6
    }

    output {
        File maf_funcotated_output_file = "~{output_file}"
        File maf_funcotated_output_file_index = "~{output_file_index}"
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
