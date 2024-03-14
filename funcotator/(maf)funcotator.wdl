version 1.0


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

workflow Funcotator {
    input {
      File? intervals
      File ref_fasta
      File ref_fai
      File ref_dict
      File variant_vcf_to_funcotate
      File variant_vcf_to_funcotate_index
      String output_basename
      File? pon
      File? pon_idx
      File? gnomad
      File? gnomad_idx
      File? variants_for_contamination
      File? variants_for_contamination_idx
      File? realignment_index_bundle
      String? realignment_extra_args
      Boolean? run_orientation_bias_mixture_model_filter
      String? m2_extra_args
      String? m2_extra_filtering_args
      String? getpileupsummaries_extra_args
      String? split_intervals_extra_args
      Boolean? make_bamout
      Boolean? compress_vcfs
      File? gga_vcf
      File? gga_vcf_idx
      String? gcs_project_for_requester_pays

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

      # runtime
      String gatk_docker
      File? gatk_override
      String basic_bash_docker = "ubuntu:16.04"
      Boolean? filter_funcotations

      Int? preemptible
      Int? max_retries
      Int small_task_cpu = 2
      Int small_task_mem = 4
      Int small_task_disk = 100
      #Int boot_disk_size = 12
      Int learn_read_orientation_mem = 8000
      Int filter_alignment_artifacts_mem = 9000

      # Use as a last resort to increase the disk given to every task in case of ill behaving data
      Int? emergency_extra_disk

      # These are multipliers to multipler inputs by to make sure we have enough disk to accommodate for possible output sizes
      # Large is for Bams/WGS vcfs
      # Small is for metrics/other vcfs
      Float large_input_to_output_multiplier = 2.25
      Float small_input_to_output_multiplier = 2.0
      Float cram_to_bam_multiplier = 6.0
    }

    Int preemptible_or_default = select_first([preemptible, 2])
    Int max_retries_or_default = select_first([max_retries, 2])

    Boolean compress = select_first([compress_vcfs, false])
    Boolean run_ob_filter = select_first([run_orientation_bias_mixture_model_filter, false])
    Boolean make_bamout_or_default = select_first([make_bamout, false])
    Boolean run_funcotator_or_default = select_first([run_funcotator, false])
    Boolean filter_funcotations_or_default = select_first([filter_funcotations, true])


    # If no tar is provided, the task downloads one from broads ftp server
    Int funco_tar_size = if defined(funco_data_sources_tar_gz) then ceil(size(funco_data_sources_tar_gz, "GB") * 3) else 100
    Int gatk_override_size = if defined(gatk_override) then ceil(size(gatk_override, "GB")) else 0

    # This is added to every task as padding, should increase if systematically you need more disk for every call
    Int disk_pad = 10 + gatk_override_size + select_first([emergency_extra_disk,0])

    # logic about output file names -- these are the names *without* .vcf extensions
    String unfiltered_name = output_basename + "-unfiltered"
    String filtered_name = output_basename + "-filtered"
    String funcotated_name = output_basename + "-funcotated"

    String output_vcf_name = output_basename + ".vcf"

    Runtime standard_runtime = {"gatk_docker": gatk_docker, "gatk_override": gatk_override,
            "max_retries": max_retries_or_default, "preemptible": preemptible_or_default, "cpu": small_task_cpu,
            "machine_mem": small_task_mem * 1000, "command_mem": small_task_mem * 1000 - 500,
            "disk": small_task_disk + disk_pad}



    call Funcotate {
        input:
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            input_vcf = variant_vcf_to_funcotate, 
            input_vcf_idx = variant_vcf_to_funcotate_index,
            reference_version = select_first([funco_reference_version, "hg38"]),
            output_file_base_name = basename(variant_vcf_to_funcotate, ".vcf") + ".annotated",
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
            disk_space = ceil(size(variant_vcf_to_funcotate, "GB") * large_input_to_output_multiplier)  + funco_tar_size + disk_pad
    }

    call MafFuncotate {
        input:
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            input_vcf = variant_vcf_to_funcotate, 
            input_vcf_idx = variant_vcf_to_funcotate_index,
            reference_version = select_first([funco_reference_version, "hg38"]),
            output_file_base_name = basename(variant_vcf_to_funcotate, ".vcf") + ".annotated",
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
            disk_space = ceil(size(variant_vcf_to_funcotate, "GB") * large_input_to_output_multiplier)  + funco_tar_size + disk_pad
    }

    output {
        File? funcotated_file = Funcotate.funcotated_output_file
        File? funcotated_file_index = Funcotate.funcotated_output_file_index
        File? maf_funcotated_file = MafFuncotate.maf_funcotated_output_file
        File? maf_funcotated_file_index = MafFuncotate.maf_funcotated_output_file_index
        #File? bamout = MergeBamOuts.merged_bam_out
        #File? bamout_index = MergeBamOuts.merged_bam_out_index
        #File? maf_segments = CalculateContamination.maf_segments
        #File? read_orientation_model_params = LearnReadOrientationModel.artifact_prior_table
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
       String? case_id
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
             --annotation-default tumor_barcode:~{default="Unknown" case_id} \
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
       String? case_id
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
             --annotation-default tumor_barcode:~{default="Unknown" case_id} \
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
