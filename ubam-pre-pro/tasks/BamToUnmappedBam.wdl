version 1.0
## Copyright Broad Institute, 2018
## 
## This WDL converts BAM  to unmapped BAMs
##
## Requirements/expectations :
## - BAM file
##
## Outputs :
## - Sorted Unmapped BAMs
##
## Cromwell version support
## - Successfully tested on v47
## - Does not work on versions < v23 due to output syntax
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
## Updated by Aparicio Lab (BC Cancer Research Centre) May 2022 to optimize runtime 
## parameters for Cromwell on Azure implementation (instead of Google Cloud).


# WORKFLOW DEFINITION
workflow BamToUnmappedBams {
  input {
    File input_bam

    Int preemptible_tries = 2

    Int additional_disk_size = 20
    String gatk_docker = "broadinstitute/gatk:latest"
    String gatk_path = "/gatk/gatk"
  }
    
  Float input_size = size(input_bam, "GB")
    
  call RevertSam {
    input:
      input_bam = input_bam,
      disk_size = ceil(input_size * 4) + additional_disk_size,
      docker = gatk_docker,
      gatk_path = gatk_path,
      preemptible_tries = preemptible_tries
  }

  scatter (reverted_bam in RevertSam.reverted_bams) {
    String output_basename = basename(reverted_bam, ".coord.sorted.unmapped.bam")
    Float reverted_bam_size = size(reverted_bam, "GB")
  
    call SortSam {
      input:
        input_bam = reverted_bam,
        sorted_bam_name = output_basename + ".unmapped.bam",
        disk_size = ceil(reverted_bam_size * 5) + additional_disk_size,
        docker = gatk_docker,
        gatk_path = gatk_path,
        preemptible_tries = preemptible_tries
    }
  }

  output {
    Array[File] unmapped_bams = SortSam.sorted_bam
  }
}

task RevertSam {
  input {
    File input_bam
    String gatk_path

    Int disk_size
    String docker
    Int machine_mem_gb = 2 # 2 -> 20
    Int preemptible_tries
  }
    Int command_mem_gb = machine_mem_gb - 1    ####Needs to occur after machine_mem_gb is set 
    Int increased_disk_size = disk_size * 3

  command <<< 
    ~{gatk_path} --java-options "-Xmx~{command_mem_gb}g" \
    RevertSam \
    --INPUT ~{input_bam} \
    --OUTPUT ./ \
    --OUTPUT_BY_READGROUP true \
    --VALIDATION_STRINGENCY LENIENT \
    --ATTRIBUTE_TO_CLEAR FT \
    --ATTRIBUTE_TO_CLEAR CO \
    --SORT_ORDER coordinate
  >>>

  runtime {
    docker: docker
    disk: increased_disk_size + " GB"
    memory: machine_mem_gb + " GB"
    preemptible: true
    maxRetries: preemptible_tries
  }

  output {
    Array[File] reverted_bams = glob("*.bam")
  }
}

task SortSam {
  input {
    File input_bam
    String sorted_bam_name

    String gatk_path
    Int disk_size
    String docker
    Int machine_mem_gb = 15 # 4 -> 15
    Int preemptible_tries
  }
    Int command_mem_gb = machine_mem_gb - 1    ####Needs to occur after machine_mem_gb is set 
    Int increased_disk_size = disk_size * 3

  command <<<
    ~{gatk_path} --java-options "-Xmx~{command_mem_gb}g" \
    SortSam \
    --INPUT ~{input_bam} \
    --OUTPUT ~{sorted_bam_name} \
    --SORT_ORDER queryname \
    --MAX_RECORDS_IN_RAM 1000000
  >>>
  
  runtime {
    docker: docker
    disk: increased_disk_size + " GB"
    memory: machine_mem_gb + " GB"
    preemptible: true
    maxRetries: preemptible_tries
  }
  
  output {
    File sorted_bam = "~{sorted_bam_name}"
  }
}
