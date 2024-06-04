version 1.0

# Mark duplicate reads to avoid counting non-independent observations
workflow MarkDuplicates {
    input {
        File input_bam
        String sample_name
    }

    call MarkDuplicates {
        input:
            input_bam = input_bam,
            sample_name = sample_name
    }

    call SortBam {
        input:
            marked_bam = MarkDuplicates.output_bam,
            output_file_name = sample_name + ".no_duplicates.bam"
    }

    call IndexBam {
        input: 
            sorted_bam = SortBam.sorted_bam,
            output_file_name = sample_name + ".no_duplicates.bai"
    }

    output {
        File outfile = SortBam.sorted_bam
        File index = IndexBam.bai
    }
}

task MarkDuplicates {
  input {
    File input_bam
    String sample_name
    Int compression_level = 2
    Int preemptible_tries = 3

    # The program default for READ_NAME_REGEX is appropriate in nearly every case.
    # Sometimes we wish to supply "null" in order to turn off optical duplicate detection
    # This can be desirable if you don't mind the estimated library size being wrong and optical duplicate detection is taking >7 days and failing
    String? read_name_regex
    Int memory_multiplier = 3
    Int additional_disk = 20							
  }

  # The merged bam will be smaller than the sum of the parts so we need to account for the unmerged inputs and the merged output.
  # Mark Duplicates takes in as input readgroup bams and outputs a slightly smaller aggregated bam. Giving .25 as wiggleroom
  Float md_disk_multiplier = 3
  Int disk_size = ceil(md_disk_multiplier * size(input_bam, "GB")) + additional_disk

  Float memory_size = 7.5 * memory_multiplier
  Int java_memory_size = (ceil(memory_size) - 2)

  # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly
  # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
  # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"

  command <<<
    java -Dsamjdk.compression_level=~{compression_level} -Xms~{java_memory_size}g -jar /usr/gitc/picard.jar \
      MarkDuplicates \
      INPUT=~{sep=' INPUT=' input_bam} \
      OUTPUT=~{sample_name + "_unsorted"}.bam \
      METRICS_FILE=~{sample_name} \
      VALIDATION_STRINGENCY=SILENT \
      ~{"READ_NAME_REGEX=" + read_name_regex} \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      ASSUME_SORT_ORDER="queryname" \
      CLEAR_DT="false" \
      ADD_PG_TAG_TO_READS=false \
      REMOVE_DUPLICATES=true
  >>>
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330"
    preemptible: true
    maxRetries: preemptible_tries
    memory: "~{memory_size} GB"
    disk: disk_size + " GB"
  }
  output {
    File output_bam = "~{sample_name + "_unsorted"}.bam"
  }
}

task SortBam {
    input {
        File marked_bam
        String output_file_name
    }

    # based these disk numbers on GATK pre-processing tasks
    Float sort_sam_disk_multiplier = 3.25
    Int disk_size = ceil(sort_sam_disk_multiplier * size(marked_bam, "GB")) + 20

    command <<<
        samtools sort -o ~{output_file_name} ~{marked_bam} -@ 4
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/wasp-mapping:latest"
        disk: disk_size + " GB"
        cpu: 8
        memory: "10 GB"
        preemptible: true
        maxRetries: 3
    }

    output {
        File sorted_bam = "~{output_file_name}" # keep this
    }
}

task IndexBam {
    input {
        File sorted_bam
        String output_file_name
    }

    Float multiplier = 2.5
    Int disk_size = ceil(size(sorted_bam, "GB") * multiplier) + 20 

    command <<<
        samtools index -b ~{sorted_bam} ~{output_file_name}
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/wasp-mapping:latest"
        disk: disk_size + " GB"
        cpu: 12
        memory: "15 GB"
        preemptible: true
        maxRetries: 3
    }

    output {
        File bai = "~{output_file_name}" # keep this
    }
}