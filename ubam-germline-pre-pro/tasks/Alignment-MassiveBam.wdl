version 1.0

## Copyright Broad Institute, 2018
##
## This WDL defines tasks used for alignment of human whole-genome or exome sequencing data.
##
## Runtime parameters are often optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.


import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/ubam-germline-pre-pro/tasks/GermlineStructs.wdl"


# Get version of BWA
task GetBwaVersion {
  command {
    # not setting set -o pipefail here because /bwa has a rc=1 and we dont want to allow rc=1 to succeed because
    # the sed may also fail with that error and that is something we actually want to fail on.
    /usr/gitc/bwa 2>&1 | \
    grep -e '^Version' | \
    sed 's/Version: //'
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330"
    memory: "5 GB" # 1 -> 5
  }
  output {
    String bwa_version = read_string(stdout())
  }
}

# Read unmapped BAM, convert on-the-fly to FASTQ and stream to BWA MEM for alignment, then stream to MergeBamAlignment
task SamToFastqAndBwaMemAndMba {
  input {
    File input_bam
    String bwa_commandline
    String bwa_version
    String output_bam_basename

    # reference_fasta.ref_alt is the .alt file from bwa-kit
    # (https://github.com/lh3/bwa/tree/master/bwakit),
    # listing the reference contigs that are "alternative".
    ReferenceFasta reference_fasta

    Int compression_level
    Int preemptible_tries
  }

  Float unmapped_bam_size = size(input_bam, "GB")
  Float ref_size = size(reference_fasta.ref_fasta, "GB") + size(reference_fasta.ref_fasta_index, "GB") + size(reference_fasta.ref_dict, "GB")
  Float bwa_ref_size = ref_size + size(reference_fasta.ref_alt, "GB") + size(reference_fasta.ref_amb, "GB") + size(reference_fasta.ref_ann, "GB") + size(reference_fasta.ref_bwt, "GB") + size(reference_fasta.ref_pac, "GB") + size(reference_fasta.ref_sa, "GB")
  # Sometimes the output is larger than the input, or a task can spill to disk.
  # In these cases we need to account for the input (1) and the output (1.5) or the input(1), the output(1), and spillage (.5).
  Float disk_multiplier = 2.5
  Int disk_size = ceil(unmapped_bam_size + bwa_ref_size + (disk_multiplier * unmapped_bam_size) + 60) # + 20 -> + 60

  command <<<
    set -o pipefail
    set -e

    # set the bash variable needed for the command-line
    bash_ref_fasta=~{reference_fasta.ref_fasta}
    # if reference_fasta.ref_alt has data in it,
    if [ -s ~{reference_fasta.ref_alt} ]; then
      java -Xms1000m -Xmx1000m -jar /usr/gitc/picard.jar \
        SamToFastq \
        INPUT=~{input_bam} \
        FASTQ=/dev/stdout \
        INTERLEAVE=true \
        NON_PF=true | \
      /usr/gitc/~{bwa_commandline} /dev/stdin - 2> >(tee ~{output_bam_basename}.bwa.stderr.log >&2) | \
      java -Dsamjdk.compression_level=~{compression_level} -Xms1000m -Xmx1000m -jar /usr/gitc/picard.jar \
        MergeBamAlignment \
        VALIDATION_STRINGENCY=SILENT \
        EXPECTED_ORIENTATIONS=FR \
        ATTRIBUTES_TO_RETAIN=X0 \
        ATTRIBUTES_TO_REMOVE=NM \
        ATTRIBUTES_TO_REMOVE=MD \
        ALIGNED_BAM=/dev/stdin \
        UNMAPPED_BAM=~{input_bam} \
        OUTPUT=~{output_bam_basename}.bam \
        REFERENCE_SEQUENCE=~{reference_fasta.ref_fasta} \
        PAIRED_RUN=true \
        SORT_ORDER="unsorted" \
        IS_BISULFITE_SEQUENCE=false \
        ALIGNED_READS_ONLY=false \
        CLIP_ADAPTERS=false \
        MAX_RECORDS_IN_RAM=2000000 \
        ADD_MATE_CIGAR=true \
        MAX_INSERTIONS_OR_DELETIONS=-1 \
        PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
        PROGRAM_RECORD_ID="bwamem" \
        PROGRAM_GROUP_VERSION="~{bwa_version}" \
        PROGRAM_GROUP_COMMAND_LINE="~{bwa_commandline}" \
        PROGRAM_GROUP_NAME="bwamem" \
        UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
        ALIGNER_PROPER_PAIR_FLAGS=true \
        UNMAP_CONTAMINANT_READS=true \
        ADD_PG_TAG_TO_READS=false

      grep -m1 "read .* ALT contigs" ~{output_bam_basename}.bwa.stderr.log | \
      grep -v "read 0 ALT contigs"

    # else reference_fasta.ref_alt is empty or could not be found
    else
      exit 1;
    fi
  >>>
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330"
    preemptible: true
    maxRetries: preemptible_tries
    memory: "30 GB" # 14 GB -> 20 GB -> 30 GB
    cpu: "24" # 16 -> 20 -> 24
    disk: disk_size + " GB"
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File bwa_stderr_log = "~{output_bam_basename}.bwa.stderr.log"
  }
}

task SamSplitter {
  input {
    File input_bam
    Int n_reads
    Int preemptible_tries
    Int compression_level
  }

  Float unmapped_bam_size = size(input_bam, "GB")
  # Since the output bams are less compressed than the input bam we need a disk multiplier that's larger than 2.
  Float disk_multiplier = 2.5
  Int disk_size = ceil(disk_multiplier * unmapped_bam_size + 20)
  Int increased_disk_size = disk_size * 5

  command {
    set -e
    mkdir output_dir

    total_reads=$(samtools view -c ~{input_bam})

    java -Dsamjdk.compression_level=~{compression_level} -Xms3000m -jar /usr/gitc/picard.jar SplitSamByNumberOfReads \
      INPUT=~{input_bam} \
      OUTPUT=output_dir \
      SPLIT_TO_N_READS=~{n_reads} \
      TOTAL_READS_IN_INPUT=$total_reads
  }
  output {
    Array[File] split_bams = glob("output_dir/*.bam")
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330"
    preemptible: true
    maxRetries: preemptible_tries
    memory: "28 GB" # 3.75 -> 7.5 -> 14 -> 28
    disk: increased_disk_size + " GB" # multiplied disk size by 2, then 3, now 5
  }
}
