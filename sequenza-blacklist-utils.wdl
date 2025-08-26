## This WDL file calls the Sequenza Pype workflow and run_sequenza.py via a docker for Sequenza.
## The source information follows:
##
## ***************************************************************************
## run_sequenza.py: 
## 2017, Francesco Favero <favero.francesco@gmail.com>
## Center for Biological Sequence Analysis, Technical University of Denmark
## Reference Paper: https://www.sciencedirect.com/science/article/pii/S0923753419313237 
## 
## ---------------------------------------------------------------------------
## Last modified: 19 August, 2019
## ---------------------------------------------------------------------------
## ***************************************************************************
##
## ***************************************************************************
## Sequenza Pype: 
## 2017, Francesco Favero <favero.francesco@gmail.com>
## Center for Biological Sequence Analysis, Technical University of Denmark
## https://bitbucket.org/sequenzatools/sequenza_pype_modules/src/master/
## 
## ---------------------------------------------------------------------------
## Last modified: 16 November, 2020
## ---------------------------------------------------------------------------
## ***************************************************************************
##
## ***************************************************************************
## Sequenza Docker:
## Center for Biological Sequence Analysis, Technical University of Denmark
## https://hub.docker.com/r/sequenza/sequenza 
## ---------------------------------------------------------------------------
## Last pushed: 20 July, 2018
## ---------------------------------------------------------------------------
## ***************************************************************************
##
## This WDL pipeline is optimized for Cromwell on Azure. It takes in a tumour BAM
## and its matched normal BAM, to output allele-specific copy numbers. The WDL
## calls run_sequenza.py, which runs the Sequenza Pype workflow installed on the docker. 
##
## Expected inputs:
## - Tumour BAM file 
## - Matched Normal BAM file
## - Reference Fasta (can be gzipped or unzipped)
## - Reference GC in WIG format (generated from Fasta using sequenza-utils gc_wiggle)
##
## Outputs: 
## - <sample_ID>_logs.tar.gz
## - <sample_ID>_seqz_bin.tar.gz
## - <sample_ID>_sequenza.tar.gz
##
## All Sequenza source code, including this WDL file, can be found at: https://sequenzatools.bitbucket.io/#/home 

task SequenzaTask {
  
  String sampleName
  File normalBam
  File tumorBam
  File referenceFasta
  File referenceFai
  
  command {

    sequenza-utils gc_wiggle -w 50 --fasta ${referenceFasta} -o output.gc50Base.wig.gz

    sequenza-utils bam2seqz -n ${normalBam} -t ${tumorBam} --fasta ${referenceFasta} \
      -gc output.gc50Base.wig.gz -o "${sampleName}.seqz.gz"

    sequenza-utils seqz_binning --seqz "${sampleName}.seqz.gz" -w 50 -o "${sampleName}_small.seqz.gz"

    mkdir temp
    Rscript /run_sequenza.R --seqz_file "${sampleName}_small.seqz.gz" --output_dir "./temp" --sample_name ${sampleName}
    tar -czvf ${sampleName}_outputs.tar.gz ./temp/
  }

  runtime {
    docker: "apariciobioinformaticscoop/sequenza-outputs"
    memory: "32 GB" 
    cpu: 8 # 24 -> 8
    disk: "250 GB" # 500 -> 250
    maxRetries: 3
  }

  output {
    File outputs = "${sampleName}_outputs.tar.gz"
  }
}

task BlacklistFilter {
    File input_bam

    String output_file_name
    String sample_name

    File blacklist

    Float multiplier = 15
    Int disk_size = ceil(size(input_bam, "GB") * multiplier) + 150 

    command {
        bedtools intersect -v -abam ${input_bam} -b ${blacklist} > ${output_file_name}
        samtools index ${output_file_name}
    }

    runtime {
        docker: "apariciobioinformaticscoop/wasp-mapping:latest"
        disk: disk_size + " GB" 
        cpu: 40 # changed from 24
        memory: "64 GB" # changed from 24
        preemptible: true
        maxRetries: 3
    }

    output {
        File outfile = "${output_file_name}"
        File index = "${output_file_name}.bai"
    }
}

workflow SequenzaWorkflow {
  String study # metadata for clean-up automation 
  String sampleName
  String normalName # metadata for cleanup 
  File normalBam
  String tumorName # metadata for cleanup
  File tumorBam
  File referenceFasta
  File referenceFai
  File blacklist

  call BlacklistFilter {
      input:  input_bam = tumorBam,
              output_file_name = sampleName + ".postblacklist.bam",
              sample_name = sampleName,
              blacklist = blacklist
  }

  call SequenzaTask {
    input: sampleName = sampleName,
           normalBam=normalBam, 
           tumorBam=BlacklistFilter.outfile,
           referenceFasta=referenceFasta,
           referenceFai=referenceFai
  }
}
