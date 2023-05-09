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
  File ReferenceFastaGz
  File ReferenceGcWig
  command {
  sequenza-pipeline \
      --sample-id ${sampleName} \
      --normal-bam ${normalBam} \
      --tumor-bam ${tumorBam} \
      --reference-gz ${ReferenceFastaGz} \
      --gc_wig ${ReferenceGcWig} 
  }

  runtime {
    docker: "sequenza/sequenza"
    memory: "8 GB" # 24 -> 8
    cpu: 5 # 10 -> 5
    disk: "200 GB" 
    }

  output {
    File logs = "${sampleName}_logs.tar.gz"
    File binSeqz = "${sampleName}_seqz_bin.tar.gz"
    # Removed as output, these are partial sequenza outputs that aren't necessary; 
    # would cause workflow to fail b/c failed to output
    # File fullSeqz = "${sampleName}_parts_seqz.tar.gz" 
    File scnaRes = "${sampleName}_sequenza.tar.gz"
  }
}

workflow SequenzaWorkflow {
  File normalBam
  File tumorBam
  File ReferenceFastaGz
  File ReferenceGcWig

  call SequenzaTask {
    input: normalBam=normalBam, tumorBam=tumorBam,
           ReferenceFastaGz=ReferenceFastaGz,
           ReferenceGcWig=ReferenceGcWig
  }
}
