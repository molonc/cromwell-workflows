version 1.0
## This WDL file calls readCounter.cpp and runIchorCNA.R via a docker for ichorCNA.
## The source information follows:
##
## ***************************************************************************
## readCounter.cpp 
## (c) 2011 Daniel Lai <jujubix@cs.ubc.ca>, Gavin Ha <gha@bccrc.ca>
## Shah Lab, Department of Molecular Oncology, BC Cancer Research Agency
## All rights reserved.
## ---------------------------------------------------------------------------
## Last modified: 14 December, 2012 
## ---------------------------------------------------------------------------
## ***************************************************************************
##
## ***************************************************************************
## ichorCNA Docker:
## Aparicio Lab, Vancouver, B.C.
## https://hub.docker.com/r/apariciobioinformaticscoop/ichorcna-updated-packages
## ---------------------------------------------------------------------------
## Last pushed: 24 March, 2023
## ---------------------------------------------------------------------------
## ***************************************************************************
##
## This WDL pipeline is optimized for Cromwell on Azure. It is able to take in 
## one or more samples at a time and first runs readCounter.cpp on them to 
## produce WIG output files for each sample. With each of the outputted WIG 
## files, it will use runIchorCNA.R to estimate the tumor fraction in cell-free 
## DNA from ULP-WGS.
##
## Expected inputs:
## - Sample BAM file(s), with the associated index file (must have suffix ".bam.bai")
## and sample ID
##
## Outputs (for each sample BAM specified in inputs.json list)
## - <sample_ID>.wig
## - <sample_ID>.cna.seg
## - <sample_ID>.correctedDepth.txt
## - <sample_ID>.params.txt
## - <sample_ID>.RData
## - <sample_ID>.seg
## - <sample_ID>.seg.txt
## - <sample_ID>/<sample_ID>_bias.pdf
## - <sample_ID>/<sample_ID>_correct.pdf
## - <sample_ID>/<sample_ID>_genomeWide.pdf
## - <sample_ID>/<sample_ID>_genomeWide_all_sols.pdf
## - <sample_ID>/<sample_ID>_tpdf.pdf
##
## This pipeline has been built upon read-counter-multi.wdl, which can be found at 
## https://github.com/jliebe-bccrc/cromwell-workflows/tree/main/samtools-and-oneshots/read-counter. 
## ichorCNA has been integrated into the same wdl to run both ReadCounter and ichorCNA in succession. 
##
## The pipeline was developed by Kelly Zhang and Jenna Liebe at the Aparicio Lab (BC Cancer
## Research Centre) in April 2023. 

struct ReadCounterNormal {
	String normal_id
	File normal_bam
	File normal_index
}

struct TumorBamIndex {
	File bam
	File index
}

struct IchorReferences {
  File gc_wig 
  File map_wig
  File centromere
  File input_r_script
}

struct IchorInputInfo {
  String plot_file_type
}

# WORKFLOW DEFINITION
workflow ReadCounterIchor {
	input {
		ReadCounterNormal rc_normal_inputs
		Array[Pair[String, TumorBamIndex]] rc_tumor_inputs
		IchorReferences ichor_references
		IchorInputInfo ichor_inputs

		String read_docker = "apariciobioinformaticscoop/ichorcna-updated-packages:latest"
	}

	meta {allowNestedInputs: true}

	call ReadCounter as normalRc {
		input: 
			input_bam = rc_normal_inputs.normal_bam, 
			input_index = rc_normal_inputs.normal_index,
			sample_ID = rc_normal_inputs.normal_id,  
			docker = read_docker
	}

	scatter(input_bam_pair in rc_tumor_inputs) {
		String sample_id = input_bam_pair.left
		TumorBamIndex bam_index = input_bam_pair.right 

		call ReadCounter as tumorRc {
			input: 
				input_bam = bam_index.bam, 
				input_index = bam_index.index, 
				sample_ID = sample_id, 
				docker = read_docker
		}

		call ichorCNA {
			input: 
				library_ID = sample_id, 
				tumor_wig = tumorRc.output_wig,
				normal_wig = normalRc.output_wig, 
				plot_file_type = ichor_inputs.plot_file_type, 
				gc_wig = ichor_references.gc_wig, 
				map_wig = ichor_references.map_wig, 
				centromere = ichor_references.centromere,
				input_r_script = ichor_references.input_r_script, 
				docker = read_docker 
		}
	}

	output {
		File normal_wig = normalRc.output_wig
		Array[File] tumor_wig = tumorRc.output_wig 

		Array[File] bias = ichorCNA.bias
		Array[File] correct = ichorCNA.correct
		Array[File] genome_wide = ichorCNA.genome_wide
		Array[File] genome_wide_all_sols = ichorCNA.genome_wide_all_sols
		Array[File] tpdf = ichorCNA.tpdf
		
		Array[File] corrected_depth = ichorCNA.corrected_depth
		Array[File] cna = ichorCNA.cna
		Array[File] segTxt = ichorCNA.segTxt
		Array[File] seg = ichorCNA.seg
		Array[File] rdata = ichorCNA.rdata
	}
}

task ReadCounter {
	input {
		File input_bam
		File input_index
		String sample_ID
		String docker
	}

	Int disk_size = 200

	command <<<
		readCounter -w 1000000 -q 20 -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY ~{input_bam} > ~{sample_ID}.wig
	>>>

	runtime {
    	docker: docker
    	disk: disk_size + " GB"
    	cpu: 4
    	preemptible: true
	}

	output {
		File output_wig = "~{sample_ID}.wig"
	}
}

task ichorCNA {
	input {
		File tumor_wig
		File normal_wig
		File gc_wig
		File map_wig 
		File centromere
		File input_r_script 
		String plot_file_type
		String library_ID
		String docker
	}

	Int disk_size = 400

	command <<<
		Rscript ~{input_r_script} \
		--WIG ~{tumor_wig} \
		--NORMWIG ~{normal_wig} \
		--id ~{library_ID} \
		--gcWig ~{gc_wig} \
		--mapWig ~{map_wig} \
		--centromere ~{centromere} \
		--normal "c(0.5,0.6,0.7,0.8,0.9)" \
		--maxCN "5" \
		--scStates "c(1,3)" \
		--txnE "0.9999" \
		--txnStrength "10000"
	>>>

	runtime {
    	docker: docker
    	disk: disk_size + " GB"
		memory: "10 GB"
    	cpu: 6
    	preemptible: true
	}

	output {
		File genome_wide = "~{library_ID}/~{library_ID}_genomeWide.pdf"
		File genome_wide_all_sols = "~{library_ID}/~{library_ID}_genomeWide_all_sols.pdf"
		File tpdf = "~{library_ID}/~{library_ID}_tpdf.pdf"
		File bias = "~{library_ID}/~{library_ID}_bias.pdf"
		File correct = "~{library_ID}/~{library_ID}_correct.pdf"

		File corrected_depth = "~{library_ID}.correctedDepth.txt"
		File cna = "~{library_ID}.cna.seg"
		File params = "~{library_ID}.params.txt"
		File segTxt = "~{library_ID}.seg.txt"
		File seg = "~{library_ID}.seg"
		File rdata = "~{library_ID}.RData" 
	}
}

