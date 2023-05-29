version 1.0
## This WDL file calls runIchorCNA.R via a docker for ichorCNA.
## The source information follows:
##
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
## one or more WIG files at a time and runs runIchorCNA.R for copy number analysis. 
##
## Expected inputs:
## - Sample WIG file(s)
##
## Outputs (for each sample WIG specified in inputs.json list)
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
##
## The pipeline was developed by Kelly Zhang at the Aparicio Lab (BC Cancer
## Research Centre) in May 2023. 

struct IchorNormal {
	String normal_id
	File normal_wig
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
		IchorNormal normal_inputs
		Array[Pair[String, File]] tumor_inputs
		IchorReferences ichor_references
		IchorInputInfo ichor_inputs

		String read_docker = "apariciobioinformaticscoop/ichorcna-updated-packages:latest"
	}

	meta {allowNestedInputs: true}

	scatter(input_wig in tumor_inputs) {
		String sample_id = input_wig.left
		String tumor_wig = input_wig.right 

		call ichorCNA {
			input: 
				library_ID = sample_id, 
				tumor_wig = tumor_wig,
				normal_wig = normal_inputs.normal_wig, 
				plot_file_type = ichor_inputs.plot_file_type, 
				gc_wig = ichor_references.gc_wig, 
				map_wig = ichor_references.map_wig, 
				centromere = ichor_references.centromere,
				input_r_script = ichor_references.input_r_script, 
				docker = read_docker 
		}
	}

	output {
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
		File genome_wide = "~{library_ID}_genomeWide.pdf"
		File genome_wide_all_sols = "~{library_ID}_genomeWide_all_sols.pdf"
		File tpdf = "~{library_ID}_tpdf.pdf"
		File bias = "~{library_ID}_bias.pdf"
		File correct = "~{library_ID}_correct.pdf"

		File corrected_depth = "~{library_ID}.correctedDepth.txt"
		File cna = "~{library_ID}.cna.seg"
		File params = "~{library_ID}.params.txt"
		File segTxt = "~{library_ID}.seg.txt"
		File seg = "~{library_ID}.seg"
		File rdata = "~{library_ID}.RData" 
	}
}

