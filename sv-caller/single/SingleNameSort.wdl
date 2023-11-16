version 1.0

import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/ubam-pre-pro/tasks/GermlineStructs.wdl"

workflow SingleNameSort {
	input {
		String tumor_name
		File tumor_bam
		File tumor_bai

        ReferenceFasta references
	}

    File reference_fasta = references.ref_fasta

	meta {allowNestedInputs: true}

	call SingleSort {
		input:
			tumor_bam = tumor_bam,
			tumor_bai = tumor_bai,
			output_tumor_name = tumor_name + ".sort.bam",
	}

	output {
		File out_tumor_bam = SingleSort.out_bam_tumor 
	}
}

task SingleSort {
	input {
		File tumor_bam
		File tumor_bai
		String output_tumor_name
	}

	Int disk_size = ceil(size(tumor_bam, "GB") * 4)

	command <<<
		samtools sort -n -o ~{output_tumor_name} ~{tumor_bam}
	>>>

	runtime {
        docker: "apariciobioinformaticscoop/wasp-mapping:latest"
        disk: disk_size + " GB"
        cpu: 8
        memory: "32 GB"
        preemptible: true
        maxRetries: 0
    }

	output {
        File out_bam_tumor =  "~{output_tumor_name}"
	}
}