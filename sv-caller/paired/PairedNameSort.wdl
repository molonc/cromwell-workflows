version 1.0

import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/ubam-pre-pro/tasks/GermlineStructs.wdl"

workflow PairedNameSort {
	input {
		String tumor_name
		File tumor_bam
		File tumor_bai

		String normal_name
		File normal_bam
		File normal_bai

        ReferenceFasta references
	}

    File reference_fasta = references.ref_fasta

	meta {allowNestedInputs: true}

	call PairedSortOne {
		input:
			tumor_bam = tumor_bam,
			tumor_bai = tumor_bai,
			output_tumor_name = tumor_name + ".sort.bam"
	}

	call PairedSortTwo {
		input:
			normal_bam = normal_bam,
			normal_bai = normal_bai,
			output_normal_name = normal_name + ".sort.bam"
	}

	output {
		File out_tumor_bam = PairedSortOne.out_bam_tumor 
		File out_normal_bam = PairedSortTwo.out_bam_normal
	}
}

task PairedSortOne {
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
        cpu: 16
        memory: "64 GB"
        preemptible: true
        maxRetries: 0
    }

	output {
        File out_bam_tumor =  "~{output_tumor_name}"
	}
}

task PairedSortTwo {
	input {
		File normal_bam
		File normal_bai
		String output_normal_name
	}

	Int disk_size = ceil(size(normal_bam, "GB") * 4)

	command <<<
		samtools sort -n -o ~{output_normal_name} ~{normal_bam}
	>>>

	runtime {
        docker: "apariciobioinformaticscoop/wasp-mapping:latest"
        disk: disk_size + " GB"
        cpu: 16
        memory: "64 GB"
        preemptible: true
        maxRetries: 0
    }

	output {
		File out_bam_normal =  "~{output_normal_name}"
	}
}