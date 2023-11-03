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
        String docker
	}

    File reference_fasta = references.ref_fasta

	meta {allowNestedInputs: true}

	call PairedSort {
		input:
			tumor_bam = tumor_bam,
			tumor_bai = tumor_bai,
			output_tumor_name = tumor_name + ".sort.bam",
			normal_bam = normal_bam,
			normal_bai = normal_bai,
			output_normal_name = normal_name + ".sort.bam"
	}

	output {
		File out_tumor_bam = PairedSort.out_bam_tumor 
		File out_normal_bam = PairedSort.out_bam_normal
	}
}

task PairedSort {
    input {
		File tumor_bam
		File tumor_bai
        String output_tumor_name

		File normal_bam
		File normal_bai
		String output_normal_name
    }

    Int disk_size = ceil(size(tumor_bam, "GB") * 4)

    command <<<
        samtools sort -n -o ~{output_tumor_name} ~{tumor_bam}
		samtools sort -n -o ~{output_normal_name} ~{normal_bam}
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
		File out_bam_normal =  "~{output_normal_name}"
	}
}
