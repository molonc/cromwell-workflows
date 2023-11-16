version 1.0

import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/ubam-pre-pro/tasks/GermlineStructs.wdl"

workflow SingleGridss {
	input {
		Boolean exclude_regions

		String tumor_name
		File tumor_bam
		File tumor_bai

        ReferenceFasta references
    }

    File reference_fasta = references.ref_fasta

	meta {allowNestedInputs: true}

    if (exclude_regions) {
        String param = "BLACKLIST={}"
    }

    call Single {
            input:
                exclude_regions = exclude_regions,
                tumor_name = tumor_name,
                tumor_bam = tumor_bam,
                tumor_bai = tumor_bai,
                reference_fasta = reference_fasta,
                reference_fai = references,
                param = param
        }

        call SingleFilter {
            input:
                input_file = Single.gridss_out,
                tumor_name = tumor_name
        }

    output {
        File filtered = SingleFilter.gridss_out
        File unfiltered = Single.gridss_unfil
    }

}

task Single {
	input {
        Boolean exclude_regions

		String tumor_name
		File tumor_bam
        File tumor_bai

		File reference_fasta
        ReferenceFasta reference_fai
        String? param
	}

    Int disk_size = ceil(size(tumor_bam, "GB") * 4)  

    command <<<
        gridss -r ~{reference_fasta} -o ~{tumor_name + "_gridss.unfiltered.vcf"} -a ./gridss_assembly.bam ~{tumor_bam}
    >>>

    runtime {
        docker: "gridss/gridss:latest"
        disk: disk_size + " GB"
        cpu: 8
        memory: "31 GB"
        preemptible: true
        maxRetries: 0
    }

    output {
        File gridss_out = '~{tumor_name + "_gridss.unfiltered.vcf"}'
        File gridss_unfil = '~{tumor_name + "_gridss.unfiltered.vcf"}'
    }
}

task SingleFilter {
	input {
        File input_file

		String tumor_name
	}

    command <<<
        bcftools filter -O v -o ~{tumor_name + "_gridss.vcf"} -i "FILTER == '.'" ~{input_file}
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/wasp-mapping:latest"
        cpu: 8
        memory: "16 GB"
        preemptible: true
        maxRetries: 0
    }

    output {
        File gridss_out = '~{tumor_name + "_gridss.vcf"}'
    }
}