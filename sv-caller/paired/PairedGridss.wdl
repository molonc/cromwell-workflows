version 1.0

import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/ubam-pre-pro/tasks/GermlineStructs.wdl"

workflow PairedGridss {
	input {
		Boolean exclude_regions

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

    if (exclude_regions) {
        String param = "BLACKLIST={}"
    }

    call Paired {
        input: 
            exclude_regions = exclude_regions,
            tumor_name = tumor_name,
            tumor_bam = tumor_bam,
            tumor_bai = tumor_bai,
            normal_name = normal_name,
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            reference_fasta = reference_fasta,
            reference_fai = references,
            param = param
    }

    call PairedFilter {
        input:
            input_file = Paired.gridss_out,
            tumor_name = tumor_name,
            normal_name = normal_name  
    }

    output {
        File filtered = PairedFilter.gridss_out
        File unfiltered = Paired.gridss_out
    }

}

task Paired {
	input {
        Boolean exclude_regions

		String tumor_name
		File tumor_bam
        File tumor_bai

        String normal_name
		File normal_bam
		File normal_bai

		File reference_fasta
        ReferenceFasta reference_fai
        String? param
	}

    Int disk_size = ceil(size(tumor_bam, "GB") * 4)  

    command <<<
        gridss -r ~{reference_fasta} -o ~{tumor_name + "_"+ normal_name + "_gridss.unfiltered.vcf"} -a ./gridss_assembly.bam ~{normal_bam} ~{tumor_bam}
    >>>

    runtime {
        docker: "gridss/gridss:latest"
        disk: disk_size + " GB"
        cpu: 24
        memory: "64 GB"
        preemptible: true
        maxRetries: 0
    }

    output {
        File gridss_out = '~{tumor_name + "_"+ normal_name + "_gridss.unfiltered.vcf"}'
    }
}

task PairedFilter {
	input {
        File input_file

		String tumor_name
        String? normal_name
	}

    command <<<
        bcftools filter -O v -o ~{tumor_name + "_" + normal_name + "_gridss.vcf"} -i "FILTER == 'PASS'" ~{input_file}
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/wasp-mapping:latest"
        cpu: 8
        memory: "16 GB"
        preemptible: true
        maxRetries: 0
    }

    output {
        File gridss_out = '~{tumor_name + "_" + normal_name + "_gridss.vcf"}'
    }
}