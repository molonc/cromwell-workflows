version 1.0

import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/ubam-pre-pro/tasks/GermlineStructs.wdl"

# bams need to be name-sorted
workflow PairedLumpy {
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
        String param = '"%s"'
    }

    call Paired {
        input:
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

    output {
        File filtered = Paired.lumpy_out
        File unfiltered = Paired.lumpy_unfil
    }
}

task Paired {
	input {
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

    Int disk_size = ceil((size(tumor_bam, "GB") + size(normal_bam, "GB")) * 4)

    command <<<
        source /opt/conda/bin/activate py2

        mkdir temp

        lumpyexpress \
            -B ~{tumor_bam},~{normal_bam} \
            ~{"-x" + param} -o ~{tumor_name + "_" + normal_name + "_lumpy.unfiltered.vcf"} \
            -m 4 \
            -r 0 \
            -k \
            -v \
            -T ./temp
        # somatic + SV quality filtering
        #   'normal' sample assumes index 1
        bcftools filter -O v -o ~{tumor_name + "_" + normal_name + "_lumpy.vcf"} -i "FORMAT/SU[1] == 0 && FILTER == '.'" ~{tumor_name + "_" + normal_name + "_lumpy.unfiltered.vcf"}
    >>>

    runtime {
    	docker: "apariciobioinformaticscoop/sv-caller-c:latest"
    	disk: disk_size + " GB"
    	cpu: 16
    	memory: "64 GB"
    	preemptible: true
    	maxRetries: 1
    }

    output {
        File lumpy_out = '~{tumor_name + "_" + normal_name + "_lumpy.vcf"}'
        File lumpy_unfil = '~{tumor_name + "_" + normal_name + "_lumpy.unfiltered.vcf"}'
    }
}
