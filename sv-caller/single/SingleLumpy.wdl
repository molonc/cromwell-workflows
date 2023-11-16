version 1.0

import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/ubam-pre-pro/tasks/GermlineStructs.wdl"

# bams need to be name-sorted
workflow SingleLumpy {
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
        String param = '"%s"'
    }

    call Single {
        input:
            tumor_name = tumor_name,
            tumor_bam = tumor_bam,
            tumor_bai = tumor_bai,
            reference_fasta = reference_fasta,
            reference_fai = references,
            param = param
	}

    output {
        File filtered = Single.lumpy_out
        File unfiltered = Single.lumpy_unfil
    }
}

task Single {
	input {
        String tumor_name
		File tumor_bam
		File tumor_bai

		File reference_fasta
        ReferenceFasta reference_fai

        String? param
    }

    Int disk_size = ceil(size(tumor_bam, "GB") * 4)

    command <<<
        source /opt/conda/bin/activate py2

        mkdir temp

        lumpyexpress \
            -B ~{tumor_bam} \
            ~{"-x" + param} -o ~{tumor_name + "_lumpy.unfiltered.vcf"} \
            -m 4 \
            -r 0 \
            -k \
            -v \
            -T ./temp
        # SV quality filtering
        bcftools filter \
            -O v \
            -o ~{tumor_name + "_lumpy.vcf"} \
            -i "FILTER == '.'" \
            ~{tumor_name + "_lumpy.unfiltered.vcf"}
    >>>
    
    runtime {
        docker: "apariciobioinformaticscoop/sv-caller-c:latest"
        disk: disk_size + " GB"
        cpu: 24
        memory: "64 GB"
        preemptible: true
        maxRetries: 0
    }

    output {
        File lumpy_out =  '~{tumor_name + "_lumpy.vcf"}'
        File lumpy_unfil = '~{tumor_name + "_lumpy.unfiltered.vcf"}'
    }
}