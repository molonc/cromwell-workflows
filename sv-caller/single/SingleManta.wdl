version 1.0

import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/ubam-pre-pro/tasks/GermlineStructs.wdl"

workflow SingleManta {
	input {
		String tumor_name
		File tumor_bam
		File tumor_bai

        ReferenceFasta references
        String docker
	}

    File reference_fasta = references.ref_fasta

	meta {allowNestedInputs: true}


    call Single {
        input:
            tumor_name = tumor_name,
            tumor_bam = tumor_bam,
            tumor_bai = tumor_bai,
            reference_fasta = reference_fasta,
            reference_fai = references,
            docker = docker
		}

    output {
        File outfile = Single.manta_out
    }

}

task Single {
	input {
		String tumor_name
		File tumor_bam
		File tumor_bai

		File reference_fasta
        ReferenceFasta reference_fai
        String docker
	}

	command <<<
        source /opt/conda/bin/activate py2

        configManta.py \
            --runDir $PWD \
            --reference ~{reference_fasta} \
            --tumorBam ~{tumor_bam}
        ./runWorkflow.py \
            -m local \
            -j 20 \
            -g 60 \
        # SV quality filtering
        bcftools filter \
            -O v \
            -o ~{tumor_name + "_manta.vcf"} \
            -i "FILTER == 'PASS'" ./results/variants/tumorSV.vcf.gz
    >>>

    runtime {
        docker: docker
        cpu: 24
        memory: "64 GB"
        preemptible: true
        maxRetries: 0
    }

    output {
        File manta_out = '~{tumor_name + "_manta.vcf"}'
    }
}