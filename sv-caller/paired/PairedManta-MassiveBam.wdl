version 1.0

import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/ubam-pre-pro/tasks/GermlineStructs.wdl"

workflow PairedManta {
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


    call Paired {
        input:
            tumor_name = tumor_name,
            tumor_bam = tumor_bam,
            tumor_bai = tumor_bai,
            normal_name = normal_name,
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            reference_fasta = reference_fasta,
            reference_fai = references
    }

    output {
        File filtered = Paired.manta_out
        File unfiltered = Paired.manta_unfil
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
	}

    Int disk_size = ceil((size(tumor_bam, "GB") + size(normal_bam, "GB")) * 6)

	command <<<
        source /opt/conda/bin/activate py2

        configManta.py \
            --runDir $PWD \
            --reference ~{reference_fasta} \
            --tumorBam ~{tumor_bam} \
            --normalBam ~{normal_bam}
        ./runWorkflow.py \
            -m local \
            -j 24 \
            -g 60
        cp ./results/variants/somaticSV.vcf.gz ./
        # SV quality filtering
        bcftools filter -O v -o ~{tumor_name + "_" + normal_name + "_manta.vcf"} -i "FILTER == 'PASS'" ./somaticSV.vcf.gz

        gunzip ./somaticSV.vcf.gz
        mv ./somaticSV.vcf ~{tumor_name + "_" + normal_name + "_manta.unfiltered.vcf"}
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/sv-caller-c:latest"
        disk: disk_size + " GB"
        cpu: 32
        memory: "80 GB"
        preemptible: true
        maxRetries: 1
    }

    output {
        File manta_out = '~{tumor_name + "_" + normal_name + "_manta.vcf"}'
        File manta_unfil = '~{tumor_name + "_" + normal_name + "_manta.unfiltered.vcf"}'
    }
}
