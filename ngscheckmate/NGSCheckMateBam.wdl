version 1.0

# Mark duplicate reads to avoid counting non-independent observations
workflow NGSCheckMate {
    input {
        Array[File] bams

        Array[String] sample_names
        String patient_name

        File ID

        File ngscheckmate_script
        File reference_bed
        File reference
    }


    scatter (i in range(length(bams))) {
        call GetVCF {
            input:
                bam = bams[i],
                sample_name = sample_names[i],
                reference_bed = reference_bed,
                reference = reference
        }
    }

    call NGSCheckMate {
        input:
            vcfs = GetVCF.vcf,
            patient_name = patient_name,
            ID = ID,
            ngscheckmate_script = ngscheckmate_script,
            reference_bed = reference_bed
    }

    output {
        File outfile = NGSCheckMate.outfile
    }
}

task GetVCF {
    input {
        File bam

        String sample_name

        File reference_bed
        File reference
    }

    Float multiplier = 6
    Int disk_size = ceil(size(bam, "GB") * multiplier) + 50

    command <<<
        samtools mpileup -I -uf ~{reference} -l ~{reference_bed} ~{bam} | bcftools call -c - > ./~{sample_name}.vcf
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/wasp-mapping:latest"
        disk: disk_size + " GB"
        cpu: 12
        memory: "15 GB"
        preemptible: true
        maxRetries: 0
    }

    output {
        File vcf = "~{sample_name}.vcf"
    }
}

task NGSCheckMate {
    input {
        Array[File] vcfs
        
        String patient_name

        File ID

        File ngscheckmate_script
        File reference_bed

        Int len = length(vcfs)
    }

    Float multiplier = 6
    Int disk_size = ceil(size(vcfs[0], "GB") * multiplier) + 4

    command <<<
        mkdir outputs
        mkdir data

        ARRAY=(~{sep=" " vcfs})
        for (( c = 0; c < ~{len}; c++ ))
        do
            mv ${ARRAY[${c}]} ./data
        done

        python ~{ngscheckmate_script} -V -l ~{ID} -bed ~{reference_bed} -O outputs

        tar -czvf ~{patient_name}.tar.gz ./outputs/
    >>>

    runtime {
        docker: "quay.io/biocontainers/ngscheckmate:1.0.1--py27pl5321r40hdfd78af_1"
        disk: disk_size + " GB"
        cpu: 12
        memory: "15 GB"
        preemptible: true
        maxRetries: 0
    }

    output {
        File outfile = "~{patient_name}.tar.gz" # keep this
    }
}