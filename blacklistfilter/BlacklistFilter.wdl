version 1.0

# Mark duplicate reads to avoid counting non-independent observations
workflow BlacklistFilter {
    input {
        File input_bam
        String sample_name
        File blacklist
    }

    call BlacklistFilter {
        input:
            input_bam = input_bam,
            output_file_name = sample_name + ".postblacklist.bam",
            sample_name = sample_name,
            blacklist = blacklist
    }
}

task BlacklistFilter {
    input {
        File input_bam
    
        String output_file_name
        String sample_name

        File blacklist
    }

    Float multiplier = 15
    Int disk_size = ceil(size(input_bam, "GB") * multiplier) + 100 

    command <<<
        bedtools intersect -v -abam ~{input_bam} -b ~{blacklist} > ~{output_file_name}
        samtools index ~{output_file_name}
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/wasp-mapping:latest"
        disk: disk_size + " GB" 
        cpu: 40 # changed from 24
        memory: "64 GB" # changed from 24
        preemptible: true
        maxRetries: 3
    }


    output {
        File outfile = "~{output_file_name}"
        File index = "~{output_file_name + '.bai'}"
    }
}
