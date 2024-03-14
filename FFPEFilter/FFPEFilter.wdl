version 1.0

# Filter FFPE noise
workflow FFPEFilter {
    input {
        File vcf
        String sample_name

        String dp
        String ad
        String af

        File filter_script
    }

    call FFPEFilter {
        input:
            vcf = vcf,
            sample_name = sample_name,
            dp = dp,
            ad = ad,
            af = af,
            filter_script = filter_script
    }

    call FFPEFilterIndex {
        input:
            vcf = FFPEFilter.outfile,
            sample_name = sample_name
    }

    output {
        File outfile = FFPEFilter.outfile
        File index = FFPEFilterIndex.outfile
    }
}

task FFPEFilter {
    input {
        File vcf
        
        String sample_name

        String dp
        String ad
        String af

        File filter_script
    }

    command <<<
        python ~{filter_script} -vcf ~{vcf} -o ~{sample_name + "-refiltered.vcf"} -dp ~{dp} -af ~{af} -ad ~{ad}
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/ffpe-filter:latest"
        disk: "4 GB"
        cpu: 4
        memory: "4 GB"
        preemptible: true
        maxRetries: 0
    }

    output {
        File outfile = '~{sample_name + "-refiltered.vcf"}' # keep this
    }
}

task FFPEFilterIndex {
    input {
      File vcf
      String sample_name
    }

    command {
        gatk IndexFeatureFile --feature-file ~{vcf} --output ~{sample_name + "-refiltered.vcf.idx"} 
    }

    runtime {
        docker: "broadinstitute/gatk:4.1.3.0"
        #bootDiskSizeGb: runtime_params.boot_disk_size
        memory: 24 + " GB"
        disk: 200 + " GB"
        preemptible: true
        maxRetries: 0
        cpu: 14
    }

    output {
        File outfile = '~{sample_name + "-refiltered.vcf.idx"}'
    }
}