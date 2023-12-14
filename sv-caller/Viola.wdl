version 1.0

workflow Viola {

    input {
		File manta_input
		File lumpy_input
		File gridss_input

        File viola_script_m
        File viola_script_l
        File viola_script_g
    }

    call Viola_M {
        input:
            v_input = manta_input,
            viola_script = viola_script_m
    }

    call Viola_G {
        input:
            v_input = gridss_input,
            viola_script = viola_script_g
    }

    call Viola_L {
        input:
            v_input = lumpy_input,
            viola_script = viola_script_l
    }

    output {
        File manta_out = Viola_M.outfile
        File lumpy_out = Viola_L.outfile
        File gridss_out = Viola_G.outfile
    }
}

task Viola_M {
	input {
        File v_input
        File viola_script
	}

    Int disk_size = ceil(size(v_input, "GB") * 6)  
    
    command <<<
        python ~{viola_script} ~{v_input}
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/sv-caller-p:latest"
        disk: disk_size + " GB"
        cpu: 24
        memory: "64 GB"
        preemptible: true
        maxRetries: 0
    }

    output {
        File outfile = 'manta_out.vcf'
    }
}

task Viola_G {
	input {
        File v_input
        File viola_script
	}

    Int disk_size = ceil(size(v_input, "GB") * 6)  
    
    command <<<
        python ~{viola_script} ~{v_input}
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/sv-caller-p:latest"
        disk: disk_size + " GB"
        cpu: 24
        memory: "64 GB"
        preemptible: true
        maxRetries: 0
    }

    output {
        File outfile = 'gridss_out.vcf'
    }
}

task Viola_L {
	input {
        File v_input
        File viola_script
	}

    Int disk_size = ceil(size(v_input, "GB") * 6)  
    
    command <<<
        python ~{viola_script} ~{v_input}
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/sv-caller-p:latest"
        disk: disk_size + " GB"
        cpu: 24
        memory: "64 GB"
        preemptible: true
        maxRetries: 0
    }

    output {
        File outfile = 'lumpy_out.vcf'
    }
}

