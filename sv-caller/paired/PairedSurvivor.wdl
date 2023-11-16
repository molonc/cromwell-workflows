version 1.0

workflow PairedSurvivor {

    input {
		File manta_input
		File lumpy_input
		File gridss_input

        String tumor_name
        String normal_name
    }

    Array[File] inputs = [manta_input, lumpy_input, gridss_input]
    Array[String] callers = ['manta', 'lumpy', 'gridss']

    scatter (i in range(length(inputs))){
        call Filter {
            input:
                survivor_input = inputs[i],
                caller = callers[i]
        }
    }

    call Merge1 {
        input:
            survivor_input = Filter.outfile,
            tumor_name = tumor_name, 
            normal_name = normal_name
    }

    call Merge2 {
        input:
            survivor_input = Filter.outfile,
            tumor_name = tumor_name, 
            normal_name = normal_name
    }

    call Merge3 {
        input:
            survivor_input = Filter.outfile,
            tumor_name = tumor_name, 
            normal_name = normal_name
    }

    output {
        File merged1 = Merge1.outfile
        File merged2 = Merge2.outfile
        File merged3 = Merge3.outfile
    }


}

task Filter {
	input {
        File survivor_input
        String caller
	}
    
    command <<<
        SURVIVOR filter ~{survivor_input} NA 30 -1 0 -1 ~{caller + "_survivor.vcf"}
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/sv-caller-p:latest"
        cpu: 8
        memory: "24 GB"
        preemptible: true
        maxRetries: 0
    }

    output {
        File outfile = '~{caller + "_survivor.vcf"}'
    }
}

# Merges all lines together
task Merge1 {
	input {
        Array[File] survivor_input
        String tumor_name
        String normal_name
	}

    File input_1 = survivor_input[0]
    File input_2 = survivor_input[1]
    File input_3 = survivor_input[2]
    
    command <<<
        echo ~{input_1} > inputs.txt
        echo ~{input_2} >> inputs.txt
        echo ~{input_3} >> inputs.txt
   
        SURVIVOR merge inputs.txt 100 1 1 1 0 30 ~{tumor_name + "_" + normal_name + "_merged1.vcf"}
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/sv-caller-p:latest"
        cpu: 8
        memory: "24 GB"
        preemptible: true
        maxRetries: 0
    }

    output {
        File outfile = '~{tumor_name + "_" + normal_name + "_merged1.vcf"}'
    }
}

# Merges all lines that agree with at least one other caller
task Merge2 {
	input {
        Array[File] survivor_input
        String tumor_name
        String normal_name
	}

    File input_1 = survivor_input[0]
    File input_2 = survivor_input[1]
    File input_3 = survivor_input[2]
    
    command <<<
        echo ~{input_1} > inputs.txt
        echo ~{input_2} >> inputs.txt
        echo ~{input_3} >> inputs.txt
   
        SURVIVOR merge inputs.txt 100 2 1 1 0 30 ~{tumor_name + "_" + normal_name + "_merged2.vcf"}
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/sv-caller-p:latest"
        cpu: 8
        memory: "24 GB"
        preemptible: true
        maxRetries: 0
    }

    output {
        File outfile = '~{tumor_name + "_" + normal_name + "_merged2.vcf"}'
    }
}

# Merges all lines that agree with at least two other callers
task Merge3 {
	input {
        Array[File] survivor_input
        String tumor_name
        String normal_name
	}

    File input_1 = survivor_input[0]
    File input_2 = survivor_input[1]
    File input_3 = survivor_input[2]
    
    command <<<
        echo ~{input_1} > inputs.txt
        echo ~{input_2} >> inputs.txt
        echo ~{input_3} >> inputs.txt
   
        SURVIVOR merge inputs.txt 100 3 1 1 0 30 ~{tumor_name + "_" + normal_name + "_merged3.vcf"}
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/sv-caller-p:latest"
        cpu: 8
        memory: "24 GB"
        preemptible: true
        maxRetries: 0
    }

    output {
        File outfile = '~{tumor_name + "_" + normal_name + "_merged3.vcf"}'
    }
}



