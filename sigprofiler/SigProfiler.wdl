version 1.0

# Filter FFPE noise
workflow SigProfiler {
    input {

        File run_script
        File matrix
    }

    scatter (idx in range(10)) {
        call SigProfiler {
            input:
                run_script = run_script,
                matrix = matrix,
                idx = idx
        }
    }
}

task SigProfiler {
    input {    
        File matrix

        File run_script

        Int idx
    }

    command <<<
        python ~{run_script} -o outputs --matrix ~{matrix} -n 20 --max_sigs 10

        cp ./outputs/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities.txt ./COSMIC_SBS96_Activities_~{idx}.txt
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/sigprofiler:feb29"
        disk: "24 GB"
        cpu: 8
        memory: "8 GB"
        preemptible: true
        maxRetries: 0
    }

    output {
        File outfile = 'COSMIC_SBS96_Activities_~{idx}.txt' # keep this
    }
}