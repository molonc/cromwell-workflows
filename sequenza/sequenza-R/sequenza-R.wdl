version 1.0
workflow SequenzaR {
    input {
        File seqz
        String sample_ID
        File script
    }

    call Untar {
        input:
        seqz = seqz,
        sample_ID = sample_ID
    }

    call SequenzaR{
        input:
        seqz = Untar.outfile,
        sample_ID = sample_ID,
        script = script
    }
}

task Untar {
    input {
        File seqz
        String sample_ID
    }
  
    command <<<
        tar -xzvf ~{seqz}
    >>>

    runtime {
        docker: "sequenza/sequenza"
        memory: "24 GB" # 24 -> 8 -> 16
        cpu: 8 # 10 -> 5 -> 8
        disk: "64 GB" 
    }

    output {
        File outfile = "~{sample_ID}_bin50.seqz.gz"
    }
}

task SequenzaR {
    input {
        File seqz
        String sample_ID
        File script
    }
  
    command <<<

        mkdir TEST

        Rscript --vanilla ~{script} ~{seqz} ~{sample_ID}
        tar -czvf ~{sample_ID}.tar.gz ./TEST/
    >>>

    runtime {
        docker: "sequenza/sequenza"
        memory: "24 GB" # 24 -> 8 -> 16
        cpu: 8 # 10 -> 5 -> 8
        disk: "64 GB" 
    }

    output {
        File outfile = "~{sample_ID}.tar.gz"
    }
}