version 1.0

workflow CalcBafSignal {
    input {
        String study
        String sample_name
        String normal_name

        File plasma
        File tumor
        File normal
        File normal_og

        File merge_bafs_py
        File calc_baf_signal_py
        File joint_filters_py
    }

    call Merge {
        input:
            sample_name = sample_name,
            plasma = plasma, 
            tumor = tumor,
            normal = normal,
            normal_og = normal_og,
            merge_bafs_py = merge_bafs_py
    }

    call CalcBafSignal {
        input:
            sample_name = sample_name,
            combined = Merge.out,
            calc_baf_signal_py = calc_baf_signal_py,
            joint_filters_py = joint_filters_py
    }

    output {
        File merged = Merge.out
        File signal = CalcBafSignal.out
    }
}

task Merge {
    input {
        String sample_name

        File plasma
        File tumor
        File normal
        File normal_og

        File merge_bafs_py
    }

    command <<<
        python ~{merge_bafs_py} -p ~{plasma} -t ~{tumor} -n ~{normal} -o ~{sample_name} -g ~{normal_og}
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/mergebafs:latest"
        disk: "12 GB"
        cpu: 4
        memory: "16 GB"
        preemptible: true
        maxRetries: 0
    }

    output {
        File out = "~{sample_name}.tum.cnv.unfilt.txt"
    }
}

task CalcBafSignal {
    input {
        String sample_name

        File combined

        File calc_baf_signal_py
        File joint_filters_py
    }

    command <<<
        python ~{calc_baf_signal_py} -c ~{combined} -o ~{sample_name}
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/mergebafs:latest"
        dick: "12 GB"
        cpu: 4
        memory: "16 GB"
        preemptible: true
        maxRetries: 0
    }

    output {
        File out = "~{sample_name}.cnv.grp.snpconst.txt"
    }
}
