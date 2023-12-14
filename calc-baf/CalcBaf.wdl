version 1.0

import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/ubam-pre-pro/tasks/GermlineStructs.wdl"

# WORKFLOW DEFINITION
workflow CalcBaf {
    input {
        String study # metadata for clean-up automation 

        File input_bam
        File input_bai

        File input_vcf
        File input_vcf_index

        File cnv

        String sample_name
        String normal_name
        String haplotype_caller_version

        ReferenceFasta references
        File seg_script
        File calc_baf
        File remove_indels_script
    }

    File reference_fasta = references.ref_fasta

    call GenSeg {
        input:
            cnv = cnv,
            seg_script = seg_script,
            sample_name = sample_name,
            normal_name = normal_name
    }

    call AlleleCount {
        input: 
            bam = input_bam, 
            bai = input_bai,
            vcf = input_vcf,
            vcf_index = input_vcf_index,
            s_name = sample_name,
            n_name = normal_name,
            h_call_ver = haplotype_caller_version,
            ref = references,
            ref_fasta = reference_fasta
    }

    # call MergeCounts {
    #     input:
    #         input_txts = AlleleCount.output_txts,
    #         s_name = sample_name,
    #         n_name = normal_name,
    #         h_call_ver = haplotype_caller_version
    # }

    call CalcBaf {
        input: 
            het_tar = AlleleCount.het_tar,
            vcf = input_vcf,
            s_name = sample_name,
            n_name = normal_name,
            h_call_ver = haplotype_caller_version,
            calc_baf = calc_baf
    }

    call RemoveIndels{
        input:
            script = remove_indels_script, 
            input_baf = CalcBaf.baf_txt
    }

    call MergeAlleleCounts {
        input:
            input_baf = RemoveIndels.outfile,
            haplotypecalls = input_vcf,
            vcf_index = input_vcf_index,
            cnv_calls = GenSeg.seg_file,
            s_name = sample_name,
            n_name = normal_name,
            h_call_ver = haplotype_caller_version
    }
}

task GenSeg {
    input {
        File cnv
        File seg_script

        String sample_name
        String normal_name
    }

    Int disk_size = ceil(size(cnv, "GB") * 4)

    command <<<
        ~{seg_script} -i ~{cnv} -o ~{sample_name + "_" + normal_name +".seg.txt"} -t sequenza
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/genseg:latest"
        disk: disk_size + " GB"
        cpu: 12
        memory: "32 GB"
        preemptible: true
        maxRetries: 3
    }

    output {
        File seg_file = '~{sample_name + "_" + normal_name +".seg.txt"}'
    }
}

# A. Calculate Allele Counts by Chromosome
task AlleleCount {
    input {
        File bam
        File bai

        File vcf
        File vcf_index

        String s_name
        String n_name
        String h_call_ver

        ReferenceFasta ref
        File ref_fasta
    }

    Int disk_size = ceil(size(bam, "GB") * 4)

    command <<<
        mkdir het

        for chr in {1..22}; do
            samtools mpileup -f ~{ref_fasta} -q 10 -Q 10 --ff UNMAP --ff SECONDARY --ff QCFAIL --ff DUP -r chr${chr} -l ~{vcf} -B ~{bam} | sequenza-utils pileup2acgt -p - -o ./het/~{s_name}--~{n_name}.haplotypeCalls.~{h_call_ver}.het.alleles."$chr".txt
        done

        head -1 ./het/~{s_name}--~{n_name}.haplotypeCalls.~{h_call_ver}.het.alleles.1.txt > ./het/~{s_name}--~{n_name}.haplotypeCalls.~{h_call_ver}.het.alleles.txt
        for i in `seq 1 22`; do
            grep -w -v "^chr" ./het/~{s_name}--~{n_name}.haplotypeCalls.~{h_call_ver}.het.alleles."$i".txt >> ./het/~{s_name}--~{n_name}.haplotypeCalls.~{h_call_ver}.het.alleles.txt
        done

        tar -czvf het.tar.gz ./het/
    >>>

    runtime {
        docker: "sequenza/sequenza" # Sequenza docker
        disk: disk_size + " GB"
        cpu: 12
        memory: "32 GB"
        preemptible: true
        maxRetries: 3
    }

    output {
        File het_tar = "het.tar.gz"
    }
}

# B. Merge allele counts across chromosomes
# task MergeCounts {
#     input {
#         Array[File] input_txts

#         String s_name
#         String n_name
#         String h_call_ver
#     }

#     command <<<
#         head -1 ~{s_name}--~{n_name}_haplotypeCalls_~{h_call_ver}_het_alleles_1.txt > ~{s_name}--~{n_name}_haplotypeCalls_~{h_call_ver}_het_alleles.txt
#         for i in `seq 1 22`; do
#             grep -w -v "^chr" ~{s_name}--~{n_name}_haplotypeCalls_~{h_call_ver}_het_alleles_"$i".txt >> ~{s_name}--~{n_name}_haplotypeCalls_~{h_call_ver}_het_alleles.txt
#         done
#     >>>

#     runtime {
#         docker: "apariciobioinformaticscoop/wasp-mapping:latest"
#         disk: "2 GB"
#         cpu: 1
#         memory: "1 GB"
#         preemptible: true
#         maxRetries: 3
#     }

#     output {
#         File merged_output = "~{s_name}--~{n_name}_haplotypeCalls_~{h_call_ver}_het_alleles.txt"
#     }
# }

# C. Run calc_baf.py on all samples to format allele counts
task CalcBaf {
    input {
        File het_tar
        File vcf

        String s_name
        String n_name
        String h_call_ver
        
        File calc_baf
    }

    command <<<
        tar -xzvf ~{het_tar}
        python ~{calc_baf} -i ./het -o $PWD -he ~{vcf} -t ~{s_name} -n ~{n_name} -v ~{h_call_ver}
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/wasp-mapping:latest"
        disk: "12 GB"
        cpu: 4
        memory: "4 GB"
        preemptible: true
        maxRetries: 3
    }

    output {
        File baf_txt = "~{s_name}--~{n_name}.haplotypeCalls.~{h_call_ver}.het.baf.txt"
    }
}

task RemoveIndels {
    input {
        File script 
        File input_baf
    }

    command <<<
        Rscript --vanilla ~{script} ~{input_baf}
    >>>
    runtime {
        docker: "apariciobioinformaticscoop/ichorcna-updated-packages:latest"
        disk: "4 GB"
        cpu: 1
        memory: "4 GB"
        preemptible: true
        maxRetries: 0
    }

    output {
        File outfile = "tmp.baf.txt"
        File metrics = "metrics.txt"
    }


}

task MergeAlleleCounts {
    input {
        File input_baf
        File haplotypecalls
        File vcf_index
        File cnv_calls

        String s_name
        String n_name
        String h_call_ver
    }

    command <<<
        grep -v "^CHROM" ~{input_baf} | awk -F"\t" '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > tmp.bed
        bcftools view -O v ~{haplotypecalls} | awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' > tmp_hap.bed
        sed '1d' ~{cnv_calls} > cnv.bed

        echo -ne "chr\tstart\tend\tref\talt\tref_cnt\talt_cnt\tbaf\thc_qual\thc_filter\thc_info\thc_format\thc_quals\twin_chr\twin_start\twin_end\tTotal_cn\tMinor_cn\tCall\tSize\toverlap\n" > ~{s_name}.haplotypeCalls.~{h_call_ver}.het.baf.blacklist.bed
        bedtools intersect -wb -a tmp.bed -b tmp_hap.bed | bedtools intersect -wo -a stdin -b cnv.bed | cut -f '1-8,15-25,27-28' >> ~{s_name}.haplotypeCalls.~{h_call_ver}.het.baf.blacklist.bed
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/wasp-mapping:latest"
        disk: "16 GB"
        cpu: 4
        memory: "8 GB"
        preemptible: true
        maxRetries: 0
    }

    output {
        File bed = "~{s_name}.haplotypeCalls.~{h_call_ver}.het.baf.blacklist.bed"
    }
}


