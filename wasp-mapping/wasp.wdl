version 1.0

# this is taken from the working version (august 29, 2023), but also has an added task to actually get the SA ids, which hasn't been tested yet. 
# for find_intersecting_snps.py, all the imported python scripts (util + snptable) has been compiled into one big file.  
# WORKFLOW DEFINITION
workflow WaspMapping {
    input {
        String study # metadata for clean-up automation 
        
        File input_vcf
        File input_bam
        String sample_name

        String read_group
        File reference_fasta
        File chrom_info

        File intersecting_snps_script
        File filter_remapped_script
    }

    call compressVcf {
        input: 
            vcf = input_vcf, 
            sample_name = sample_name, 
            chrom_info = chrom_info, 
            snp_index_name = sample_name + ".snp_index.h5",
            snp_tab_name = sample_name + ".snp_tab.h5", 
            haplotype_name = sample_name + ".h5"
    }

    # call splitVcfByChr {
    #     input: 
    #         compressed_vcf = compressVcf.compressed_vcf, 
    #         sample_name = sample_name
    # }

    # call VcfToH5 {
    #     input: 
    #         chrom_info = chrom_info, 
    #         vcfs = compressVcf.output_vcfs, 
    #         # vcfs = splitVcfByChr.output_vcfs, 
    #         snp_index_name = sample_name + ".snp_index.h5",
    #         snp_tab_name = sample_name + ".snp_tab.h5", 
    #         haplotype_name = sample_name + ".h5"
    # }

    call MakeSampleNameTxt {
        input: 
            chr_1_vcf = compressVcf.output_vcfs[0]
    }

    call FindIntersectingSnpsPairedEnd {
        input: 
            input_bam = input_bam, 
            snp_index = compressVcf.snp_index, 
            snp_tab = compressVcf.snp_tab, 
            haplotype = compressVcf.haplotype, 
            samples_txt = MakeSampleNameTxt.txt_file, 
            fastq1_name = sample_name + ".to.remap.fq1.gz", 
            fastq2_name = sample_name + ".to.remap.fq2.gz", 
            keep_bam = sample_name + ".keep.bam", 
            remap_bam = sample_name + ".to.remap.bam",
            intersecting_snps_script = intersecting_snps_script
    }

    call BwaAlignment {
        input: 
            input_fastq1 = FindIntersectingSnpsPairedEnd.fastq1,
            input_fastq2 = FindIntersectingSnpsPairedEnd.fastq2,
            read_group = read_group,
            reference_fasta = reference_fasta, 
            output_file_name = sample_name + ".realigned.bam"
    }

    call FilterRemappedReads {
        input: 
            realigned_bam = BwaAlignment.realigned_bam, 
            remapped_bam = FindIntersectingSnpsPairedEnd.remap_bam,
            output_file_name = sample_name + ".kept.bam", 
            filter_remapped_script = filter_remapped_script
    }

    call MergeBam {
        input: 
            keep_bam = FindIntersectingSnpsPairedEnd.keep_bam,
            kept_bam = FilterRemappedReads.kept_bam,
            output_file_name = sample_name + ".keep.merge.bam"
    }

    call SortBam {
        input: 
            merged_bam = MergeBam.merged_bam, 
            output_file_name = sample_name + "keep.merge.sort.bam"
    }

    call IndexBam {
        input: 
            sorted_bam = SortBam.sorted_bam, 
            output_file_name = sample_name + ".keep.merge.sort.bam.bai"
    }

    output {
        File keep_bam = FindIntersectingSnpsPairedEnd.keep_bam
        File kept_bam = FilterRemappedReads.kept_bam
        File fastq1 = FindIntersectingSnpsPairedEnd.fastq1
        File fastq2 = FindIntersectingSnpsPairedEnd.fastq2
        File sorted_bam = SortBam.sorted_bam
        File sorted_bam_index = IndexBam.bai
    }
}

# Preliminary task: compress VCF file with bgzip 
# edited to combine the compress + split tasks 
task compressVcf {
    input {
        File vcf
        String sample_name
        String file_name = basename(vcf) + ".gz"
        String tbi_name = basename(vcf) + ".gz.tbi"

        File chrom_info 
        String snp_index_name
        String snp_tab_name
        String haplotype_name 
    }

    Int disk_size = ceil(size(vcf, "GB") * 4)

    command <<<
        # compress and index vcf
        bgzip -c ~{vcf} > ~{file_name}
        tabix -p vcf ~{file_name}
        # split vcf by chromosome 
        mkdir vcf_by_chr
        for i in {1..22} 
        do
            bcftools view ~{file_name} --regions chr$i > vcf_by_chr/~{sample_name}.$i.het.vcf.gz 
        done
        # convert split vcfs to a HDF5 SNP file
        snp2h5 \
        --chrom ~{chrom_info} \
        --format vcf \
        --snp_index ~{snp_index_name} \
        --snp_tab ~{snp_tab_name} \
        --haplotype ~{haplotype_name} \
        vcf_by_chr/*.het.vcf.gz
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/wasp-mapping:latest"
        disk: disk_size + " GB"
        cpu: 12
        memory: "12 GB"
        preemptible: true
        maxRetries: 3
    }

    output {
        File compressed_vcf = "~{file_name}"
        File vcf_index = "~{tbi_name}"
        Array[File] output_vcfs = glob("vcf_by_chr/*.het.vcf.gz") 
        File snp_index = "~{snp_index_name}" # keep this output
        File snp_tab = "~{snp_tab_name}" # keep this output
        File haplotype = "~{haplotype_name}" # keep this output
    }
}


# split pre-processing outputted VCF into multiple VCFs, sorted by chromosome 
# A. CHECK: syntax of the regex after --regions against the mrdedge example. 
# task splitVcfByChr {
#     input {
#         File compressed_vcf
#         String sample_name
#     }

#     Int disk_size = ceil(size(compressed_vcf, "GB") * 4)

#     command <<<
#         mkdir vcf_by_chr
#         for i in {1..22} 
#         do
#             bcftools view ~{compressed_vcf} --regions chr$i > vcf_by_chr/~{sample_name}.$i.het.vcf.gz 
#         done
#     >>>

#     runtime {
#         docker: "apariciobioinformaticscoop/wasp-mapping:latest"
#         disk: disk_size + " GB"
#         cpu: 12
#         memory: "12 GB"
#         preemptible: true
#         maxRetries: 3
#     }

#     output {
#         Array[File] output_vcfs = glob(vcf_by_chr/*.het.vcf.gz) # unsure of this
#     }
# }

# B. Convert input VCF to a HDF5 SNP file
# task VcfToH5 {
#     input {
#         File chrom_info # stored as a reference in inputs. 
#         Array[File] vcfs
#         String snp_index_name
#         String snp_tab_name 
#         String haplotype_name
#     }

#     Int disk_size = ceil(size(vcfs, "GB") * 2)

#     command <<<
#         # mv ~{vcfs} .
#         snp2h5 \
#         --chrom ~{chrom_info} \
#         --format vcf \
#         --snp_index ~{snp_index_name} \
#         --snp_tab ~{snp_tab_name} \
#         --haplotype ~{haplotype_name} \
        
#     >>>

#     runtime {
#         docker: "apariciobioinformaticscoop/wasp-mapping:latest"
#         disk: disk_size + " GB"
#         cpu: 16
#         memory: "16 GB"
#         preemptible: true
#         maxRetries: 3
#     }

#     output {
#         File snp_index = "~{snp_index_name}" # keep this output
#         File snp_tab = "~{snp_tab_name}" # keep this output
#         File haplotype = "~{haplotype_name}" # keep this output
#     }
# }

task MakeSampleNameTxt {
    input {
        File chr_1_vcf
        String output_txt_file = "samples.txt"
    }

    # command <<<
    #     echo ~{sample_name} > ~{output_txt_file}
    #     grep "\S" ~{output_txt_file}
    # >>>

    # heredoc command: 
    # command <<<
    #     echo "SA1145N" > ~{output_txt_file}
    #     echo "SA1145X1" >> ~{output_txt_file}
    # >>>

    command <<<
        awk '{FS="="}; /##normal_sample/ {print $2}' ~{chr_1_vcf} >> ~{output_txt_file}
        awk '{FS="="}; /##tumor_sample/ {print $2}' ~{chr_1_vcf} >> ~{output_txt_file}
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/wasp-mapping:latest"
        disk: "6 GB"
        cpu: 6
        memory: "6 GB"
        preemptible: true
        maxRetries: 3
    }

    output {
        File txt_file = "~{output_txt_file}"
    }
}

# C. find intersecting SNPs using find_intersecting_snps.py 
task FindIntersectingSnpsPairedEnd {
    input {
        File input_bam
        File snp_index
        File snp_tab
        File haplotype
        File samples_txt
        String fastq1_name
        String fastq2_name
        String keep_bam
        String remap_bam
        String output_dir = "find_intersecting_snps"
        File intersecting_snps_script
    }

    Float multiplier = 2.5
    Int disk_size = ceil(size(input_bam, "GB") * multiplier) + 100 # also bumped up runtime metrics

    # BAM_PATH=$(pwd ~{input_bam})/$(ls ~{input_bam})
    # echo $BAM_PATH
    
    command <<<
        python ~{intersecting_snps_script} \
        --is_paired_end \
        --is_sorted \
        --output_dir ~{output_dir} \
        --snp_tab ~{snp_tab} \
        --snp_index ~{snp_index} \
        --samples ~{samples_txt} \
        --haplotype ~{haplotype} ~{input_bam}
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/wasp-mapping:latest"
        disk: disk_size + " GB" 
        cpu: 24
        memory: "24 GB"
        preemptible: true
        maxRetries: 3
    }


    output {
        # File sorted_bam = "~{sorted_bam_name}"
        # File sorted_bai = "~{sorted_bai_name}"
        File fastq1 = "~{fastq1_name}"
        File fastq2 = "~{fastq2_name}"
        File keep_bam = "~{keep_bam}" # reads that didn't intersect with SNPs
        File remap_bam = "~{remap_bam}" # reads that intersected with SNPs, needs to be flipped + remapped
    }
}

# D. Realign remap.bam reads using BWA aligner again
task BwaAlignment {
    input {
        File input_fastq1
        File input_fastq2
        String read_group
        File reference_fasta
        String output_file_name
    }

    Float multiplier = 2.5
    Int disk_size = ceil(size(input_fastq1, "GB") * multiplier) + 200

    command <<<
        bwa mem -Y -K 100000000 -t 8 -R ~{read_group} ~{reference_fasta} ~{input_fastq1} ~{input_fastq2} | samtools view -Shb -o ~{output_file_name}
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/wasp-mapping:latest"
        disk: disk_size + " GB"
        cpu: 18
        memory: "24 GB"
        preemptible: true
        maxRetries: 3
    }

    output {
        File realigned_bam = "~{output_file_name}"
    }
}

# E. compare realigned reads (realigned bam) vs. original reads (remap.bam) using filter_remapped_reads
task FilterRemappedReads {
    input {
        File realigned_bam 
        File remapped_bam 
        String output_file_name
        File filter_remapped_script
    }

    Float multiplier = 2.5
    Int disk_size = ceil(size(realigned_bam, "GB") * multiplier) + 50
    

    command <<<
        python ~{filter_remapped_script} \
        ~{remapped_bam} \
        ~{realigned_bam} ~{output_file_name}
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/wasp-mapping:latest"
        disk: disk_size + " GB"
        cpu: 18
        memory: "24 GB"
        preemptible: true
        maxRetries: 3
    }


    output {
        File kept_bam = "~{output_file_name}"
    }
}

# F. merge kept.bam + keep.bam
task MergeBam {
    input {
        File keep_bam 
        File kept_bam
        String output_file_name
    }

    Float multiplier = 2.5
    Int disk_size = ceil(size(keep_bam, "GB") * multiplier) + 50

    command <<<
        samtools merge ~{output_file_name} ~{keep_bam} ~{kept_bam} -@ 4
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/wasp-mapping:latest"
        disk: disk_size + " GB"
        cpu: 15
        memory: "20 GB"
        preemptible: true
        maxRetries: 3
    }


    output {
        File merged_bam = "~{output_file_name}"
    }
}

# G. sort merged bam 
task SortBam {
    input {
        File merged_bam
        String output_file_name
    }

    # based these disk numbers on GATK pre-processing tasks
    Float sort_sam_disk_multiplier = 3.25
    Int disk_size = ceil(sort_sam_disk_multiplier * size(merged_bam, "GB")) + 20

    command <<<
        samtools sort -o ~{output_file_name} ~{merged_bam} -@ 4
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/wasp-mapping:latest"
        disk: disk_size + " GB"
        cpu: 8
        memory: "10 GB"
        preemptible: true
        maxRetries: 3
    }

    output {
        File sorted_bam = "~{output_file_name}" # keep this
    }
}

# H. create bam index 
task IndexBam {
    input {
        File sorted_bam
        String output_file_name
    }

    Float multiplier = 2.5
    Int disk_size = ceil(size(sorted_bam, "GB") * multiplier) + 20 

    command <<<
        samtools index ~{sorted_bam}
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/wasp-mapping:latest"
        disk: disk_size + " GB"
        cpu: 12
        memory: "15 GB"
        preemptible: true
        maxRetries: 3
    }

    output {
        File bai = "~{output_file_name}" # keep this
    }
}


