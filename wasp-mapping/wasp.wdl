version 1.0

# for find_intersecting_snps.py, all the imported python scripts (util + snptable) has been 
# WORKFLOW DEFINITION
import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/ubam-pre-pro/tasks/GermlineStructs.wdl"

workflow WaspMapping {
    input {
        String study # metadata for clean-up automation 
        
        File input_vcf
        File input_vcf_index
        File input_bam
        File input_bai
        String sample_name
        String normal_name # metadata for clean-up automation 
        String output_dir
        String sa_id_tumor
        String sa_id_normal

        String read_group
        ReferenceFasta references
        File chrom_info

        File intersecting_snps_script
        File filter_remapped_script
    }

    File reference_fasta = references.ref_fasta

    call compressVcf {
        input: 
            vcf = input_vcf,
            vcf_index = input_vcf_index, 
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
            sample_name = sample_name,
            sa_id_tumor = sa_id_tumor,
            sa_id_normal = sa_id_normal          
    }

    call FindIntersectingSnpsPairedEnd {
        input: 
            input_bam = input_bam,
            input_bai = input_bai, 
            snp_index = compressVcf.snp_index, 
            snp_tab = compressVcf.snp_tab, 
            haplotype = compressVcf.haplotype, 
            samples_txt = MakeSampleNameTxt.txt_file, 
            fastq1_name = sample_name + ".remap.fq1.gz", 
            fastq2_name = sample_name + ".remap.fq2.gz", 
            keep_bam = sample_name + ".keep.bam", 
            remap_bam = sample_name + ".to.remap.bam",
            output_dir = output_dir,
            intersecting_snps_script = intersecting_snps_script
    }

    call BwaAlignment {
        input: 
            input_fastq1 = FindIntersectingSnpsPairedEnd.fastq1,
            input_fastq2 = FindIntersectingSnpsPairedEnd.fastq2,
            read_group = read_group,
            reference_fasta = reference_fasta,
            references = references, 
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
            output_file_name = sample_name + ".keep.merge.sort.bam"
    }

    call IndexBam {
        input: 
            sorted_bam = SortBam.sorted_bam, 
            output_file_name = sample_name + ".keep.merge.sort.bam.bai"
    }

    output {
        File snp_index = compressVcf.snp_index
        File snp_tab = compressVcf.snp_tab
        File haplotype = compressVcf.haplotype
        File sorted_bam = SortBam.sorted_bam
        File sorted_bam_index = IndexBam.bai
    }
}

# Preliminary task: compress VCF file with bgzip 
# edited to combine the compress + split tasks 
task compressVcf {
    input {
        File vcf
        File vcf_index
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
#         docker: "apariciobioinformaticscoop/wasp-mapping:sept12"
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
#         docker: "apariciobioinformaticscoop/wasp-mapping:sept12"
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
        String sample_name
        String sa_id_tumor
        String sa_id_normal
        String output_txt_file = "samples.txt"
    }

    # command <<<
    #     echo ~{sample_name} > ~{output_txt_file}
    #     grep "\S" ~{output_txt_file}
    # >>>

    # heredoc command: 
    command <<<
        echo ~{sa_id_normal} > ~{output_txt_file}
        echo ~{sa_id_tumor} >> ~{output_txt_file}
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
        File input_bai
        File snp_index
        File snp_tab
        File haplotype
        File samples_txt
        String fastq1_name
        String fastq2_name
        String keep_bam
        String remap_bam
        String output_dir
        File intersecting_snps_script
    }

    Float multiplier = 15
    Int disk_size = ceil(size(input_bam, "GB") * multiplier) + 100 # also bumped up runtime metrics

    # BAM_PATH=$(pwd ~{input_bam})/$(ls ~{input_bam})
    # echo $BAM_PATH
    
    command <<<
        python ~{intersecting_snps_script} \
        --is_paired_end \
        --is_sorted \
        --snp_tab ~{snp_tab} \
        --snp_index ~{snp_index} \
        --samples ~{samples_txt} \
        --haplotype ~{haplotype} ~{input_bam}
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/wasp-mapping:sept12"
        disk: disk_size + " GB" 
        cpu: 40 # changed from 24
        memory: "64 GB" # changed from 24
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
        ReferenceFasta references
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
        docker: "apariciobioinformaticscoop/wasp-mapping:sept12"
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
    Float sort_sam_disk_multiplier = 4
    Int disk_size = ceil(sort_sam_disk_multiplier * size(merged_bam, "GB")) + 60

    command <<<
        samtools sort -o ~{output_file_name} ~{merged_bam} -@ 4
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/wasp-mapping:latest"
        disk: disk_size + " GB"
        cpu: 14
        memory: "24 GB"
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
        samtools index -b ~{sorted_bam} ~{output_file_name}
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


