version 1.0

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
    }

    call splitVcfByChr {
        input: 
            vcf = input_vcf, 
            sample_name = sample_name
    }

    call VcfToH5 {
        input: 
            chrom_info = chrom_info, 
            vcfs = splitVcfByChr.output_vcfs, 
            snp_index_name = sample_name + ".snp_index.h5",
            snp_tab_name = sample_name + ".snp_tab.h5", 
            haplotype_name = sample_name + ".h5"
    }

    call MakeSampleNameTxt {
        input: 
            sample_name = sample_name
    }

    call FindIntersectingSnpsPairedEnd {
        input: 
            input_bam = input_bam, 
            snp_index = VcfToH5.snp_index, 
            snp_tab = VcfToH5.snp_tab, 
            haplotype = VcfToH5.haplotype, 
            sample_txt = MakeSampleNameTxt.txt_file, 
            fastq1_name = sample_name + ".to.remap.fq1.gz", 
            fastq2_name = sample_name + ".to.remap.fq2.gz", 
            keep_bam = sample_name + ".keep.bam", 
            remap_bam = sample_name + ".to.remap.bam"
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
            output_file_name = sample_name + ".kept.bam"
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

# split pre-processing outputted VCF into multiple VCFs, sorted by chromosome 
# A. CHECK: syntax of the regex after --regions against the mrdedge example. 
task splitVcfByChr {
    input {
        File vcf
        String sample_name
    }

    Int disk_size = ceil(size(vcf, "GB") * 4)

    command <<<
        mkdir vcf_by_chr
        for i in {1..22} 
        do
            bcftools view ~{vcf} --regions chr$i > vcf_by_chr/~{sample_name}.$i.het.vcf.gz 
        done
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
        Array[File] output_vcfs = "/vcf_by_chr/*.het.vcf.gz" # unsure of this
    }
}

# B. Convert input VCF to a HDF5 SNP file
task VcfToH5 {
    input {
        File chrom_info # stored as a reference in inputs. 
        Array[File] vcfs
        String snp_index_name
        String snp_tab_name 
        String haplotype_name
    }

    Int disk_size = ceil(size(vcfs, "GB") * 2)

    command <<<
        python snp2h5.py \
        --chrom ~{chrom_info} \
        --format vcf \
        --snp_index ~{snp_index_name} \
        --snp_tab ~{snp_tab_name} \
        --haplotype ~{haplotype_name} \
        ~{vcfs}
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/wasp-mapping:latest"
        disk: disk_size + " GB"
        cpu: 16
        memory: "16 GB"
        preemptible: true
        maxRetries: 3
    }

    output {
        File snp_index = "~{snp_index_name}"
        File snp_tab = "~{snp_tab_name}"
        File haplotype = "~{haplotype_name}"
    }
}

task MakeSampleNameTxt {
    input {
        String sample_name
        String output_txt_file = "sample_names.txt"
    }

    command <<<
        echo ~{sample_name} > ~{output_txt_file}
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/wasp-mapping:latest"
        disk: "10 GB"
        cpu: 4
        memory: "4 GB"
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
        File sample_txt 
        String fastq1_name
        String fastq2_name
        String keep_bam
        String remap_bam
        String output_dir = "find_intersecting_snps"
    }

    Float multiplier = 2.5
    Int disk_size = ceil(size(input_bam, "GB") * multiplier) + 50

    command <<<
        python find_intersecting_snps.py \
        --is_paired_end \
        --is_sorted \
        --output_dir ~{output_dir} \
        --snp_tab ~{snp_tab} \
        --snp_index ~{snp_index} \
        --haplotype ~{haplotype} \
        --samples ~{sample_txt} \ 
        ~{input_bam}
    >>>

    runtime {
        docker: "apariciobioinformaticscoop/wasp-mapping:latest"
        disk: disk_size + " GB" 
        cpu: 16
        memory: "16 GB"
        preemptible: true
        maxRetries: 3
    }


    output {
        File fastq1 = fastq1_name
        File fastq2 = fastq2_name
        File keep_bam = keep_bam # reads that didn't intersect with SNPs
        File remap_bam = remap_bam # reads that intersected with SNPs, needs to be flipped + remapped
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
    }

    Float multiplier = 2.5
    Int disk_size = ceil(size(realigned_bam, "GB") * multiplier) + 50
    

    command <<<
        python filter_remapped_reads.py \
        ~{remapped_bam} \
        ~{realigned_bam} \
        ~{output_file_name}
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
        File sorted_bam = "~{output_file_name}"
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
        File bai = "~{output_file_name}"
    }
}


