version 1.0

import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/ubam-pre-pro/tasks/GermlineStructs.wdl"
import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/sv-caller/paired/PairedGridss.wdl" as Gridss
import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/sv-caller/paired/PairedManta.wdl" as Manta
import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/sv-caller/paired/PairedLumpy.wdl" as Lumpy
import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/sv-caller/paired/PairedNameSort.wdl" as Sort
import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/sv-caller/Viola.wdl" as Viola
import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/sv-caller/paired/PairedSurvivor.wdl" as Survivor

workflow PairedSVCaller {

    input {
        String study

        Boolean exclude_regions

		String tumor_name
		File tumor_bam
		File tumor_bai

		String normal_name
		File normal_bam
		File normal_bai

        ReferenceFasta references

        File viola_script_m
        File viola_script_l
        File viola_script_g
    }

    call Sort.PairedNameSort {
        input:
            tumor_name = tumor_name,
            tumor_bam = tumor_bam,
            tumor_bai = tumor_bai,
            normal_name = normal_name,
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            references = references
    }

    call Lumpy.PairedLumpy {
        input:
            exclude_regions = exclude_regions,
            tumor_name = tumor_name,
            tumor_bam = PairedNameSort.out_tumor_bam,
            tumor_bai = tumor_bai,
            normal_name = normal_name,
            normal_bam = PairedNameSort.out_normal_bam,
            normal_bai = normal_bai,
            references = references
    }

    call Manta.PairedManta {
        input:
            tumor_name = tumor_name,
            tumor_bam = tumor_bam,
            tumor_bai = tumor_bai,
            normal_name = normal_name,
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            references = references
    }

    call Gridss.PairedGridss {
        input:
            exclude_regions = exclude_regions,
            tumor_name = tumor_name,
            tumor_bam = tumor_bam,
            tumor_bai = tumor_bai,
            normal_name = normal_name,
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            references = references
    }

    call Viola.Viola {
        input:
            manta_input = PairedManta.filtered,
            lumpy_input  = PairedLumpy.filtered,
            gridss_input = PairedGridss.filtered,
            viola_script_m = viola_script_m,
            viola_script_l = viola_script_l,
            viola_script_g = viola_script_g
    }

    call Survivor.PairedSurvivor {
        input:
            manta_input = Viola.manta_out,
            lumpy_input = Viola.lumpy_out,
            gridss_input = Viola.gridss_out,
            tumor_name = tumor_name,
            normal_name = normal_name
    }

    output {
        File lumpy_filtered = PairedLumpy.filtered
        File manta_filtered = PairedManta.filtered
        File gridss_filtered = PairedGridss.filtered
        File lumpy_unfiltered = PairedLumpy.unfiltered
        File manta_unfiltered = PairedManta.unfiltered
        File gridss_unfiltered = PairedGridss.unfiltered
        File merged1 = PairedSurvivor.merged1
        File merged2 = PairedSurvivor.merged2
        File merged3 = PairedSurvivor.merged3
    }
}
