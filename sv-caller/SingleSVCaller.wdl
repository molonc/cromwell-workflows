version 1.0

import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/ubam-pre-pro/tasks/GermlineStructs.wdl"
import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/sv-caller/single/SingleGridss.wdl" as Gridss
import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/sv-caller/single/SingleManta.wdl" as Manta
import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/sv-caller/single/SingleLumpy.wdl" as Lumpy
import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/sv-caller/single/SingleNameSort.wdl" as Sort
import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/sv-caller/Viola.wdl" as Viola
import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/sv-caller/single/SingleSurvivor.wdl" as Survivor

workflow SingleSVCaller {

    input {
        String study

        Boolean exclude_regions

		String tumor_name
		File tumor_bam
		File tumor_bai

        ReferenceFasta references

        File viola_script_m
        File viola_script_l
        File viola_script_g
    }

    call Sort.SingleNameSort {
        input:
            tumor_name = tumor_name,
            tumor_bam = tumor_bam,
            tumor_bai = tumor_bai,
            references = references
    }

    call Lumpy.SingleLumpy {
        input:
            exclude_regions = exclude_regions,
            tumor_name = tumor_name,
            tumor_bam = SingleNameSort.out_tumor_bam,
            tumor_bai = tumor_bai,
            references = references
    }

    call Manta.SingleManta {
        input:
            tumor_name = tumor_name,
            tumor_bam = tumor_bam,
            tumor_bai = tumor_bai,
            references = references
    }

    call Gridss.SingleGridss {
        input:
            exclude_regions = exclude_regions,
            tumor_name = tumor_name,
            tumor_bam = tumor_bam,
            tumor_bai = tumor_bai,
            references = references
    }

    call Viola.Viola {
        input:
            manta_input = SingleManta.filtered,
            lumpy_input  = SingleLumpy.filtered,
            gridss_input = SingleGridss.filtered,
            viola_script_m = viola_script_m,
            viola_script_l = viola_script_l,
            viola_script_g = viola_script_g
    }

    call Survivor.SingleSurvivor {
        input:
            manta_input = Viola.manta_out,
            lumpy_input = Viola.lumpy_out,
            gridss_input = Viola.gridss_out,
            tumor_name = tumor_name
    }

    output {
        File lumpy_filtered = SingleLumpy.filtered
        File manta_filtered = SingleManta.filtered
        File gridss_filtered = SingleGridss.filtered
        File lumpy_unfiltered = SingleLumpy.unfiltered
        File manta_unfiltered = SingleManta.unfiltered
        File gridss_unfiltered = SingleGridss.unfiltered
        File merged1 = SingleSurvivor.merged1
        File merged2 = SingleSurvivor.merged2
        File merged3 = SingleSurvivor.merged3
    }
}
