version 1.0

import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/sv-caller/Viola.wdl" as Viola
import "https://raw.githubusercontent.com/aparicio-bioinformatics-coop/cromwell-workflows/main/sv-caller/paired/PairedSurvivor.wdl" as Survivor

workflow SVCleanup {

    input {
        String study

        File manta
        File lumpy
        File gridss

		String tumor_name
        String normal_name


        File viola_script_m
        File viola_script_l
        File viola_script_g
    }

    call Viola.Viola {
        input:
            manta_input = manta,
            lumpy_input  = lumpy,
            gridss_input = gridss,
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
        File merged1 = PairedSurvivor.merged1
        File merged2 = PairedSurvivor.merged2
        File merged3 = PairedSurvivor.merged3
    }
}
