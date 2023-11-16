import argparse
import os
import shutil
import viola

def parse_options():
    parser = argparse.ArgumentParser(description="runs viola.py on the caller outputs")
    
    parser.add_argument("vcf_path", help="vcf path")

    return parser.parse_args()

def viola_run(vcf_path):
    vcf_in = str(vcf_path)
    vcf_org = vcf_in + '.org'
    vcf_out = 'gridss_out.vcf'

    os.rename(vcf_in, vcf_org)
    with open(vcf_in, 'w') as new:
        with open(vcf_org, 'r') as org:
            for line in org:
                # FIX INFO field: change PARID to MATEID
                new.write(line.replace('PARID', 'MATEID'))

    try:
        sv = viola.read_vcf(vcf_in, variant_caller='gridss').breakend2breakpoint()
        sv.to_vcf(vcf_out)
    except Exception:
        shutil.copyfile(vcf_in, vcf_out)       

if __name__ == "__main__":
    options = parse_options()
    viola_run(options.vcf_path)
