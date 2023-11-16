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
    vcf_out = 'lumpy_out.vcf'

    try:
        sv = viola.read_vcf(vcf_in, variant_caller='lumpy').breakend2breakpoint()
        sv.to_vcf(vcf_out)
    except Exception:
        shutil.copyfile(vcf_in, vcf_out)       

if __name__ == "__main__":
    options = parse_options()
    viola_run(options.vcf_path)
