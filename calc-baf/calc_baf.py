#!/nfs/sw/python/python-2.7.8/bin/python

########################################
# Minita Shah (mshah@nygenome.org)
# New York Genome Center
########################################

import os
import sys
import argparse
import re
import pandas as pd
import numpy as np
import gzip
import io

########## Classes and Functions ##########

class ArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        self.print_help(sys.stderr)
        self.exit(2, '\nERROR: %s\n\n' % (message))

################### Main ###################
parser = ArgumentParser(prog='parsefragcounter.robustZscore.py', description='', epilog='')
parser.add_argument('-i', '--inp-dir', help='Input directory', required=True)
parser.add_argument('-o', '--out-dir', help='Output Directory', required=True)
parser.add_argument('-he', '--het-vcf', help='Gzipped het vcf', required=True)
parser.add_argument('-t', '--tumor-name', help='Tumor name', required=True)
parser.add_argument('-n', '--normal-name', help='Normal name', required=True)
parser.add_argument('-v', '--version', help='HaplotypeCaller Version', required=True)

args = parser.parse_args()
inp_dir = args.inp_dir
out_dir = args.out_dir
het_vcf = args.het_vcf
tumor_name = args.tumor_name
normal_name = args.normal_name
version = args.version

with gzip.open(het_vcf, 'rb') as f:
    lines = [l for l in f if not l.startswith(b'##')]
df_vcf = pd.read_csv(
    io.BytesIO(b''.join(lines)),
    dtype={'#CHROM': str, 'POS': str, 'ID': str, 'REF': str, 'ALT': str,
        'QUAL': str, 'FILTER': str, 'INFO': str, 'FORMAT': str, normal_name: str},
    sep='\t'
).rename(columns={'#CHROM': 'CHROM'})

#print(df_vcf.count()[0])

inp_file = os.path.join(inp_dir, tumor_name+"--"+normal_name+".haplotypeCalls." + version + ".het.alleles.txt")
df_allele = pd.read_csv(inp_file, header=0, sep="\t")
df_allele.columns = ['CHROM', 'POS', 'REF', 'RD', 'A', 'C', 'G', 'T', 'STRAND']
df_allele['CHROM'] = df_allele.CHROM.astype(str)
df_allele['POS'] = df_allele.POS.astype(str)

#print(df_allele.count()[0])

# Merge het vcf and sample het alleles 
df_combined = pd.merge(df_vcf, df_allele, on=['CHROM', 'POS'])
#print(df_combined.head())

# Get ref, alt allele count
for i in ['A', 'C', 'G', 'T']:
    df_combined.loc[df_combined['REF_x']==i,'REF_COUNT'] = df_combined[i].astype(float)
    df_combined.loc[df_combined['ALT']==i,'ALT_COUNT'] = df_combined[i].astype(float)

#print(df_combined.head())

# Calc BAF
df_combined.loc[(df_combined['ALT_COUNT']==0) & (df_combined['REF_COUNT']==0), 'BAF'] = 0.0
df_combined.loc[(df_combined['ALT_COUNT']!=0) | (df_combined['REF_COUNT']!=0), 'BAF'] = df_combined['ALT_COUNT']/(df_combined['ALT_COUNT'] + df_combined['REF_COUNT'])

# Print output
df_combined.rename(columns={'REF_x':'REF'}, inplace=True)
out_file = os.path.join(out_dir, tumor_name+"--"+normal_name+".haplotypeCalls." + version + ".het.baf.txt")
df_combined.to_csv(out_file, sep='\t', index=False, columns=['CHROM', 'POS', 'REF', 'ALT', 'REF_COUNT', 'ALT_COUNT', 'BAF'])
