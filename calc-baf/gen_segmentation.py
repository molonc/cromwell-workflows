#!/usr/bin/env python3

########################################
# Converts ascat/sequenza segmentation to bed
# Minita Shah (mshah@nygenome.org)
# New York Genome Center
########################################

import os
import sys
import argparse
import re
import pandas as pd
import numpy as np
import pyranges as pr

########## Classes and Functions ##########

class ArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        self.print_help(sys.stderr)
        self.exit(2, '\nERROR: %s\n\n' % (message))

def mask_chroms(df):
    chrs = ['chr{}'.format(i) for i in range(1, 23) if i != 19]
    #mask = df.Chromosome.apply(lambda x: any(item in x for item in chrs))
    mask = df.Chromosome.isin(chrs)
    df = df[mask]
    return df

def read_cnv_seg_ascat(inp):
    #Chromosome Start   End SegmentedLogR   MinorCopyNumber MinorAllele
    #1   13116   638771  0.6000449017429913  4   0
    df = pd.read_csv(inp, sep="\t", header=0)
    df['Chromosome'] = df['Chromosome'].apply(lambda x: f"chr{x}")
    
    df.columns = ['Chromosome', 'Start', 'End', 'TN_ratio', 'Total_cn', 'Minor_cn']
    df = mask_chroms(df)
    df['Start'] = df['Start']-1
    df['Call'] = '-'
    df['Size'] = (df['End']-df['Start'])/float(1000000)

    return df

def read_cnv_seg_ascat_mean(inp):
    #Chromosome Start   End SegmentedLogR   SegmentedBAF    SegmentedLogR_trans mean_transformed    Call
    #chr1    13116   3093405 -0.0249493962849521 0.5 0.9828550723906235  2.198007527744699   NEU
    df = pd.read_csv(inp, sep="\t", header=0)
    df = mask_chroms(df)
    df['Start'] = df['Start']-1
    df = df[['Chromosome', 'Start', 'End', 'SegmentedLogR_trans', 'mean_transformed', 'SegmentedBAF', 'Call', 'Minor_cn']]
    df.rename(columns={'SegmentedLogR_trans': 'TN_ratio', 'mean_transformed':'Total_cn', 'Call':'Call_by_TN_ratio'}, inplace=True)
    df['Call'] = '-'
    df['Size'] = (df['End']-df['Start'])/float(1000000)
    df['LOH'] = np.where(((df['Minor_cn']==0) & (df['Size']>1.5)), 'Yes', 'NA')
    df['Minor_cn'] = np.where((df['Minor_cn']!=0), 'NA', df['Minor_cn'])
    return df


def read_cnv_seg_sequenza(inp):
    # #chromosome   start.pos   end.pos Bf  N.BAF   sd.BAF  depth.ratio N.ratio sd.ratio    CNt A   B   LPP
    df = pd.read_csv(inp, sep="\t", header=0)
    df = df.drop(['Bf', 'N.BAF', 'sd.BAF', 'N.ratio', 'sd.ratio', 'A', 'LPP'], axis=1)
    df.columns = ['Chromosome', 'Start', 'End', 'TN_ratio', 'Total_cn', 'Minor_cn']
    df['Call'] = '-'

    df = df[['Chromosome', 'Start', 'End', 'Total_cn', 'Minor_cn', 'Call', 'TN_ratio']]

    df = mask_chroms(df)
    df['Start'] = df['Start']-1
    df['Size'] = (df['End']-df['Start'])/float(1000000)

    return df

def read_cnv_seg_ichor(inp):
    #AD-05_A/ichor_filtered/AD-05_A.seg.txt
    #ID      chrom   start   end     num.mark        seg.median.logR copy.number     call    subclone.status logR_Copy_Number        Corrected_Copy_Number   Corrected_Call
    df = pd.read_csv(inp, sep="\t", header=0)
    df = df.drop(['ID', 'num.mark', 'seg.median.logR', 'subclone.status', 'logR_Copy_Number', 'Corrected_Copy_Number', 'Corrected_Call'], axis=1)
    df['Minor_cn'] = "-"
    df.columns = ['Chromosome', 'Start', 'End', 'Total_cn', 'Call', 'Minor_cn']
    df['TN_ratio'] = '-'

    col = 'Call'
    conditions = [ df[col] == "NEUT", ((df[col] == 'GAIN') | (df[col] == 'AMP') | (df[col] == 'HLAMP')), df[col] == 'HETD' ]
    choices = [ "NEU", 'DUP', 'DEL' ]
    
    df['Call_by_TN_ratio'] = np.select(conditions, choices, default=np.nan)

    df = mask_chroms(df)
    df['Start'] = df['Start']-1
    df['Size'] = (df['End']-df['Start'])/float(1000000)
    df['LOH'] = '-'
    df['SegmentedBAF'] = '-'

    df = df[['Chromosome', 'Start', 'End', 'TN_ratio', 'Total_cn', 'SegmentedBAF', 'Call_by_TN_ratio', 'Minor_cn', 'Call', 'Size', 'LOH']]

    return df
    
    
################### Main ###################
parser = ArgumentParser(prog='ascat_sequenza_to_bed.py', description='', epilog='')
parser.add_argument('-i', '--inp-file', help='Input ascat/sequenza/ichor/bicseq segmentation.', required=True)
parser.add_argument('-o', '--out-file', help='Output segmentation file', required=True)
parser.add_argument('-t', '--tool', help='ascat/ascat_mean/sequenza/ichor/bicseq', required=True)

args = parser.parse_args()
inp_file = args.inp_file
out_file = args.out_file
tool = args.tool

pd.set_option('display.max_columns', None)

## Read segmentation file
if tool == 'ascat':
    df = read_cnv_seg_ascat(inp_file)
elif tool == 'ascat_mean':
    df = read_cnv_seg_ascat_mean(inp_file)
elif tool == 'sequenza':
    df = read_cnv_seg_sequenza(inp_file)
elif tool == 'ichor':
    df = read_cnv_seg_ichor(inp_file)

## Print df to file
## <sample_name>.segmentation.txt
df.to_csv(out_file, sep='\t', index=False, na_rep='NA')
