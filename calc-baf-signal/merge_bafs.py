#!/usr/bin/env python3

########################################
# Daniel Halmos
# New York Genome Center
########################################

import os
import sys
import argparse
import pandas as pd
import numpy as np

###############################################################################################################
### Classes and Functions ###
###############################################################################################################


class ArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        self.print_help(sys.stderr)
        self.exit(2, '\nERROR: %s\n\n' % (message))


def read_arguments() -> argparse.Namespace:
    """ Reads the arguments passed to the script.
    """
    parser = ArgumentParser(prog='analyze_baf_wins.py',
                            description='', epilog='')
    parser.add_argument('-p', '--plasma-tumor-baf',
                        help='Plasma tumor BAF file.', required=True)
    parser.add_argument('-t', '--tumor-baf',
                        help='Tumor BAF file to get regions.', required=True)
    parser.add_argument('-n', '--normal-baf',
                        help='normal BAF file to get regions.', required=True)
    parser.add_argument('-o', '--out-prefix',
                        help='Output directory.', required=True)
    parser.add_argument(
        '-g', '--ogfile', help='Original File if using WASP', required=False)

    args = parser.parse_args()

    return args

###############################################################################################################
### Main ###
###############################################################################################################


def main():
    """ Runs the whole script.
    """

    args = read_arguments()

    plasma_tumor_baf = args.plasma_tumor_baf
    tumor_baf = args.tumor_baf
    normal_baf = args.normal_baf
    outprefix = args.out_prefix
    ogfile = args.ogfile

    # Setting parameters

    global chunksize
    chunksize = 50

    cols = ['chr', 'start', 'end', 'ref', 'alt', 'hc_qual', 'hc_filter', 'hc_info', 'hc_format', 'hc_quals', 'win_chr', 'win_start', 'win_end', 'Total_cn', 'Minor_cn', 'Size', 'overlap']

    # Loading dataframes
    df_plasma_tumor = pd.read_csv(plasma_tumor_baf, sep="\t", header=0)

    if normal_baf != 'None':
        
        df_normal = pd.read_csv(normal_baf, sep="\t", header=0)
        # WASP filtration
        # We load a WASP-filtered file and discard all the SNPs where the WASP filtered-count and the no filtration count differ
        df_normal_og = pd.read_csv(ogfile, sep="\t", header=0)

        df_normals = pd.merge(df_normal, df_normal_og, how='left', on=['chr', 'start', 'end', 'ref', 'alt', 'hc_qual', 'hc_filter', 'hc_info', 'hc_format', 'hc_quals', 'win_chr', 'win_start', 'win_end', 'Total_cn', 'Minor_cn', 'Size'] if 'AA' in outprefix else [
            'chr', 'start', 'end', 'ref', 'alt', 'hc_qual', 'hc_filter', 'hc_info', 'hc_format', 'hc_quals', 'win_chr', 'win_start', 'win_end', 'Total_cn', 'Minor_cn', 'Size', 'overlap'], suffixes=['', '_OG'])

        # We filter df_normal
        df_normal = df_normal.loc[np.abs(df_normals['ref_cnt_OG'] + df_normals['alt_cnt_OG']) -
                                  np.abs(df_normals['ref_cnt'] + df_normals['alt_cnt']) < 1]
    
    else:
        df_plasma_tumor_og = pd.read_csv(ogfile, sep="\t", header=0)
        df_plasmas = pd.merge(df_plasma_tumor, df_plasma_tumor_og, how='left', on=['chr', 'start', 'end', 'ref', 'alt', 'hc_qual', 'hc_filter', 'hc_info', 'hc_format', 'hc_quals', 'win_chr', 'win_start', 'win_end', 'Total_cn', 'Minor_cn', 'Size', 'overlap'], suffixes=['', '_OG'])
        
        df_plasma_tumor = df_plasma_tumor.loc[np.abs(df_plasmas['ref_cnt_OG'] + df_plasmas['alt_cnt_OG']) -
                                  np.abs(df_plasmas['ref_cnt'] + df_plasmas['alt_cnt']) < 1]
        
        

    if tumor_baf != 'None':
        df_tumor = pd.read_csv(tumor_baf, sep="\t", header=0)
        
        
        df_tumor['maj_all_col'] = df_tumor.apply(lambda x: np.where(
        x['ref_cnt'] > x['alt_cnt'], 'ref_cnt', 'alt_cnt'), axis=1)
        df_tumor['maj_all_baf'] = np.where(df_tumor['maj_all_col'] == 'ref_cnt', (df_tumor['ref_cnt']/(
        df_tumor['ref_cnt']+df_tumor['alt_cnt'])), (df_tumor['alt_cnt']/(df_tumor['ref_cnt']+df_tumor['alt_cnt'])))
        

        # Merge plasma and tumor data frames
        df_combined = pd.merge(df_plasma_tumor, df_tumor,
                           on=cols, suffixes=("", "_tumor"))

    # Get counts and Baf of major allele (in according to biopsy) in plasma

    # Plasma maj_all cnt and baf -- First part of tumor-informed
    else:
        df_plasma_tumor['alt_ref'] = np.random.randint(0, 2, df_plasma_tumor.shape[0])
        df_plasma_tumor['maj_all_col'] = np.where(df_plasma_tumor.alt_ref==1, 'ref_cnt', 'alt_cnt')
        df_combined = df_plasma_tumor
    
    df_combined['maj_all_cnt_plasma'] = np.where(
        df_combined['maj_all_col'] == 'ref_cnt', df_combined['ref_cnt'], df_combined['alt_cnt'])
    df_combined['maj_all_baf_plasma'] = np.where(df_combined['maj_all_col'] == 'ref_cnt', (df_combined['ref_cnt']/(
        df_combined['ref_cnt']+df_combined['alt_cnt'])), (df_combined['alt_cnt']/(df_combined['ref_cnt']+df_combined['alt_cnt'])))

    if normal_baf != 'None':
        # Merge Combined with plasma
        df_combined = pd.merge(df_combined, df_normal, on=cols,
                           suffixes=("", "_normal"))

        # Get counts and Baf of major allele in normal
        df_combined['maj_all_cnt_normal'] = np.where(
        df_combined['maj_all_col'] == 'ref_cnt_normal', df_combined['ref_cnt_normal'], df_combined['alt_cnt_normal'])
        df_combined['maj_all_baf_normal'] = np.where(df_combined['maj_all_col'] == 'ref_cnt_normal', (df_combined['ref_cnt_normal']/(
        df_combined['ref_cnt_normal']+df_combined['alt_cnt_normal'])), (df_combined['alt_cnt_normal']/(df_combined['ref_cnt_normal']+df_combined['alt_cnt_normal'])))

    else:
        pass
    
    df_combined['plasma_cov'] = df_combined['ref_cnt']+df_combined['alt_cnt']

    # Get copy number rounded to int (from tumor biopsy)
    df_combined['Total_cn'] = np.round(pd.to_numeric(
        df_combined['Total_cn'], errors='coerce'))
    df_combined['Minor_cn'] = np.round(pd.to_numeric(
        df_combined['Minor_cn'], errors='coerce'))

    df_combined['Minor_cn'] = df_combined['Minor_cn'].fillna(1)
    df_combined['Minor_cn'] = df_combined.apply(lambda x: 0 if (x['Total_cn'] < 2) else x['Minor_cn'], axis = 1)

    # Get CNV regions (either DEL/DUP or CNLOH)
    cnvregion = np.logical_and(((df_combined['Total_cn'] != 2) | ((df_combined['Total_cn'] == 2) & (
        df_combined['Minor_cn'] == 0))), df_combined['Total_cn']/df_combined['Minor_cn'] != 2)

    # Second part of tumor-informed -- we keep only the CNV bearing regions
    df_cnv = df_combined.loc[cnvregion]
    df_cnv = df_cnv.reset_index(drop=True)

    # This is the file used downstream
    out_file = os.path.join(outprefix + f".{'ctrl' if normal_baf == 'None' else 'tum'}.cnv.unfilt.txt")
    df_cnv.to_csv(out_file, sep='\t', index=False)


if __name__ == '__main__':
    main()
