#!/usr/bin/env python3

########################################
# Daniel Halmos, Minita Shah
# New York Genome Center
########################################

import os
import sys
import argparse
import pandas as pd
import numpy as np
from math import floor
from joint_filters import *

###############################################################################################################
### Classes and Functions ###
###############################################################################################################
global chunksize

class ArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        self.print_help(sys.stderr)
        self.exit(2, '\nERROR: %s\n\n' % (message))


def alt_agg(df):
    """ Aggregate

    df: input dataframe
    """
    chunksize = 50

    bafs = pd.DataFrame()

    for cn in list(set(df['win_chr'].astype(str) + ' ' + df['Total_cn'].astype(str) + ' ' + df['Minor_cn'].astype(str))):
        piece = pd.DataFrame()

        segs = df.loc[df['win_chr'].astype(
            str) + ' ' + df['Total_cn'].astype(str) + ' ' + df['Minor_cn'].astype(str) == cn]

        if len(segs['maj_all_cnt_plasma']) >= chunksize:
            cnks = np.array_split(segs[['maj_all_cnt_plasma', 'ref_cnt', 'alt_cnt', 'win_chr']], floor(
                len(segs['maj_all_cnt_plasma'])/chunksize))

            piece['maj_all_cnt_plasma'] = [
                np.sum(cnk['maj_all_cnt_plasma']) for cnk in cnks]
            piece['ref_cnt_plasma'] = [np.sum(cnk['ref_cnt']) for cnk in cnks]
            piece['alt_cnt_plasma'] = [np.sum(cnk['alt_cnt']) for cnk in cnks]

            piece['chunk_len'] = [len(cnk) for cnk in cnks]

            piece['maj_all_baf_plasma'] = piece['maj_all_cnt_plasma'] / \
                (piece['ref_cnt_plasma'] + piece['alt_cnt_plasma'])

            piece['Total_cn'] = segs.iloc[0]['Total_cn']
            piece['Minor_cn'] = segs.iloc[0]['Minor_cn']
            piece['win_chr'] = segs.iloc[0]['win_chr']

            bafs = bafs._append(piece)


    return bafs


def read_arguments() -> argparse.Namespace:
    """ Reads the arguments passed to the script.
    """
    parser = ArgumentParser(prog='analyze_baf_wins.py',
                            description='', epilog='')
    parser.add_argument('-c', '--compiled-baf',
                        help='Pre-Filtered Merged BAF File', required=True)
    parser.add_argument('-o', '--out-prefix',
                        help='Output directory.', required=True)

    args = parser.parse_args()

    return args


###############################################################################################################
### Main ###
###############################################################################################################


def main():
    """ Runs the whole script.
    """

    args = read_arguments()

    df_cnv_file = args.compiled_baf
    outprefix = args.out_prefix
    
    chunksize = 50

    df_cnv = pd.read_csv(df_cnv_file, sep="\t", header=0)

    df_cnv['Minor_cn'] = df_cnv.apply(lambda x: 0 if (
        x['Total_cn'] <= 2 or x['Minor_cn'] == 0) else 1, axis=1)

    df_cnv = filter_bafs(df_cnv, True if 'tum' in df_cnv_file else False)

    alt_wins = alt_agg(df_cnv)

    out_file = os.path.join(outprefix + ".cnv.grp.snpconst.txt")
    if len(alt_wins) > 0:
        alt_wins.to_csv(out_file, sep='\t', index=False)


if __name__ == '__main__':
    main()
