import numpy as np
import pandas as pd

###############################################################################################################
### Constants ###
###############################################################################################################

plasmabafmean = 0.5
lowthreshbafplasma = 0.4
highthreshbafplasma = 0.4

normalbafmean = 0.5
lowthreshbafnormal = 0.3
highthreshbafnormal = 0.3


###############################################################################################################
### Functions ###
###############################################################################################################

def plasma_baf_filter(df_cnv):
    """ Filters input df by plasma BAF value

    df_cnv: input dataframe
    """

    return df_cnv.loc[np.logical_and((df_cnv['maj_all_baf_plasma'] < (plasmabafmean + highthreshbafplasma)), (df_cnv['maj_all_baf_plasma'] > (plasmabafmean - lowthreshbafplasma)))].reset_index(drop=True)


def normal_baf_filter(df_cnv):
    """ Filters input df by normal BAF value

    df_cnv: input dataframe
    """

    return df_cnv.loc[np.logical_and((df_cnv['maj_all_baf_normal'] < (normalbafmean + highthreshbafnormal)), (df_cnv['maj_all_baf_normal'] > (normalbafmean - lowthreshbafnormal)))].reset_index(drop=True)


def homozygous_plasma_filter(df_cnv):
    """ Removes regions where baf is completely homozygous in plasma sample

    df_cnv: input dataframe
    """

    return df_cnv.loc[~(df_cnv['maj_all_baf_plasma'] == 1) & ~(df_cnv['maj_all_baf_plasma'] == 0)]


def homozygous_normal_filter(df_cnv):
    """ Removes regions where baf is completely homozygous in normal sample

    df_cnv: input dataframe
    """

    return df_cnv.loc[~(df_cnv['maj_all_baf_normal'] == 1) & ~(df_cnv['maj_all_baf_normal'] == 0)]


def filter_primary(df_cnv):
    """ Filter input dataframe if a matching normal is present

    df_cnv: input dataframe
    """
    df_cnv = homozygous_normal_filter(df_cnv)
    df_cnv = homozygous_plasma_filter(df_cnv)

    df_cnv = plasma_baf_filter(df_cnv)
    df_cnv = normal_baf_filter(df_cnv)

    return df_cnv


def filter_no_primary(df_cnv):
    """ Filter input dataframe if no matching normal is present

    df_cnv: input dataframe
    """

    df_cnv = homozygous_plasma_filter(df_cnv)
    df_cnv = plasma_baf_filter(df_cnv)
    return df_cnv


def filter_bafs(df_cnv, is_primary=False):
    """ Filters input dataframe on several parameters (plasma and normal baf values for instance)

    df_cnv: input dataframe
    """

    return filter_primary(df_cnv) if is_primary else filter_no_primary(df_cnv)
