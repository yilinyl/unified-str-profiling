#!/usr/bin/env python3

import sys
import vcf
import subprocess
import pandas as pd
import Select as st
from Bio.PopGen.GenePop.EasyController import EasyController
from decimal import *
import time

def rul_filter(df_raw, rul='2_3_4'):
    """
    Initial selection of STR loci: apply repeat unit length filter
    :param df_raw:
    :param rul:
    :return:
    """
    print('========================================================================================')
    print('Processing begins')
    print('Total loci:', len(df_raw))
    df_raw = df_raw[~df_raw.CHROM.isin(['chrX', 'chrY'])]  # drop chrX and chrY
    print('========================================================================================')
    print('After drop chrX and chrY, remaining:', len(df_raw))
    df_raw = df_raw[df_raw.RUL.isin([int(i) for i in rul.split('_')])]  # select repeat length
    df_raw = df_raw.sort_values(by=['PD', 'GT_FREQ'], ascending=[False, False]).reset_index(drop=True)
    return df_raw

def para_filter(df, info, filter_type='default'):
    """

    :param df:
    :param info:
    :param filter_type:
    :return:
    """
    # loose filter
    if filter_type == 'loose':
        df = df[df.EXP_HE >= crit.MIN.HE]
        print('=======================================================================================')
        print('After filtering with HE, remaining:', len(df))
        info += len(df),

        df = df[df.PM <= crit.MAX.PM]
        print('=======================================================================================')
        print('After filtering with PM, remaining:', len(df))
        info += len(df),

        df = df[df.PE >= crit.MIN.PE]
        print('=======================================================================================')
        print('After filtering with PE, remaining:', len(df), '\n')
        info += len(df),

    # tight filter
    elif filter_type == 'tight':
        df = df[df.EXP_HE >= crit.MEAN.HE]
        print('=======================================================================================')
        print('After filtering with HE, remaining:', len(df))
        info += len(df),

        df = df[df.PM <= crit.MEAN.PM]
        print('=======================================================================================')
        print('After filtering with PM, remaining:', len(df))
        info += len(df),

        df = df[df.PE >= crit.MEAN.PE]
        print('=======================================================================================')
        print('After filtering with PE, remaining:', len(df), '\n')
        info += len(df),
    # Default
    else:
        # df = df[df.EXP_HE >= 0.777]
        df = df[df.EXP_HE >= 0.1]
        print('=======================================================================================')
        print('After filtering with HE, remaining:', len(df))
        info += len(df),

        # df = df[df.PM <= 0.087]
        df = df[df.PM <= 0.9]
        print('=======================================================================================')
        print('After filtering with PM, remaining:', len(df))
        info += len(df),

        # df = df[df.PE >= 0.579]
        df = df[df.PE >= 0.1]
        print('=======================================================================================')
        print('After filtering with PE, remaining:', len(df), '\n')
        info += len(df),

def main():
    ################### Phrase Parameters ######################
    import argparse
    Argparser = argparse.ArgumentParser()
    Argparser.add_argument('-f', '--file', action='store', dest='data_in', type=str, help='Input file [.txt]')
    Argparser.add_argument('-v', '--vcf-in', action='store', dest='vcf_in', type=str, default=None)
    Argparser.add_argument('-n', '--n-sample', action='store', dest='n_sample', type=int, default=30)
    Argparser.add_argument('-t', '--thres', action='store', dest='called_threshold', type=float)
    Argparser.add_argument('-p', '--prob-thres', action='store', dest='prob_threshold', type=str, default='0.000001')
    Argparser.add_argument('--tag', action='store', dest='species_name', type=str, default='Unknown')
    Argparser.add_argument('--filter', action='store', dest='filter_type', type=str, default='default')
    Argparser.add_argument('-r', '--rul', action='store', dest='rul', type=str, default='2_3_4')
    Argparser.add_argument('-o', '--out-dir', action='store', dest='output_path', type=str)

    para, unknown = Argparser.parse_known_args()
