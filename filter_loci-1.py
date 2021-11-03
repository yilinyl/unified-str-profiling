#!/usr/bin/python

import sys
import vcf
import pandas as pd
import Select as st
from decimal import *
import time


# df_origin = pd.read_table(data_in, dtype={'ID': str})
# vcf_canis = vcf.Reader(filename= vcf_in)
# #crit = pd.read_table(crit_file, index_col=0)

def initial_selection(df_origin, rul='2_3_4'):
    """
    Initial selection of STR loci: apply repeat unit length filter
    :param df_origin:
    :param rul:
    :return:
    """
    print('========================================================================================')
    print('Processing begins')
    print('Total loci:', len(df_origin))
    df_origin = df_origin[~df_origin.CHROM.isin(['chrX', 'chrY'])]  # drop chrX and chrY
    print('========================================================================================')
    print('After drop chrX and chrY, remaining:', len(df_origin))
    df_origin = df_origin[df_origin.RUL.isin([int(i) for i in rul.split('_')])]  # select repeat length
    df_origin = df_origin.sort_values(by=['PD', 'GT_FREQ'], ascending=[False, False]).reset_index(drop=True)


def para_selection(called_threshold, n_sample, df_origin):
    """
    Select loci with required HE, PM, PE
    :param called_threshold:
    :param df_origin:
    :param n_sample:
    :return:
    """
    # called_threshold = float(called_threshold)
    # isFound = False
    # filter_info = []
    # final_selected_df = None

    step_filter_info = Decimal(called_threshold).quantize(Decimal('0.01'), rounding=ROUND_UP),  # turple
    # step_filter_info = round(called_threshold, 2)
    print('=======================================================================================')
    print('After selecting repeat length, remaining:', len(df_origin))
    step_filter_info += len(df_origin),  # append filter info

    if called_threshold < 0:
        print('No eligible loci')

    else:
        df = st.filterCalled(df_origin, int(n_sample), called_threshold)  # filter with called number
    print('=======================================================================================')
    print('After filtering nonentity, remaining:', len(df))
    step_filter_info += len(df),

    # select filter type #
    # loose filter

    # Default

    # df = df[df.EXP_HE >= 0.777]
    df = df[df.EXP_HE >= 0.1]
    print('=======================================================================================')
    print('After filtering with HE, remaining:', len(df))
    step_filter_info += len(df),

    # df = df[df.PM <= 0.087]
    df = df[df.PM <= 0.9]
    print('=======================================================================================')
    print('After filtering with PM, remaining:', len(df))
    step_filter_info += len(df),

    # df = df[df.PE >= 0.579]
    df = df[df.PE >= 0.1]
    print('=======================================================================================')
    print('After filtering with PE, remaining:', len(df), '\n')
    step_filter_info += len(df),

    return df


# HE  0.864	0.777
# PD  0.967	0.913
# PE  0.729	0.579
# PM  0.198	0.087

def main():
    
    ################### Phrase Parameters ######################
    import argparse
    Argparser = argparse.ArgumentParser()
    Argparser.add_argument('-f', '--file', action='store', dest='data_in', type=str, help='Input file [.txt]')
    Argparser.add_argument('-v', '--vcf-in', action='store', dest='vcf_in', type=str, default=None)
    Argparser.add_argument('-n', '--n-sample', action='store', dest='n_sample', type=int, default=30)
    Argparser.add_argument('-t', '--thres', action='store', dest='called_threshold', type=float)
    Argparser.add_argument('-p', '--prob-thres', action='store', dest='prob_threshold', type=float, default=0.000001)
    Argparser.add_argument('--tag', action='store', dest='species_name', type=str, default='Unknown')
    Argparser.add_argument('--filter', action='store', dest='filter_type', type=str, default='default')
    Argparser.add_argument('-r', '--rul', action='store', dest='rul', type=str, default='2_3_4')
    Argparser.add_argument('-o', '--out-dir', action='store', dest='output_path', type=str)

    para, unknown = Argparser.parse_known_args()

    ################### Main Body ######################
    print('Selection info.:')
    print('\tNumber of sample:', para.n_sample)
    print('\tCalled threshold:', para.called_threshold)
    print('\tProb threshold:', para.prob_threshold)
    print('\tFilter:', para.filter_type)
    print('\tRUL:', para.rul, '\n\n')
    df_origin = pd.read_table(para.data_in, dtype={'ID': str})
    #vcf_canis = vcf.Reader(filename=para.vcf_in)

    print('========================================================================================')
    print('Processing begins')

    called_threshold = float(para.called_threshold)

    initial_selection(df_origin, rul='2_3_4')
    print('Current called threshold:', called_threshold)
    df = para_selection(called_threshold, para.n_sample, df_origin)
    df.to_csv(para.output_path + para.filter_type + '_' + para.prob_threshold + '_info', sep='\t', index=False)

    return 0


if __name__ == '__main__':
    start_time = time.time()
    return_code = main()
    print("--- %s seconds ---" % (time.time() - start_time))
    sys.exit(return_code)
