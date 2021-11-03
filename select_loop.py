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
    :return: Loci remain after repeat unit length selection
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


def para_filter(df, info):
    """
    Filter of forensic parameters: PM, PE, HE
    :param df:
    :param info:
    :return:
    """

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


def hwe_test(file_in, list_in, hw_df, num):
    """
    Hardy-Weinberg Equilibrium Test
    :param file_in: Input file path (Genepop input format)
    :param list_in: Input loci [list]
    :param hw_df: Input HW dataframe
    :param num: loci number
    :return:
    """
    subprocess.Popen(['Genepop settingsFile=' + file_in + ' Mode=Batch'], shell=True).wait()

    rm_loci = []
    for locus in list_in:
        #            print '===============================================================\ndeficiency hw:'
        #            print locus, hw_deficiency_df[hw_deficiency_df.locus == locus]['P-val']
        if float(hw_df[hw_df.locus == locus]['P-val']) < 0.05 / num:
            rm_loci += [locus]

    return rm_loci


def ld_test(dir_path, df_in, num, rul='2_3_4'):
    """
    Pairwise Linkage Disequilibrium Test
    :param dir_path: 
    :param df_in: 
    :param rul: 
    :param num: 
    :return: 
    """
    process = 'Genepop GameticDiseqTest=Proba settingsFile={dir_path}'
    subprocess.Popen([process.format_map(vars()) + rul + '/Genepopsettings' + rul + '.txt'], shell=True).wait()
    ld_df = st.ld_parser(dir_path + rul + '/canis' + rul + '.txt.DIS')
    remove = []
    for i, row in ld_df.iterrows():
        if row['P-Value'] < 0.05 * 2 / (num * (num - 1)):
            loci = [str(row['Locus#1']), str(row['Locus#2'])]
            cpm_df = df_in[df_in.ID.isin(loci)]
            print(list(df_in.ID))
            print(loci[0], cpm_df[cpm_df.ID == loci[0]].index, cpm_df[cpm_df.ID == loci[0]].PD)
            print(loci[1], cpm_df[cpm_df.ID == loci[1]].index, cpm_df[cpm_df.ID == loci[1]].PD)
            # pairwise check - remove locus with smaller PD value
            if float(cpm_df[cpm_df.ID == loci[0]].PD) > float(cpm_df[cpm_df.ID == loci[1]].PD):
                remove += [loci[1]]
            elif float(cpm_df[cpm_df.ID == loci[0]].PD) == float(cpm_df[cpm_df.ID == loci[1]].PD):
                if cpm_df[cpm_df.ID == loci[0]].index < cpm_df[cpm_df.ID == loci[1]].index:
                    remove += [loci[1]]
                else:
                    remove += [loci[0]]
            else:
                remove += [loci[0]]
    return remove


def main_loop(df_in, vcf_in, dir_path, called_thres, n_sample, prob_thres, rul):
    """
    Loci Selection Loop
    :param df_in: 
    :param vcf_in: 
    :param dir_path: 
    :param called_thres: 
    :param n_sample: 
    :param prob_thres: 
    :param rul: 
    :return: 
    """
    called_thres = float(called_thres)
    isFound = False
    filter_info = []
    final_selected_df = None
    df_sorted = df_in.sort_values(by=['PD', 'GT_FREQ'], ascending=[False, False]).reset_index(drop=True)

    while not isFound:
        step_filter_info = Decimal(called_thres).quantize(Decimal('0.01'), rounding=ROUND_UP),  # turple
        # step_filter_info = round(called_threshold, 2)
        print('=======================================================================================')
        print('After selecting repeat length, remaining:', len(df_sorted))
        step_filter_info += len(df_sorted),  # append filter info
        if called_thres < 0:
            print('No eligible loci')
            return
        df = st.filterCalled(df_sorted, int(n_sample), called_thres)
        para_filter(df, step_filter_info)

        selected = []
        hwe = False
        ld = False

        while not hwe or not ld:

            print('====================================================================')
            print('Loci left:', len(df))
            selected = st.selectLoci(df, float(prob_thres))
            print('====================================================================')
            print('Selected prob:', selected[0])
            if selected[0] > float(prob_thres):
                print('No eligible loci')
                break

            n_loci = len(selected[1])
            df_selected = df.loc[selected[1]]
            # Generate Genepop input file
            st.generateGenepop(df_selected, vcf_in, dir_path + rul + '/canis' + rul + '.txt')

            ctrl = EasyController(dir_path + rul + '/canis' + rul + '.txt')
            IDs = list(df_selected.ID)
            deficiency = []
            excess = []
            for ID in IDs:
                t_fis = ctrl.get_fis(0, ID)
                if t_fis[-1][1] > 0:
                    deficiency += [ID]
                else:
                    excess += [ID]
            df_deficiency = df_selected[df_selected.ID.isin(deficiency)]
            df_excess = df_selected[df_selected.ID.isin(excess)]
            st.generateGenepop(df_deficiency, vcf_in, dir_path + rul + '/canis' + rul + '_deficiency.txt')
            st.generateGenepop(df_excess, vcf_in, dir_path + rul + '/canis' + rul + '_excess.txt')

            remove = []

            if deficiency:
                #        print deficiency
                file = dir_path + rul + '/hwdefict' + rul + '.txt'
                hw_deficiency_df = st.hw_parser(dir_path + rul + '/canis' + rul + '_deficiency.txt.D')
                remove += hwe_test(file, list_in=deficiency, hw_df=hw_deficiency_df, num=n_loci)

            else:
                hw_deficiency_df = None

            if excess:
                file = dir_path + rul + '/hwexcess' + rul + '.txt'
                hw_excess_df = st.hw_parser(dir_path + rul + '/canis' + rul + '_excess.txt.E')
                remove += hwe_test(file, list_in=excess, hw_df=hw_excess_df, num=n_loci)

            else:
                hw_excess_df = None

            if remove:
                hwe = False
                df = df[~df.ID.isin(remove)]
                continue

            else:
                ''' Passed HWE Test '''
                hwe = True
                hw_df = []

                if hw_deficiency_df is not None:
                    hw_df += [hw_deficiency_df]
                if hw_excess_df is not None:
                    hw_df += [hw_excess_df]

                hw_df = pd.concat(hw_df, ignore_index=True)
                colname = list(hw_df.columns)
                colname[0] = 'ID'
                colname[3] = 'Fis'
                hw_df.columns = colname

                remove = ld_test(dir_path, df_selected, n_loci)

                if remove:
                    ld = False
                    df = df[~df.ID.isin(remove)]
                    continue
                else:
                    ld = True
                    isFound = True
                    # print hw_dict, hw_deficiency, hw_excess, deficiency, excess
                    # hw_df = pd.DataFrame({'ID' : [i for i in hw_dict], 'Fis' : [hw_dict[i][2] for i in hw_dict], \
                    # 'p_HWE' : [hw_dict[i][0] for i in hw_dict]})
                    final_selected_df = pd.merge(df_selected, hw_df[['ID', 'P-val', 'Fis']], on=['ID'])
                    print('Final loci left:', len(df))

        filter_info += [step_filter_info]
        if not isFound:
            called_thres -= 0.01  # Reset call threshold
        print('filter info:', filter_info)
        print('Current called threshold:', called_thres)

    return final_selected_df, filter_info, selected


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
    Argparser.add_argument('--type', action='store', dest='type', type=str, default='default')
    Argparser.add_argument('-r', '--rul', action='store', dest='rul', type=str, default='2_3_4')
    Argparser.add_argument('-o', '--out-dir', action='store', dest='out_dir', type=str)

    para, unknown = Argparser.parse_known_args()
    ################### Main Body ######################
    print('Selection info.:')
    print('\tNumber of sample:', para.n_sample)
    print('\tCalled threshold:', para.called_threshold)
    print('\tProb threshold:', para.prob_threshold)
    print('\tType:', para.type)
    print('\tRUL:', para.rul, '\n\n')
    df_origin = pd.read_table(para.data_in, dtype={'ID': str})
    vcf_canis = vcf.Reader(filename=para.vcf_in)
    print('========================================================================================')
    print('Processing begins')

    called_threshold = float(para.called_threshold)

    df_selected = rul_filter(df_origin)
    print('Current called threshold:', called_threshold)
    step_filter_info = Decimal(called_threshold).quantize(Decimal('0.01'), rounding=ROUND_UP),  # turple
    # step_filter_info = round(called_threshold, 2)
    print('=======================================================================================')
    print('After selecting repeat length, remaining:', len(df_origin))
    step_filter_info += len(df_selected),  # append filter info
    final_loci, info, select = main_loop(df_selected, vcf_canis, para.out_dir,
                                         called_threshold, para.n_sample, para.prob_threshold, para.rul)

    info = pd.DataFrame(info, columns=['Called_Threshold', 'RUL_Filtered', 'Called_Filtered', 'HE_Filtered',
                                       'PM_Filtered', 'PE_Filtered', 'Loci_Left', 'Prob'])
    info.to_csv(para.out_dir + '/procedure_' + para.type + '_' + para.prob_threshold + '.info', sep='\t', index=False)
    if final_loci is not None:
        final_loci.to_csv(para.out_dir + '/selected_' + para.type + '_' + para.prob_threshold, sep='\t', index=False)
        with open(para.out_dir + '/selected_' + para.type + '_' + para.prob_threshold, 'a') as fout:
            fout.write('\nProb: ' + str(select[0]) + '\n')
    return 0


if __name__ == '__main__':
    start_time = time.time()
    return_code = main()
    print("--- %s seconds ---" % (time.time() - start_time))
    sys.exit(return_code)
