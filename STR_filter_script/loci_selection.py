#!/usr/bin/python

import sys, vcf, subprocess
# sys.path.insert(0, '/home/liangxue/Desktop/Horse/process/filter')
#from _pydecimal import Decimal
#from typing import Any, Tuple, Union

import pandas as pd
import Select as st
from Bio.PopGen.GenePop.EasyController import EasyController
from decimal import *
import time

start_time = time.time()

data_in = sys.argv[1]  # txt input
# file_out = sys.argv[2]
vcf_in = sys.argv[2]
n_sample = sys.argv[3]  # total number of training samples
called_threshold = sys.argv[4]  # minimum frequency of certain locus
prob_threshold = sys.argv[5]    # 1/S, where S stands for the largest posible population size
crit_file = sys.argv[6]
filter_type = sys.argv[7]  # [loose, tight, default], if loose or tight filter is chosen, crit_file is required
rul = sys.argv[8]          # repeat unit lengths, seperate the numbers with'_' E.g. 3_4_5
output_path = sys.argv[9]
dir_path = sys.argv[10]

print('Selection info.:')
print('\tNumber of sample:', n_sample)
print('\tCalled threshold:', called_threshold)
print('\tProb threshold:', prob_threshold)
print('\tFilter:', filter_type)
print('\tRUL:', rul, '\n\n')

'''
   01.Initialization - 
   exclude loci from chrX, chrY. 
   select loci with required repeat unit length(RUL),
   sort remaining loci (non-ascending) by PD(power of discrimination), genotype frequency
'''
df_origin = pd.read_table(data_in, dtype={'ID': str})
vcf_canis = vcf.Reader(filename= vcf_in)
crit = pd.read_table(crit_file, index_col=0)

print('========================================================================================')
print('Processing begins')
print('Total loci:', len(df_origin))
df_origin = df_origin[~df_origin.CHROM.isin(['chrX', 'chrY'])]  # drop chrX and chrY
print('========================================================================================')
print('After drop chrX and chrY, remaining:', len(df_origin))

# if rul in ['1', '2', '3', '4']:
#     df_origin = df_origin[df_origin.RUL == int(rul)] # select repeat length
df_origin = df_origin[df_origin.RUL.isin([int(i) for i in rul.split('_')])]  # select repeat length

df_origin = df_origin.sort_values(by=['PD', 'GT_FREQ'], ascending=[False, False])
df_origin = df_origin.reset_index()
df_origin = df_origin.drop(df_origin.columns[0], axis=1)
#############################

def loci_filter(df, info):
    """

    :param df:
    :param info:
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

''' Apply Filter '''
called_threshold = float(called_threshold)
isFound = False
filter_info = []
final_selected_df = None
while not isFound:
    step_filter_info = Decimal(called_threshold).quantize(Decimal('0.01'), rounding=ROUND_UP),    # turple
    #step_filter_info = round(called_threshold, 2)
    print('=======================================================================================')
    print('After selecting repeat length, remaining:', len(df_origin))
    step_filter_info += len(df_origin),    # append filter info

    if called_threshold < 0:
        print('No eligible loci')
        break
    else:
        df = st.filterCalled(df_origin, int(n_sample), called_threshold)  # filter with called number
    print('=======================================================================================')
    print('After filtering nonentity, remaining:', len(df))
    step_filter_info += len(df),



        # HE  0.864	0.777
        # PD  0.967	0.913
        # PE  0.729	0.579
        # PM  0.198	0.087

    #########################################


    selected = []
    hwe = False
    ld = False


    ''' Hardy-Weinberg Equilibrium & Linkage Disequilibrium Testing '''

    while not hwe or not ld:

        print('====================================================================')
        print('Loci left:', len(df))
        selected = st.selectLoci(df, float(prob_threshold))
        print('====================================================================')
        print('Selected prob:', selected[0])
        if selected[0] > float(prob_threshold):
            print('No eligible loci')
            break

        n_loci = len(selected[1])
        df_selected = df.loc[selected[1]]
        # Generate Genepop input file
        st.generateGenepop(df_selected, vcf_canis, dir_path + rul + '/canis' + rul + '.txt')

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
        st.generateGenepop(df_deficiency, vcf_canis, dir_path + rul + '/canis' + rul + '_deficiency.txt')
        st.generateGenepop(df_excess, vcf_canis, dir_path + rul + '/canis' + rul + '_excess.txt')

        remove = []
        if deficiency:
            #        print deficiency
            subprocess.Popen(['Genepop settingsFile=' + dir_path + rul + '/hwdefict' + rul + '.txt Mode=Batch'],
                             shell=True).wait()
            hw_deficiency_df = st.hw_parser(dir_path + rul + '/canis' + rul + '_deficiency.txt.D')
            for locus in deficiency:
                #            print '===============================================================\ndeficiency hw:'
                #            print locus, hw_deficiency_df[hw_deficiency_df.locus == locus]['P-val']
                if float(hw_deficiency_df[hw_deficiency_df.locus == locus]['P-val']) < 0.05 / n_loci:
                    remove += [locus]
        else:
            hw_deficiency_df = None

        if excess:
            subprocess.Popen(['Genepop settingsFile=' + dir_path + rul + '/hwexcess' + rul + '.txt Mode=Batch'],
                             shell=True).wait()
            hw_excess_df = st.hw_parser(dir_path + rul + '/canis' + rul + '_excess.txt.E')
            for locus in excess:
                #            print '==============================================================\nexcess hw:'
                #            print float(hw_excess_df[hw_excess_df.locus == locus]['P-val'])
                if float(hw_excess_df[hw_excess_df.locus == locus]['P-val']) < 0.05 / n_loci:
                    remove += [locus]
        else:
            hw_excess_df = None

        #         print IDs
        #         raw_input('Continue')
        #    print 'remove:'
        #    print remove
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


            process='Genepop GameticDiseqTest=Proba settingsFile={dir_path}'
            subprocess.Popen([process.format_map(vars()) + rul + '/Genepopsettings' + rul + '.txt'], shell=True).wait()
            ld_df = st.ld_parser(dir_path + rul + '/canis' + rul + '.txt.DIS')
            for i, row in ld_df.iterrows():
                if row['P-Value'] < 0.05 * 2 / (n_loci * (n_loci - 1)):
                    loci = [str(row['Locus#1']), str(row['Locus#2'])]
                    cpm_df = df_selected[df_selected.ID.isin(loci)]
                    print(list(df_selected.ID))
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
            #                     print remove[-1]
            #                     raw_input('Continue')
            #                     weight = 0
            #                     if int(cpm_df[cpm_df.ID == loci[0]].N_CALLED) >= int(cpm_df[cpm_df.ID == loci[1]].N_CALLED):
            #                         weight += 1
            #                     if float(cpm_df[cpm_df.ID == loci[0]].EXP_HE) >= float(cpm_df[cpm_df.ID == loci[1]].EXP_HE):
            #                         weight += 1
            #                     if float(cpm_df[cpm_df.ID == loci[0]].PM) <= float(cpm_df[cpm_df.ID == loci[1]].PM):
            #                         weight += 1
            #                     if float(cpm_df[cpm_df.ID == loci[0]].PE) >= float(cpm_df[cpm_df.ID == loci[1]].PE):
            #                         weight += 1
            #                     # if float(hw_df[hw_df.ID == loci[0]]['P-val']) > float(hw_df[hw_df.ID == loci[1]]['P-val']):
            #                     #     weight += 1
            #                     if weight == 4:
            #                         remove += [loci[1]]
            #                     else:
            #                         remove += [loci[0]]
            if remove:
                ld = False
                df = df[~df.ID.isin(remove)]
                continue
            else:
                ld = True
                isFound = True
                #print hw_dict, hw_deficiency, hw_excess, deficiency, excess
                #hw_df = pd.DataFrame({'ID' : [i for i in hw_dict], 'Fis' : [hw_dict[i][2] for i in hw_dict], 'p_HWE' : [hw_dict[i][0] for i in hw_dict]})
                final_selected_df = pd.merge(df_selected, hw_df[['ID', 'P-val', 'Fis']], on=['ID'])
                print('Final loci left:', len(df))

    #step_filter_info += len(df), selected[0]
    filter_info += [step_filter_info]
    if not isFound:
        called_threshold -= 0.01
    print('filter info:', filter_info)
    print('Current called threshold:', called_threshold)

#df.to_csv(output_path + filter_type + '_' + prob_threshold + 'info', sep='\t', index=False)

print(filter_info)
filter_info = pd.DataFrame(filter_info, columns=['Called_Threshold', 'RUL_Filtered',
                                                 'Called_Filtered', 'HE_Filtered', 'PM_Filtered',
                                                 'PE_Filtered', 'Loci_Left', 'Prob'])

filter_info.to_csv(output_path + '/procedure_' + filter_type + '_' + prob_threshold + '.info', sep='\t', index=False)
if final_selected_df is not None:
    final_selected_df.to_csv(output_path + '/selected_' + filter_type + '_' + prob_threshold, sep='\t', index=False)
    with open(output_path + '/selected_' + filter_type + '_' + prob_threshold, 'a') as fout:
        fout.write('\nProb: ' + str(selected[0]) + '\n')
print("--- %s seconds ---" % (time.time() - start_time))

