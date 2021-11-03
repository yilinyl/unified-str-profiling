#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 18-10-23

@author: lyl
"""
import pandas as pd
import select_tools
import filter_loci
import vcf
from Bio.PopGen.GenePop.EasyController import EasyController


def select_para(file_in, out_dir='./', called_thres=0.1, num=300):
    df_origin = pd.read_table(file_in, dtype={'ID': str})
    out_file = out_dir + 'all_merged_info'
    df_select = filter_loci.initial_selection(df_origin)
    df_select = filter_loci.para_selection(called_thres, num, df_select)
    df_select.to_csv(out_file, sep='\t', index=False)
    return out_file


def pre_genepop(info_file, vcf_file):
    df_raw = pd.read_table(info_file, dtype={'ID': str})
    vcf_in = vcf.Reader(filename=vcf_file)
    all_file = select_tools.generate_genepop_file(df_raw, vcf_in)
    return all_file


def hw_test(all_file):
    for file_raw in all_file:
        ctrl = EasyController(file_raw)
        loci_map = ctrl.test_hw_global()
    

