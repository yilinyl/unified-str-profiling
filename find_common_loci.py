#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 15:58:23 2018

@author: lyl
"""

import glob
import os
import matplotlib.pyplot as plt
import pandas as pd


def load_data(path):
    """
    读取STR位点统计数据
    :param path:
    :return df_select:
    """

    flist = glob.glob(path)
    df = pd.DataFrame()
    for file_path in flist:
        df1 = pd.read_csv(file_path, sep='\t')
        tag = os.path.split(file_path)[1].split('_')[0]
        df1['TAG'] = tag
        df = df.append(df1)
    df_select = df.sort_values(by=['CHROM', 'POS', 'ID']).reset_index(drop=True)
    df_select = df_select.drop(labels='ID', axis=1)
    tmp = df_select.pop('TAG')
    df_select.insert(0, 'TAG', tmp)
    set_locus_id(df_select)
    return df_select


def set_locus_id(df_select):
    """
    设置位点ID(根据位点在染色体上的坐标)
    :param df_select:
    :return:
    """

    idx = 0
    tmp_chr = 'chr1'
    pos = '0'
    for i, locus in df_select.iterrows():
        if locus.CHROM == tmp_chr:
            if locus.POS != pos:
                pos = locus.POS
                idx += 1
        else:
            tmp_chr = locus.CHROM
            idx = 1
        df_select.loc[i, 'LOCUSID'] = tmp_chr.split('r')[1] + '.' + repr(idx)


def plot_distribution(df_select):
    """
    绘制STR位点分布散点图
    :param df_select:
    :return:
    """

    counts = df_select.CHROM.value_counts()
    # tmp = 0
    # set_locus_id(df_select)
    thres = 10
    plt.figure(figsize=(30, 20))
    for i in range(len(counts)):
        if counts.iloc[i] < thres:
            break
        tmp_chr = counts.keys()[i]
        # num = counts.iloc[i]
        ax = plt.subplot(3, 5, i + 1)
        ax.grid(axis='y', linestyle=':')
        ax.scatter(df_select.loc[df_select.CHROM == tmp_chr]['TAG'],
                   df_select.loc[df_select.CHROM == tmp_chr]['LOCUSID'], marker='D')
        ax.set_title(tmp_chr)
    plt.show()


def get_common_loci(df_select, outdir):
    """
    获取物种共有位点
    :param df_select:
    :param outpath:
    :return:
    """

    sum_ID = df_select.LOCUSID.value_counts(sort=False).sort_index().to_frame().reset_index()
    sum_ID = sum_ID.rename(index=str, columns={'index': 'LOCUS_ID', 'LOCUSID': 'COUNT'})
    print("多物种共有位点：{num}".format(num=repr(len(sum_ID.loc[sum_ID.COUNT > 1]))))
    sum_common_loci = sum_ID.loc[sum_ID.COUNT > 1]
    sum_common_loci.to_csv(outdir + 'Summary.csv', sep='\t', index=False)
    #    common_loci = list(sum_ID.loc[sum_ID > 1].keys())
    df_common_loci = df_select[df_select.LOCUSID.isin(sum_common_loci.LOCUS_ID)]
    #    df_common_loci.to_csv('/home/lyl/Documents/STR_select/stats/CommonLoci.csv')
    df_common_loci.to_csv(outdir + 'CommonLoci.csv', sep='\t', index=False)


if __name__ == '__main__':
    df_select = load_data('/home/lyl/Documents/STR_select/out_info/all_default*')
    get_common_loci(df_select, '/home/lyl/Documents/STR_select/stats/')
#    plot_distribution(df_select)
