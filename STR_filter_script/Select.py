import sys
import os

if sys.version_info[0] < 3:
    from StringIO import StringIO
else:
    from io import StringIO

import re
import string

import pandas as pd


def filterCalled(df, nSample, nCalledProportionThreshold):
    return df[df.N_CALLED >= nSample * nCalledProportionThreshold]


def selectLoci(sortedDf, probThreshold):
    """
    Select Loci
    :param sortedDf:
    :param probThreshold:
    :return: maximum genotype frequency of selected loci, selected loci set
    """
    prob = 1
    selected = []
    for i, locus in sortedDf.iterrows():
        if prob <= probThreshold:
            break
        selected += [i]
        prob *= locus.GT_FREQ  # maximum genotype frequency of selected loci
    return prob, selected


def generateGenepop(df_selected, vcf_canis, outputfile):
    """
    Generate Genepop input file
    :param df_selected:
    :param vcf_canis:
    :param outputfile:
    :return:
    """
    #    chroms = list(df_selected.CHROM.drop_duplicates())
    #    df_chrom = df_selected.sort()
    df_selected_sorted = df_selected.sort_index()
    loci_id = []
    loci_gt = []
    for i, locus in df_selected_sorted.iterrows():
        loci_id += [locus.ID]
        locus_gt = {}
        # record = vcf_canis.fetch(locus.CHROM, locus.POS)
        for x in vcf_canis.fetch(locus.CHROM, locus.POS - 1, locus.POS):
            record = x
        # print locus.CHROM, locus.POS
        for sample in record.samples:
            genotype = sample['GT']
            # if genotype is None:
            if genotype == './.':                         # missing genotype
                locus_gt.update({sample.sample: '0000'})
            else:
                coded_gt = ''
                genotype = genotype.split('/')
                for allele in genotype:
                    coded_gt += str(int(allele) + 1).zfill(2)
                locus_gt.update({sample.sample: coded_gt})
        loci_gt.append(locus_gt)

    if not os.path.exists(outputfile):
        my_dir, file= os.path.split(outputfile)
        if not os.path.exists(my_dir):
            os.makedirs(my_dir)
        os.mknod(outputfile)

    # Generate Genepop input file
    with open(outputfile, 'w') as fout:
        fout.write('Genepop canis selected\n' + '\n'.join([str(i) for i in loci_id]) + '\n' + 'Pop\n')
        for sample in vcf_canis.samples:
            indv = []
            for i in range(len(loci_gt)):
                indv += [loci_gt[i][sample]]
            fout.write(sample + ', ' + ' '.join(indv) + '\n')


def ld_parser(filename):
    with open(filename) as ifs:
        cnt = ifs.readlines()

    table = ""
    popu_no = 0
    loci_no = 0

    hhm_param = {}

    idx = 0
    while idx < len(cnt):
        line = cnt[idx]
        if 'Number of populations detected' in line:
            #            print 'parse population number'
            popu_no = int(line.replace('\n', '').split(':')[1])
        if 'Number of loci detected' in line:
            #            print 'parse loci number'
            loci_no = int(line.replace('\n', '').split(':')[1])
        if 'Markov chain parameters' in line:
            #            print 'parse parameters'
            idx += 1
            while idx < len(cnt) and cnt[idx] != '\n':
                key, value = map(string.strip, cnt[idx].replace('\n', '').split(':'))
                hhm_param[key] = float(value)
                idx += 1
        if re.match(r'--*', line):
            #            print 'parse table'
            table += cnt[idx - 1]  # header
            idx += 1
            while idx < len(cnt) and cnt[idx] != '\n':
                table += cnt[idx]
                idx += 1
        idx += 1

    table = re.sub(r'  *', ',', table)
    #    print table
    df = pd.read_csv(StringIO(table), sep=',', index_col=False, parse_dates=False,
                     dtype={'Locus#1': str, 'Locus#2': str})
    #    print df

    return df


def hw_parser(filename):
    with open(filename) as fin:
        cnt = fin.readlines()
    table = ""
    idx = 0
    while idx < len(cnt):
        line = cnt[idx]
        if re.match(r'\s+--*\n$', line):
            table += cnt[idx + 1]  # header
            idx += 3
            while idx < len(cnt) and cnt[idx] != '\n':
                table += cnt[idx]
                idx += 1
        idx += 1
    table = re.sub(r' switches.*', '', table)
    table = re.sub(r' +\n', '\n', table)
    table = re.sub(r' +', ',', table)
    df = pd.read_csv(StringIO(table), sep=',', index_col=False, parse_dates=False, dtype={'locus': str})
    return df
