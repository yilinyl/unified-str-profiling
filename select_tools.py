#!/usr/bin/env python3

import sys
import os
import re
import pandas as pd
import vcf
from io import StringIO


def filter_called(df, n_sample, n_called_thres):
    return df[df.N_CALLED >= n_sample * n_called_thres]


def select_prob(sorted_df, prob_thres):
    """
    Select Loci
    :param sorted_df:
    :param prob_thres:
    :return: maximum genotype frequency of selected loci, selected loci set
    """
    prob = 1
    selected = []
    for i, locus in sorted_df.iterrows():
        if prob <= prob_thres:
            break
        selected += [i]
        prob *= locus.GT_FREQ  # maximum genotype frequency of selected loci
    return prob, selected


def generate_genepop_file(df_selected, vcf_canis, out_dir='./'):
    """
    Generate Genepop input file
    :param df_selected:
    :param vcf_canis:
    :param out_dir:
    :return:
    """
    #    chroms = list(df_selected.CHROM.drop_duplicates())
    #    df_chrom = df_selected.sort()
    df_selected_sorted = df_selected.sort_index()
    species_dict = {}
    for i, locus in df_selected_sorted.iterrows():
        species_in = locus['SPECIES'].split(';')
        locus_species = []
        for j in range(len(species_in)):
            taxon = species_in[j].split(':')[0]
            if taxon not in species_dict:
                species_dict[taxon] = {'output': out_dir + taxon + '_genepop.txt', 'locus_id': [],
                                       'samples': [], 'locus_gt': {}, 'loci_gt': []}
            species_dict[taxon]['locus_id'] += [locus.ID]
            locus_species.append(taxon)

        # record = vcf_canis.fetch(locus.CHROM, locus.POS)
        for x in vcf_canis.fetch(locus.CHROM, locus.POS - 1, locus.POS):
            record = x
        # print locus.CHROM, locus.POS
        for sample in record.samples:
            taxon = sample.sample.split('_')[0]
            genotype = sample['GT']
            if taxon in locus_species:
                if sample.sample not in species_dict[taxon]['samples']:  # Generate sample list for each species
                    species_dict[taxon]['samples'].append(sample.sample)
                # if genotype is None:
                if genotype == './.':  # missing genotype
                    species_dict[taxon]['locus_gt'].update({sample.sample: '0000'})
                else:
                    coded_gt = ''
                    genotype = genotype.split('/')
                    for allele in genotype:
                        coded_gt += str(int(allele) + 1).zfill(2)
                    species_dict[taxon]['locus_gt'].update({sample.sample: coded_gt})
        for tag in species_dict:
            if tag in locus_species:
                species_dict[tag]['loci_gt'].append(species_dict[tag]['locus_gt'])
                species_dict[tag]['locus_gt'] = {}

    for tag in species_dict:
        output = species_dict[tag]['output']
        if not os.path.exists(output):
            my_dir, file = os.path.split(output)
            if not os.path.exists(my_dir):
                os.makedirs(my_dir)
            os.mknod(output)
        # Generate Genepop input file
        with open(output, 'w') as fout:
            fout.write('Genepop canis selected\n' +
                       '\n'.join([str(i) for i in species_dict[tag]['locus_id']]) + '\n' + 'Pop\n')
            # tag = 'cat'
            for sample in species_dict[tag]['samples']:
                indv = []
                # if tag != sample.split('_')[0]:
                #     fout.write('Pop\n')
                #     tag = sample.split('_')[0]
                for i in range(len(species_dict[tag]['loci_gt'])):
                    indv += [species_dict[tag]['loci_gt'][i][sample]]
                fout.write(sample + ', ' + ' '.join(indv) + '\n')
    all_output = [species_dict[tag]['output'] for tag in species_dict]
    return all_output


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
                key, value = map(lambda x: x.strip(), cnt[idx].replace('\n', '').split(':'))
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
    table = re.sub(r'matrices', '', table)
    table = re.sub(r' +\n', '\n', table)
    table = re.sub(r' +', ',', table)
    df = pd.read_csv(StringIO(table), sep=',', index_col=False, parse_dates=False, dtype={'locus': str})
    return df


if __name__ == '__main__':
    file_in = '/home/lyl/Documents/STR_select/test_merged_info'
    vcf_file = '/home/lyl/Documents/STR_select/all.merged.vcf.gz'
    df = pd.read_table(file_in, dtype={'ID': str})
    vcf_in = vcf.Reader(filename=vcf_file)
    file_list = generate_genepop_file(df, vcf_in)
    print(file_list)
