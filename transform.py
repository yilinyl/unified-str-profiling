#!/usr/bin/env python3

import sys
import vcf
import PopAnalysis as pa


def extract_info(file_in, file_out):
    """
    Extract useful information from merged VCF file
    :param file_in: VCF file
    :param file_out: TXT file
    :return:
    """
    vcf_canis = vcf.Reader(filename=file_in)
    with open(file_out, 'w') as fout:
        fout.write('\t'.join(['ID', 'CHROM', 'POS', 'LOCUS_ID', 'RUL', 'N_CALLED',
                              'N_ALLELE', 'N_OBS_HO', 'N_EXP_HO', 'EXP_HO', 'N_OBS_HE',
                              'N_EXP_HE', 'EXP_HE', 'PM', 'PD', 'PE', 'GT_FREQ', 'SPECIES',
                              'ALLELE_ID:ALLELE:ALLELE_COUNT:ALLELE_FREQ']) + '\n')
        index = 1
        tmp_chr = 'chr1'
        for record in vcf_canis:
            if record.CHROM.startswith('chrUn_'):
                continue
            if record.CHROM != tmp_chr:
                index = 1
                tmp_chr = record.CHROM
            print('{}\t\t{}'.format(record.CHROM, record.POS))
            locus_uid = tmp_chr.split('r')[1] + '.' + repr(index)
            locus_info = locus_uid, tmp_chr, record.POS, record.ID, len(record.INFO['RU']), \
                         record.num_called, record.ALT is None and 1 or len(record.ALT) + 1, \
                         record.num_called - record.num_het, \
                         record.num_called - pa.unbiased_exp_hetero_record(record)[0], \
                         1 - pa.unbiased_exp_hetero_record(record)[1], record.num_het, \
                         pa.unbiased_exp_hetero_record(record)[0], \
                         pa.unbiased_exp_hetero_record(record)[1], pa.pm_record(record), \
                         pa.pd_record(record), pa.pe1_record(record), \
                         pa.genotype_max_freq(record), pa.species_info(record), pa.all_allele_info(record)
            line_out = [isinstance(i, str) and i or repr(i) for i in locus_info[:-1]] + \
                       [';'.join([':'.join([repr(i) for i in allele]) for allele in locus_info[-1]])]
            fout.write('\t'.join(line_out) + '\n')
            index += 1


if __name__ == '__main__':
    vcf_in = '/home/lyl/Documents/STR_select/cattle.merged.vcf.gz'
    output = '/home/lyl/Documents/STR_select/cattle.merged.txt'
    extract_info(vcf_in, output)
