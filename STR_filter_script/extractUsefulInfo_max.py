#!/usr/bin/env python3

import sys
import vcf
import PopAnalysis as pa


# file_in = sys.argv[1]
# file_out = sys.argv[2]
# species_name = sys.argv[3]
def extract_info(file_in, file_out):
    """
    Exract useful information from merged VCF file
    :param file_in:
    :param file_out:
    :return:
    """
    vcf_canis = vcf.Reader(filename=file_in)
    with open(file_out, 'w') as fout:
        fout.write('\t'.join(['ID', 'CHROM', 'POS', 'LOCUS_ID', 'RUL', 'N_CALLED',
                              'N_ALLELE', 'N_OBS_HO', 'N_EXP_HO', 'EXP_HO', 'N_OBS_HE',
                              'N_EXP_HE', 'EXP_HE', 'PM', 'PD', 'PE', 'GT_FREQ', 'SPECIES',
                              'ALLELE_ID:ALLELE:ALLELE_COUNT:ALLELE_FREQ']) + '\n')
        index = 1
        cChrom = 'chr1'
        for record in vcf_canis:
            if record.CHROM.startswith('chrUn_'):
                continue
            if record.CHROM != cChrom:
                index = 1
                cChrom = record.CHROM
            print('{}\t\t{}'.format(record.CHROM, record.POS))
            locus_uid = cChrom.split('r')[1] + '.' + repr(index)
            locusInfoTuple = locus_uid, cChrom, record.POS, record.ID, len(record.INFO['RU']), \
                             record.num_called, record.ALT is None and 1 or len(record.ALT) + 1, \
                             record.num_called - record.num_het, \
                             record.num_called - pa.unbiasedExpHeterozygosity_record(record)[0], \
                             1 - pa.unbiasedExpHeterozygosity_record(record)[1], record.num_het, \
                             pa.unbiasedExpHeterozygosity_record(record)[0], \
                             pa.unbiasedExpHeterozygosity_record(record)[1], pa.pm_record(record), \
                             pa.pD_record(record), pa.pE1_record(record), \
                             pa.genotypeMaxFreq(record), pa.species_info(record), pa.allAlleleInfo(record)
            line_out = [isinstance(i, str) and i or repr(i) for i in locusInfoTuple[:-1]] + \
                       [';'.join([':'.join([repr(i) for i in allele]) for allele in locusInfoTuple[-1]])]
            fout.write('\t'.join(line_out) + '\n')
            index += 1


if __name__ == '__main__':
    vcf_in = sys.argv[1]
    output = sys.argv[2]
    extract_info(vcf_in, output)
