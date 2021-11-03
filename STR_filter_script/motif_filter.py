#!/usr/bin/python

import sys, vcf

file_in = sys.argv[1]
file_out = sys.argv[2]

vcf_input = vcf.Reader(open(file_in, 'r'))
with open(file_out, 'w') as fout:
    fout.write('\t'.join(['ID', 'CHROM', 'POS', 'LOCUSID', 'RUL', 'N_CALLED','N_ALLELE']))

