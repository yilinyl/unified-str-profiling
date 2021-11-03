#!/usr/bin/env python3

import sys
import vcf


def get_sample(file_in):
    vcf_in = vcf.Reader(filename=file_in)
    for record in vcf_in:
        #        for call in record:
        #            print(call.sample.split('_')[0])
        print(record.INFO)


def species_info(record):
    """
    Get species information at each locus
    :param record:
    :return:
    """
    species_dict = {}
    for call in record.samples:
        if call.called:
            name = call.sample.split('_')[0]
            if name not in species_dict.keys():
                species_dict[name] = 0
            else:
                species_dict[name] += 1
    return [(tag, species_dict[tag]) for tag in species_dict.keys()]


if __name__ == '__main__':
    # path = sys.argv[1]
    vcf_input = vcf.Reader(filename='/home/lyl/Documents/STR_select/cattle.merged.vcf.gz')
    file_out = '/home/lyl/Documents/STR_select/cattle.merged.txt'
    with open(file_out, 'w') as fout:
        for read in vcf_input:
            species_turple = species_info(read)
            fout.write(str(species_turple))
