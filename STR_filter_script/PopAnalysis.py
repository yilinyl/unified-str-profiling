import itertools as its
import math
import numpy as np
import pandas as pd


def allAlleleInfo(record):
    """
    Summarize information of all alleles at a specific locus
    :param record:
    :return: allele_info_list with each element contains
            ID, length, count, and frequency of the allele
    """
    n_allele = record.num_called * 2
    try:
        alleleList = [record.INFO['REF']] + record.INFO['RPA']
    except KeyError:
        alleleList = [record.INFO['REF']]

    alleleList_dict = {str(allele): 0 for allele in range(len(alleleList))}
    for call in record.samples:
        gt = call['GT']

        # if gt is not None:
        if gt != './.':
            alleles = gt.split('/')
            alleleList_dict[alleles[0]] += 1
            alleleList_dict[alleles[1]] += 1
    # return ALLELE_ID,ALLELE,ALLELE_COUNT,ALLELE_FREQ(allele_count / total number of alleles)
    return [(allele_id, alleleList[allele_id],
             alleleList_dict[str(allele_id)],
             alleleList_dict[str(allele_id)] * 1.0 / n_allele) \
            for allele_id in range(len(alleleList))]


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
                species_dict[name] = 1
            else:
                species_dict[name] += 1

    info = []
    for tag in species_dict.keys():
        info.append(str(tag).strip("'") + ':' + repr(species_dict[tag]))

    return ';'.join(info)


def biasedExpHeterozygosity_freq(alleleFreq):
    he = 1 - math.fsum([pi ** 2 for pi in alleleFreq])
    return he


def unbiasedExpHeterozygosity_record(record):
    """
    Calculate unbiased expectation of heterozygosity with record in VCF file
    :param record:
    :return:
    """

    n_called = record.num_called
    vcf_he = record.heterozygosity  # If there are i alleles with frequency p_i, H=1-sum_i(p_i^2)
    n_allele = 2 * n_called
    unbias_exp_he = vcf_he * n_allele / (n_allele - 1)  # Unbiased Expectation of HE for single locus
    unbias_exp_n_he = unbias_exp_he * n_called  # Unbiased Expectation of HE for n loci

    return unbias_exp_n_he, unbias_exp_he


def unbiasedExpHeterozygosity_freq(alleleFreq, nCalled):
    """
    Calculate unbiased expectation of heterozygosity
    with known allele frequency and the number of calls
    :param alleleFreq:
    :param nCalled:
    :return:
    """
    n_allele = 2 * nCalled
    he = 1 - math.fsum([pi ** 2 for pi in alleleFreq])
    unbias_exp_he = he * n_allele / (n_allele - 1)
    unbias_exp_n_he = unbias_exp_he * nCalled
    return unbias_exp_n_he, unbias_exp_he


def genotypeFreq_record(record, expected=False):
    """

    :param record:
    :param expected:
    :return:
    """
    alleleFreq = [allele[3] for allele in allAlleleInfo(record)]  # allele_count / total number of alleles
    n_allele = record.num_called * 2
    #     gtFreq = []
    if expected:
        gtFreq = [idx[0] == idx[1] and \
                  (alleleFreq[idx[0]] ** 2 - alleleFreq[idx[0]] / n_allele) * n_allele / (n_allele - 1) \
                  or 2 * alleleFreq[idx[0]] * alleleFreq[idx[1]] * n_allele / (n_allele - 1) \
                  for idx in its.combinations_with_replacement(range(len(alleleFreq)), 2)]
    else:
        gtFreq = [idx[0] == idx[1] and alleleFreq[idx[0]] ** 2 or 2 * alleleFreq[idx[0]] * alleleFreq[idx[1]] \
                  for idx in its.combinations_with_replacement(range(len(alleleFreq)), 2)]
    return gtFreq


def genotypeFreq_freq(alleleFreq):
    return [idx[0] == idx[1] and alleleFreq[idx[0]] ** 2 or 2 * alleleFreq[idx[0]] * alleleFreq[idx[1]] \
            for idx in its.combinations_with_replacement(range(len(alleleFreq)), 2)]


def genotypeMaxFreq(record, expected=False):
    """
    Maximum Genotype Frequency
    :param record:
    :param expected:
    :return:
    """
    gtFreq = genotypeFreq_record(record, expected)
    return max(gtFreq)


def pm_record(record, expected=False):
    """
    Matching Probability for single locus
    :param record:
    :param expected:
    :return:
    """
    gtFreq = genotypeFreq_record(record, expected)
    return math.fsum([p ** 2 for p in gtFreq])


def pm_allele_freq(alleleFreq):
    gtFreq = genotypeFreq_freq(alleleFreq)
    return math.fsum([p ** 2 for p in gtFreq])


## half random match probability
def hpm_allele_freq(alleleFreq):
    return math.fsum([p ** 2 * (2 - p) for p in alleleFreq])


def hpm_allele_freq_df(alleleFreq_df):
    return pd.DataFrame({'HPM': {locus: hpm_allele_freq([value for value in row if not np.isnan(value)]) \
                                 for locus, row in alleleFreq_df.iterrows()}})


def pD_record(record, expected=False):
    """
    Power of Discrimination for single locus
    :param record:
    :param expected:
    :return:
    """
    pm = pm_record(record, expected)
    return 1 - pm


def pD_freq(alleleFreq):
    pm = pm_allele_freq(alleleFreq)
    return 1 - pm


def pE1_record(record):
    """
    Power of Exclusion for single locus
    :param record:
    :return:
    """
    alleleFreq = [allele[3] for allele in allAlleleInfo(record)]
    #     pe = 0
    #     for i, pi in enumerate(alleleFreq):
    #         pe += pi * (1 - pi) ** 2
    #         for pj in alleleFreq[(i + 1):]:
    #             pe -= (pi ** 2) * (pj ** 2) * (4 - 3 * pi - 3 * pj)
    return math.fsum([idx[0] == idx[1] and alleleFreq[idx[0]] * (1 - alleleFreq[idx[0]]) ** 2 \
                      or -(alleleFreq[idx[0]] ** 2) * (alleleFreq[idx[1]] ** 2) * \
                      (4 - 3 * alleleFreq[idx[0]] - 3 * alleleFreq[idx[1]]) \
                      for idx in its.combinations_with_replacement(range(len(alleleFreq)), 2)])


# def pE2(record):
#     alleleFreq = [allele[3] for allele in allAlleleInfo(record)]
#     pe = 0
#     for i, pi in enumerate(alleleFreq):
#         pe += pi * (1 - pi + pi ** 2) * (1 - pi) ** 2
#         for pj in alleleFreq[(i + 1):]:
#             pe += pi * pj * (1 - pi + pj ** 2) * (pi + pj)
#     return pe
#
# def pE3(record):
#     he = record.heterozygosity
#     ho = 1 - he
#     return he ** 2 * (1 - 2 * he * ho ** 2)

def pE1_freq(alleleFreq):
    #     pe = 0
    #     for i, pi in enumerate(alleleFreq):
    #         pe = pe + pi * (1 - pi) ** 2
    #         for pj in alleleFreq[(i + 1):]:
    #             pe = pe -(pi ** 2) * (pj ** 2) * (4 - 3 * pi - 3 * pj)
    return math.fsum([idx[0] == idx[1] and alleleFreq[idx[0]] * (1 - alleleFreq[idx[0]]) ** 2 \
                      or -(alleleFreq[idx[0]] ** 2) * (alleleFreq[idx[1]] ** 2) * (
                              4 - 3 * alleleFreq[idx[0]] - 3 * alleleFreq[idx[1]]) \
                      for idx in its.combinations_with_replacement(range(len(alleleFreq)), 2)])


# pd and pe are DataFrame
def cpd(pd):
    return 1 - (1 - pd).prod()


def cpe(pe):
    return 1 - (1 - pe).prod()


def dmp(pm, nEntry):
    cpm = pm.prod()
    if cpm >= 1.0 / nEntry:
        return 1 - (1 - cpm) ** nEntry
    else:
        return cpm * nEntry
