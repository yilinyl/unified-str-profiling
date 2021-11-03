import itertools as its
import math
import numpy as np
import pandas as pd


def all_allele_info(record):
    """
    Summarize information of all alleles at a specific locus
    :param record:
    :return: allele_info_list with each element contains
            ID, length, count, and frequency of the allele
    """
    n_allele = record.num_called * 2
    try:
        allele_list = [record.INFO['REF']] + record.INFO['RPA']
    except KeyError:
        allele_list = [record.INFO['REF']]

    allele_dict = {str(allele): 0 for allele in range(len(allele_list))}
    for call in record.samples:
        gt = call['GT']

        # if gt is not None:
        if gt != './.':
            alleles = gt.split('/')
            allele_dict[alleles[0]] += 1
            allele_dict[alleles[1]] += 1
    # return ALLELE_ID,ALLELE,ALLELE_COUNT,ALLELE_FREQ(allele_count / total number of alleles)
    return [(allele_id, allele_list[allele_id],
             allele_dict[str(allele_id)],
             allele_dict[str(allele_id)] * 1.0 / n_allele) \
            for allele_id in range(len(allele_list))]


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


def biased_exp_hetero_freq(allele_freq):
    he = 1 - math.fsum([pi ** 2 for pi in allele_freq])
    return he


def unbiased_exp_hetero_record(record):
    """
    Calculate unbiased expectation of heterozygosity with record in VCF file
    :param record:
    :return:
    """

    n_called = record.num_called
    vcf_he = record.heterozygosity  # If there are i alleles with frequency p_i, H=1-sum_i(p_i^2)
    n_allele = 2 * n_called
    unbiased_exp_he = vcf_he * n_allele / (n_allele - 1)  # Unbiased Expectation of HE for single locus
    unbiased_exp_n_he = unbiased_exp_he * n_called  # Unbiased Expectation of HE for n loci

    return unbiased_exp_n_he, unbiased_exp_he


def unbiased_exp_hetero_freq(allele_freq, n_called):
    """
    Calculate unbiased expectation of heterozygosity
    with known allele frequency and the number of calls
    :param allele_freq:
    :param n_called:
    :return:
    """
    n_allele = 2 * n_called
    he = 1 - math.fsum([pi ** 2 for pi in allele_freq])
    unbiased_exp_he = he * n_allele / (n_allele - 1)
    unbiased_exp_n_he = unbiased_exp_he * n_called
    return unbiased_exp_n_he, unbiased_exp_he


def genotype_freq_record(record, expected=False):
    """

    :param record:
    :param expected:
    :return:
    """
    allele_freq = [allele[3] for allele in all_allele_info(record)]  # allele frequency list
    n_allele = record.num_called * 2
    gt_freq = []
    if expected:
        gt_freq = [idx[0] == idx[1] and
                  (allele_freq[idx[0]] ** 2 - allele_freq[idx[0]] / n_allele) * n_allele / (n_allele - 1) \
                  or 2 * allele_freq[idx[0]] * allele_freq[idx[1]] * n_allele / (n_allele - 1) \
                  for idx in its.combinations_with_replacement(range(len(allele_freq)), 2)]
    else:
        # gt_freq = [idx[0] == idx[1] and allele_freq[idx[0]] ** 2 or 2 * allele_freq[idx[0]] * allele_freq[idx[1]] \
        #           for idx in its.combinations_with_replacement(range(len(allele_freq)), 2)]
        for idx in its.combinations_with_replacement(range(len(allele_freq)), 2):
            if idx[0] == idx[1]:
                gt_freq.append(allele_freq[idx[0]] ** 2)
            else:
                gt_freq.append(2 * allele_freq[idx[0]] * allele_freq[idx[1]])

    return gt_freq


def genotype_freq(allele_freq):
    gt_freq = []
    for idx in its.combinations_with_replacement(range(len(allele_freq)), 2):
        if idx[0] == idx[1]:
            gt_freq.append(allele_freq[idx[0]] ** 2)
        else:
            gt_freq.append(2 * allele_freq[idx[0]] * allele_freq[idx[1]])

    return gt_freq


def genotype_max_freq(record, expected=False):
    """
    Maximum Genotype Frequency
    :param record:
    :param expected:
    :return:
    """
    gt_freq = genotype_freq_record(record, expected)
    return max(gt_freq)


def pm_record(record, expected=False):
    """
    Matching Probability for single locus
    :param record:
    :param expected:
    :return:
    """
    gt_freq = genotype_freq_record(record, expected)
    return math.fsum([p ** 2 for p in gt_freq])


def pm_allele_freq(allele_freq):
    gt_freq = genotype_freq(allele_freq)
    return math.fsum([p ** 2 for p in gt_freq])


# half random match probability
def hpm_allele_freq(allele_freq):
    return math.fsum([p ** 2 * (2 - p) for p in allele_freq])


def hpm_allele_freq_df(allele_freq_df):
    return pd.DataFrame({'HPM': {locus: hpm_allele_freq([value for value in row if not np.isnan(value)]) \
                                 for locus, row in allele_freq_df.iterrows()}})


def pd_record(record, expected=False):
    """
    Power of Discrimination for single locus
    :param record:
    :param expected:
    :return:
    """
    pm = pm_record(record, expected)
    return 1 - pm


def pd_freq(allele_freq):
    pm = pm_allele_freq(allele_freq)
    return 1 - pm


def pe1_record(record):
    """
    Power of Exclusion for single locus
    :param record:
    :return:
    """
    allele_freq = [allele[3] for allele in all_allele_info(record)]
    #     pe = 0
    #     for i, pi in enumerate(allele_freq):
    #         pe += pi * (1 - pi) ** 2
    #         for pj in allele_freq[(i + 1):]:
    #             pe -= (pi ** 2) * (pj ** 2) * (4 - 3 * pi - 3 * pj)
    return math.fsum([idx[0] == idx[1] and allele_freq[idx[0]] * (1 - allele_freq[idx[0]]) ** 2 \
                      or -(allele_freq[idx[0]] ** 2) * (allele_freq[idx[1]] ** 2) * \
                      (4 - 3 * allele_freq[idx[0]] - 3 * allele_freq[idx[1]]) \
                      for idx in its.combinations_with_replacement(range(len(allele_freq)), 2)])



def pe1_freq(allele_freq):
    pe_tmp = []
    for idx in its.combinations_with_replacement(range(len(allele_freq)), 2):
        if idx[0] == idx[1]:
            pe_tmp.append(allele_freq[idx[0]] * (1 - allele_freq[idx[0]]) ** 2)
        else:
            pe_tmp.append(-(allele_freq[idx[0]] ** 2) * (allele_freq[idx[1]] ** 2) *
                           (4 - 3 * allele_freq[idx[0]] - 3 * allele_freq[idx[1]]))
    return math.fsum(pe_tmp)
    #     pe = 0
    #     for i, pi in enumerate(allele_freq):
    #         pe = pe + pi * (1 - pi) ** 2
    #         for pj in allele_freq[(i + 1):]:
    #             pe = pe -(pi ** 2) * (pj ** 2) * (4 - 3 * pi - 3 * pj)


# pd and pe are DataFrame
def cpd(pd1):
    return 1 - (1 - pd1).prod()


def cpe(pe):
    return 1 - (1 - pe).prod()


def dmp(pm, n_entry):
    cpm = pm.prod()
    if cpm >= 1.0 / n_entry:
        return 1 - (1 - cpm) ** n_entry
    else:
        return cpm * n_entry
