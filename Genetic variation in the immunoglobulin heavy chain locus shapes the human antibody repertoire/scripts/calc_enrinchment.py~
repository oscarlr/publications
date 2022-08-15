#!/bin/env python
import sys
import numpy as np
from pybedtools import BedTool
from scipy.stats import fisher_exact

def get_elements(element,element_col):
    elements = set()
    with open(element,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            elements.add(line[element_col - 1])
    return elements

def filter_for_element(feature,col,e):
    return feature[col-1] == e

def get_total_snps(vcf):
    count = 0
    with open(vcf,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            if line[0] == "chr14":
                count += 1
    return count
    
guqtls_vcf = sys.argv[1]
common_vcf = sys.argv[2]
ccre = sys.argv[3]
tfbs = sys.argv[4]

data_for_test = {}

for vcf in [guqtls_vcf,common_vcf]:
    total_snps = get_total_snps(vcf)
    for type_,element,element_col,overlap_col in [["ccre",ccre,6,169],["tfbs",tfbs,4,167]]:
        elements = get_elements(element,element_col)
        overlap = BedTool(vcf).intersect(BedTool(element),wo=True,header=True)
        for e in elements:
            count = overlap.filter(filter_for_element,col=overlap_col,e=e).count()
            if e not in data_for_test:
                data_for_test[e] = [type_]
            data_for_test[e].append(count)
            data_for_test[e].append(total_snps)


for e in data_for_test:
    type_ = data_for_test[e][0]
    all_guqtl_snps = data_for_test[e][2]
    all_snps_in_e = data_for_test[e][3]
    all_snps = data_for_test[e][4]
    
    num_guqtl_snps_in_e = data_for_test[e][1]
    num_guqtl_snps_not_in_e = all_guqtl_snps - num_guqtl_snps_in_e
    num_non_guqtl_snps_in_e = all_snps_in_e - num_guqtl_snps_in_e
    num_non_guqtl_snps_not_in_e = all_snps - num_guqtl_snps_in_e - num_guqtl_snps_not_in_e - num_non_guqtl_snps_in_e
    
    table = np.array([[num_guqtl_snps_in_e, num_guqtl_snps_not_in_e],
                      [num_non_guqtl_snps_in_e, num_non_guqtl_snps_not_in_e]])
    #print(table)
    oddsr, p = fisher_exact(table, alternative='greater')
    out = [e,num_guqtl_snps_in_e,num_guqtl_snps_not_in_e,num_non_guqtl_snps_in_e,num_non_guqtl_snps_not_in_e,oddsr,p,type_]
    print("\t".join(map(str,out)))

