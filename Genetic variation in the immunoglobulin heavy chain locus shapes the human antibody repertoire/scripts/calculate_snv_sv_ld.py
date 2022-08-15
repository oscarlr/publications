#!/bin/env python
import sys
#import scipy
from scipy import stats
from collections import Counter
import numpy as np

def read_variants(fn):
    snps = {}
    with open(fn,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            if line[0] == "pos":
                samples = line[1:]
                continue
            variant = line[0]
            genotypes = line[1:]
            snps[variant] = {}
            for sample,genotype in zip(samples,genotypes):
                if genotype == "NA":
                    continue
                snps[variant][sample] = genotype #float(genotype)
    return snps

def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False

def get_cramer_v(data):
    data = np.array(data)
    X2 = stats.chi2_contingency(data, correction=False)
    n = np.sum(data)
    minDim = min(data.shape)-1

    pval = "."
    V = "."
    if minDim > 0:
        #calculate Cramer's V 
        V = np.sqrt((X2[0]/n) / minDim)
        pval = X2[1] 
    return (V,pval)
    
bi_alleleic_snps = read_variants(sys.argv[1])
svs_genotypes = read_variants(sys.argv[2])

out = ["snp","sv","slope","intercept","r_value","p_value","std_err"]
print("\t".join(map(str,out)))
        
for snp in bi_alleleic_snps:
    for sv in svs_genotypes:
        snp_gs = []
        sv_gs = []
        all_floats = True
        for sample in bi_alleleic_snps[snp]:
            if sample not in svs_genotypes[sv]:
                continue
            value = bi_alleleic_snps[snp][sample]
            is_float = isfloat(value)
            if not is_float:
                all_floats = False
            value = svs_genotypes[sv][sample]
            is_float = isfloat(value)
            if not is_float:
                all_floats = False            
        for sample in bi_alleleic_snps[snp]:
            if sample not in svs_genotypes[sv]:
                continue
            v1 = bi_alleleic_snps[snp][sample]
            v2 = svs_genotypes[sv][sample]
            if all_floats:
                v1 = float(v1)
                v2 = float(v2)
            snp_gs.append(v1)
            sv_gs.append(v2)
        if all_floats:
            slope, intercept, r_value, p_value, std_err = stats.linregress(snp_gs, sv_gs)
        else:
            occs = []
            unique_snps = set()
            unique_svs = set()
            for snp_g,sv_g in zip(snp_gs,sv_gs):
                occs.append((snp_g,sv_g))
                unique_snps.add(snp_g)
                unique_svs.add(sv_g)
            counts = Counter(occs)
            matrix = []
            for snp_g in unique_snps:
                snp_svs_counts = []
                for sv_g in unique_svs:
                    snp_svs_counts.append(counts[(snp_g,sv_g)])
                matrix.append(snp_svs_counts)
            r_value,p_value = get_cramer_v(matrix)
            slope = "."
            intercept = "."
            std_err = "."
        out = [snp,sv,slope, intercept, r_value, p_value, std_err]
        print("\t".join(map(str,out)))
