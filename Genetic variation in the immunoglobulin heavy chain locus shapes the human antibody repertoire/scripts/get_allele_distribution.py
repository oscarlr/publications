#!/bin/env python
import sys
import pandas as pd
from bioinfokit.analys import stat
from sklearn.preprocessing import LabelEncoder

def read_snps_genes(fn):
    snps_genes = []
    with open(fn,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            if line[22] != "1":
                continue
            if line[20] != "True":
                continue
            variant = line[0]            
            gene = line[1]
            snps_genes.append([variant,gene])
    return snps_genes

def read_variant_genotypes(fn):
    variant_genotypes = {}
    with open(fn,'r') as fh:
        for line in fh:
            variant_fn = line.rstrip().split('\t')[1]
            samples = None
            with open(variant_fn,'r') as variant_fh:
                for variant_line in variant_fh:
                    variant_line = variant_line.rstrip().split('\t')
                    if samples == None:
                        samples = variant_line[1:]
                        continue
                    variant = variant_line[0]
                    genotypes = variant_line[1:]
                    variant_genotypes[variant] = {}
                    for sample,genotype in zip(samples,genotypes):
                        variant_genotypes[variant][sample] = genotype
    return variant_genotypes

def read_alleles(fn,novel_alleles):
    alleles = {}
    with open(fn,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            sample = line[0]
            sample_fn = line[1]
            if sample not in alleles:
                alleles[sample] = {}
            with open(sample_fn,'r') as sample_fh:
                for line1 in sample_fh:
                    line1 = line1.rstrip().split('\t')
                    gene = line1[3]
                    hap = line1[4]
                    category = line1[5]
                    allele = line1[6]
                    if allele == "allele":
                        continue
                    if gene not in alleles[sample]:
                        alleles[sample][gene] = {}
                    if hap not in alleles[sample][gene]:
                        alleles[sample][gene][hap] = []
                    if category != "NOVEL":
                        allele = int(float(allele))
                        if allele < 10:
                            allele = "0%s" % allele
                        else:
                            allele = "%s" % allele
                        alleles[sample][gene][hap].append(allele)
                    else:
                        if gene in novel_alleles:
                            if allele in novel_alleles[gene]:
                                allele_name = novel_alleles[gene][allele]
                                novel_allele_name = allele_name.split('=')[2].split('_')[0]
                                alleles[sample][gene][hap].append(novel_allele_name)
    return alleles

def get_allele_genotype(hap_alleles):
    genotype = None
    if "0" in hap_alleles:
        if len(hap_alleles) == 1:
            if len(hap_alleles["0"]) == 0:
                return genotype
            genotype = ",".join(sorted(hap_alleles["0"]))
            genotype = "*%s/*%s" % (genotype,genotype)
            return genotype
    if "1" in hap_alleles:
        if "2" in hap_alleles:
            if len(hap_alleles["1"]) == 0:
                return genotype
            if len(hap_alleles["2"]) == 0:
                return genotype
            hap1_allele = ",".join(sorted(hap_alleles["1"]))
            hap2_allele = ",".join(sorted(hap_alleles["2"]))
            alleles = sorted([hap1_allele,hap2_allele])            
            genotype = "*%s/*%s" % (alleles[0],alleles[1])
            return genotype
    return genotype

def get_pval(matrix):
    res = stat()
    #res.tukey_hsd(df=d_melt, res_var='value', xfac_var=['Genotype','years'], anova_model='value ~ C(Genotype) + C(years) + C(Genotype):C(years)')
    #res.tukey_summary.head()
    matrix = pd.DataFrame(matrix,columns=["genotype","alleles","count"])
    print(matrix)
    print(matrix.dtypes)
    #try:
    res.tukey_hsd(df=matrix, res_var='count', xfac_var=['genotype','alleles'])
    res.tukey_summary.head()
    #except ValueError:
    #pass

def read_novel_alleles(fn):
    novel_alleles = {}
    with open(fn,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            if line[0] == 'Allele name':
                continue
            allele_name = line[0]
            gene_name = line[8]
            allele_seq = line[9]
            if gene_name not in novel_alleles:
                novel_alleles[gene_name] = {}
            novel_alleles[gene_name][allele_seq] = allele_name
    return novel_alleles

novel_alleles = read_novel_alleles(sys.argv[4])
snps_genes = read_snps_genes(sys.argv[1])
variant_genotypes = read_variant_genotypes(sys.argv[2])
alleles = read_alleles(sys.argv[3],novel_alleles)



data = {}
for variant,gene in snps_genes:
    for sample in variant_genotypes[variant]:
        if sample not in alleles:
            continue
        if gene not in alleles[sample]:
            continue
        genotype = variant_genotypes[variant][sample]
        hap_alleles = alleles[sample][gene]
        allele_genotype = get_allele_genotype(hap_alleles)
        if allele_genotype == None:
            continue
        if (variant,gene) not in data:
            data[(variant,gene)] = {}
        if (genotype,allele_genotype) not in data[(variant,gene)]:
            data[(variant,gene)][(genotype,allele_genotype)] = 0
        data[(variant,gene)][(genotype,allele_genotype)] += 1

for variant,gene in data:
    # if gene != "IGHV1-2":
    #     continue
    #matrix = []
    all_allele_genotypes = set()
    all_genotypes = set()
    for genotype,allele_genotype in data[(variant,gene)]:
        all_allele_genotypes.add(allele_genotype)
        all_genotypes.add(genotype)
    for genotype in all_genotypes:
        for allele_genotype in all_allele_genotypes:            
            if (genotype,allele_genotype) in data[(variant,gene)]:    
                count = data[(variant,gene)][(genotype,allele_genotype)]
            else:
                count = 0
            out = [variant,gene,genotype,allele_genotype,count]
            print("\t".join(map(str,out)))
        #line = [genotype,allele_genotype,count]
        #matrix.append(line)
    #pval = get_pval(matrix)
    
