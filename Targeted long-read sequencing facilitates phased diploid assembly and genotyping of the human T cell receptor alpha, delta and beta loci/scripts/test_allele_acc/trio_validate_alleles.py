#!/bin/env python
import sys

proband_alleles = sys.argv[1]
parent1_alleles = sys.argv[2]
parent2_alleles = sys.argv[3]

def read_genotypes(alleles):
    genotypes = {}
    with open(alleles,'r') as alleles_fh:
        for line in alleles_fh:
            if "#" in line:
                continue
            line = line.rstrip().split('\t')
            chrom = line[0]
            start = line[1]
            end = line[2]
            gene = line[3]
            hap = line[4]
            category = line[5]
            allele = line[6]
            allele_obs = line[7]
            if (chrom,start,end,gene) not in genotypes:
                genotypes[(chrom,start,end,gene)] = []
            genotypes[(chrom,start,end,gene)].append([hap,category,allele,allele_obs])
    return genotypes

def get_alleles(person):
    alleles = set()
    for hap in person:
        allele = hap[2]
        alleles.add(allele)
    return alleles
        
proband_genotypes = read_genotypes(proband_alleles)
parent1_genotypes = read_genotypes(parent1_alleles)
parent2_genotypes = read_genotypes(parent2_alleles)


for (chrom,start,end,gene) in proband_genotypes:
    proband_alleles = get_alleles(proband_genotypes[(chrom,start,end,gene)])
    parent1_alleles = set()
    parent2_alleles = set()
    if (chrom,start,end,gene) in parent1_genotypes:
        parent1_alleles = get_alleles(parent1_genotypes[(chrom,start,end,gene)])
    if (chrom,start,end,gene) in parent2_genotypes:        
        parent2_alleles = get_alleles(parent2_genotypes[(chrom,start,end,gene)])
    pass_mend = False
    parent1_allele = None
    parent2_allele = None
    for proband_allele1 in proband_alleles:
        for proband_allele2 in proband_alleles:
            if len(proband_alleles) == 1:
                if proband_allele1 in parent1_alleles:
                    if proband_allele1 in parent2_alleles:
                        pass_mend = True
                        parent1_allele = proband_allele1
                        parent2_allele = proband_allele1
                        break
            else:
                if proband_allele1 == proband_allele2:
                    continue
                if proband_allele1 in parent1_alleles:
                    if proband_allele2 in parent2_alleles:
                        pass_mend = True
                        parent1_allele = proband_allele1
                        parent2_allele = proband_allele2
                        break
                if proband_allele2 in parent1_alleles:
                    if proband_allele1 in parent2_alleles:
                        pass_mend = True
                        parent1_allele = proband_allele1
                        parent2_allele = proband_allele2
                        break                    
        if pass_mend:
            break
    out = [chrom,start,end,gene,",".join(list(proband_alleles)),parent1_allele,parent2_allele,pass_mend]
    print("\t".join(map(str,out)))
