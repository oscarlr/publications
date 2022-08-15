#!/bin/env python
import sys
import pandas as pd

def read_genes(fn):
    genes = set()
    with open(fn,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            if line[0] == "variable":
                continue
            gene = line[0]
            genes.add(gene)
    return genes

def read_svs(fn):
    svs = {}
    with open(fn,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            sv = line[1]
            gene = line[0]
            svs[gene] = sv
    return svs

def read_gene_coords(fn):
    gene_coords = {}
    with open(fn,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            gene = line[0]
            start = int(line[2])
            end = int(line[3])
            gene_coords[gene] = [start,end]
    return gene_coords

def get_num_sign_variants(df,gene):
    gene_df = df[df["gene"] == gene]
    gene_df = gene_df[gene_df["sign_corr"] == True]
    return len(gene_df.index)

def get_top_hits(df,gene):
    top_hits = []
    gene_df = df[df["gene"] == gene]
    #gene_df = gene_df[gene_df["sign_corr"] == True]
    #num_sign_variants = len(gene_df.index)
    #if num_sign_variants > 0:
    top_hits = gene_df[gene_df["rank_order"] == 1]
    top_hits = list(top_hits.index)
    return top_hits

def get_variant_types(fofn):
    variant_types = {}
    with open(fofn,'r') as fofh:
        for line in fofh:
            line = line.rstrip().split('\t')
            variant_type = line[0]
            if variant_type == "snps":
                variant_type = "SNV"
            if variant_type in ["large_sv_linear","large_sv_anova"]:
                variant_type = "Large SV"
            if variant_type == "indel_not_tr":
                variant_type = "Indel (not TR)"
            if variant_type == "indel_tr":
                variant_type = "Indel (TR)"
            if variant_type == "sv_not_tr":
                variant_type = "SV (not TR)"            
            eqtl_fn = line[1]
            with open(eqtl_fn,'r') as fh:
                for eqtl_line in fh:
                    eqtl_line = eqtl_line.rstrip().split('\t')
                    variant = eqtl_line[0]
                    variant_types[variant] = variant_type
    return variant_types

def get_distance_between_variant_gene(top_hit,gene_coords):
    distance = "."
    if top_hit.isdigit():
        variant = int(top_hit)
        if variant > gene_coords[0]:
            if variant < gene_coords[1]:
                distance = 0
        distance = min([abs(variant - gene_coords[0]),abs(variant - gene_coords[1])])
    elif "_" in top_hit:
        variants = map(int,top_hit.split("_"))
        for variant in variants:
            if variant > gene_coords[0]:
                if variant < gene_coords[1]:
                    distance = 0
                    break
            new_distance = min([abs(variant - gene_coords[0]),abs(variant - gene_coords[1])])
            if distance != ".":
                if new_distance < distance:
                    distance = new_distance
            else:
                distance = new_distance
    return distance

def read_sv_variant_ld(fn):
    sv_variant_ld = {}
    with open(fn,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            variant = line[0]
            ld = line[1]
            sv_variant_ld[variant] = ld
    return sv_variant_ld

def parse_sv_ld(sv,sv_lds):
    # 2-DEL (r = 0.207), 3-Complex (r = 0.769)
    ld = "."
    if sv_lds != ".":
        sv_ld_parsed = sv_lds.split(',')
        for sv_ld in sv_ld_parsed:
            corr = float(sv_ld.split("=")[1][:-1])
            if corr < .9:
                ld = False
            else:
                ld = True
                break
    return ld

def get_num_hemi(variant_type,variant,gene,df):
    num_hemi = "."
    if variant_type == "SNV":
        df = df.loc[variant]
        df = df[df["gene"] == gene]
        num_hemi = int(df["num_3"]) + int(df["num_4"])
    return num_hemi

def get_value_from_df(variant_type,variant,gene,df,key):
    value = "."
    if variant in df.index:
        df = df.loc[variant]
        df = df[df["gene"] == gene]
        value = float(df[key])
    return value

def get_bool_value_from_df(variant_type,variant,gene,df,key):
    value = "."
    if variant in df.index:
        df = df.loc[variant]
        df = df[df["gene"] == gene]
        value = df.loc[variant,key]
    return value
    
def get_anova_pval(variant_type,variant,gene,df):
    return get_value_from_df(variant_type,variant,gene,df,"pvalue")

def get_linear_beta(variant_type,variant,gene,df):
    return get_value_from_df(variant_type,variant,gene,df,"beta")

def get_linear_r(variant_type,variant,gene,df):
    return get_value_from_df(variant_type,variant,gene,df,"rsquared")

def get_if_sign(variant_type,variant,gene,df):
    return get_bool_value_from_df(variant_type,variant,gene,df,"sign")

def get_if_sign_after_corr(variant_type,variant,gene,df):
    return get_bool_value_from_df(variant_type,variant,gene,df,"sign_corr")

def get_variant_gene_usages(variant,gene,df_variant_usage):
    df = df_variant_usage.loc[variant]
    df = df[df["gene"] == gene]
    means = df.groupby('genotype')['usage'].mean()
    max_usage = means.max()
    min_usage = means.min()
    usages = []
    for i in [0,1,2]:
        usage = "."
        if i in means.index:
            usage = float(means.loc[i])
        usages.append(usage)
    usages.append(max_usage)
    usages.append(min_usage)
    return usages

def get_rank(variant,gene,df):
    rank = "."
    if variant in df.index:
        df = df.loc[variant]
        df = df[df["gene"] == gene]
        if variant in df.index:
            rank = df.loc[variant,"rank_order"]
    return rank

def read_snp_bases(vcffn):
    snp_bases = {}
    with open(vcffn,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            if line[0] != "igh":
                continue
            pos = line[1]
            ref = line[3]
            alt = line[4]
            snp_bases[pos] = [ref,alt]
    return snp_bases

def get_variant_bases(vcf_snp_bases,variant):
    ref = "."
    alt = "."
    if variant in vcf_snp_bases:
        ref,alt = vcf_snp_bases[variant]
    return (ref,alt)

def read_rsids(fn):
    rsids = {}
    with open(fn,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            igh_pos = line[0]
            rsid = line[3]
            rsids[igh_pos] = rsid
    return rsids

def get_rsid(rsids,variant):
    rsid = "."
    if variant in rsids:
        rsid = rsids[variant]
    return rsid
    
genes = read_genes(sys.argv[1])
igm = pd.read_csv(sys.argv[2],sep="\t",header=0,index_col=0)
igg = pd.read_csv(sys.argv[3],sep="\t",header=0,index_col=0)
svs_with_gene = read_svs(sys.argv[4])
variant_types = get_variant_types(sys.argv[5])
gene_coords = read_gene_coords(sys.argv[4])
variant_sv_ld = read_sv_variant_ld(sys.argv[6])
igm_linear = pd.read_csv(sys.argv[7],sep="\t",header=0,index_col=0)
igg_linear = pd.read_csv(sys.argv[8],sep="\t",header=0,index_col=0)
igm_variant_usage = pd.read_csv(sys.argv[9],sep="\t",header=None,index_col=0,names=["variant","gene","usage","sample","genotype","isotype"])
igg_variant_usage = pd.read_csv(sys.argv[10],sep="\t",header=None,index_col=0,names=["variant","gene","usage","sample","genotype","isotype"])
vcf_snp_bases = read_snp_bases(sys.argv[11])
rsids = read_rsids(sys.argv[12])

isotypes = [["M",igm,igm_linear,igm_variant_usage,igg],
            ["G",igg,igg_linear,igg_variant_usage,igm]]

header = ["isotype","gene","gene_start","gene_end","sv","num_sign_variants","top_hit","variant_rsid","variant_type","ref_base","alt_base","dist","sv_ld","high_sv_ld","num_hemi","anova_pval","linear_pval","linear_beta","linear_r","sign","sign_corr","ref_usage","het_usage","alt_usage","max_usage","min_usage","max_min_fold_changes","rank","other_isotype_rank"]
print("\t".join(map(str,header)))

for isotype,df,df_linear,df_variant_usage,df_alt in isotypes:
    for gene in genes:
        sv = svs_with_gene[gene]
        num_sign_variants = get_num_sign_variants(df,gene)
        top_hits = get_top_hits(df,gene)
        # variant_type = "."
        # dist = "."
        # sv_ld = "."
        # sv_gene_and_ld = "."
        # not_sv_gene_and_ld = "."
        # num_hemi = "."
        for top_hit in top_hits:
            variant_type = variant_types[top_hit]
            gene_start,gene_end = gene_coords[gene]
            dist = get_distance_between_variant_gene(top_hit,gene_coords[gene])
            sv_ld = variant_sv_ld[top_hit]
            #sv_gene_and_ld,not_sv_gene_and_ld = parse_sv_ld(sv,sv_ld)
            high_sv_ld = parse_sv_ld(sv,sv_ld)
            num_hemi = get_num_hemi(variant_type,top_hit,gene,df)
            anova_pval = get_anova_pval(variant_type,top_hit,gene,df)
            linear_pval = get_anova_pval(variant_type,top_hit,gene,df_linear)
            linear_beta = get_linear_beta(variant_type,top_hit,gene,df_linear)
            linear_r = get_linear_r(variant_type,top_hit,gene,df_linear)
            sign = get_if_sign(variant_type,top_hit,gene,df)
            sign_corr = get_if_sign_after_corr(variant_type,top_hit,gene,df)
            ref_usage,het_usage,alt_usage,max_usage,min_usage = get_variant_gene_usages(top_hit,gene,df_variant_usage)
            max_min_fold_changes = max_usage/min_usage
            rank = get_rank(top_hit,gene,df)            
            other_isotype_rank = get_rank(top_hit,gene,df_alt)
            ref_base,alt_base = get_variant_bases(vcf_snp_bases,top_hit)
            variant_rsid = get_rsid(rsids,top_hit)
            if gene == "IGHV1-8":
                gene_start = 247475
                gene_end = 247476
            if gene == "IGHV3-9":
                gene_start = 247477
                gene_end = 247478
            out = [isotype,gene,gene_start,gene_end,sv,num_sign_variants,top_hit,variant_rsid,variant_type,ref_base,alt_base,dist,sv_ld,high_sv_ld,num_hemi,anova_pval,
                   linear_pval,linear_beta,linear_r,sign,sign_corr,ref_usage,het_usage,alt_usage,max_usage,min_usage,max_min_fold_changes,
                   rank,other_isotype_rank]
            print("\t".join(map(str,out)))
        # if num_sign_variants == 0:
        #     out = [isotype,gene,sv,num_sign_variants,variant_type,dist,sv_ld]
        #     print("\t".join(map(str,out)))
        
# Isotype
# Gene
# SV gene is in
# Number of significant variants
# Variant
# Variant type
# Ref base
# Alt base
# Min distance between variant and gene
# Significant (p-value < 0.05) r or V to large SVs
# Gene in SV and Variant type = Large SV or r/V > 0.9
# Gene not in SV and Variant type = Large SV or r/V > 0.9
# Number of hemizygotes
# ANOVA p-value
# Linear beta
# Linear R2
# Adjusted Linear R2
# Significant
# Significant after correction
# Usage for ref genotype
# Usage for alt genotype
# Max usage
# Min usage
# Fold change
# IgM p-value rank
# IgG p-value rank
