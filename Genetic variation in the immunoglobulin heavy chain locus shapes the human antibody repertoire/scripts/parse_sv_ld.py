#!/bin/env python
import sys

variants = {}

with open(sys.argv[1],'r') as fofh:
    for line in fofh:
        line = line.rstrip().split('\t')
        sv_type = line[0]
        fn = line[1]
        with open(fn,'r') as fh:
            for line in fh:
                line = line.rstrip().split('\t')
                variant = line[0]
                sv = line[1]
                pval = line[5]
                r = line[4]
                if "." == pval:
                    continue
                if pval == "p_value":
                    continue
                pval = float(pval)
                r = float(r)
                
                if pval > 0.05:
                    continue

                if variant not in variants:
                    variants[variant] = {}

                variants[variant][(sv,sv_type)] = r

with open(sys.argv[2],'r') as fofh:
    for line in fofh:
        line = line.rstrip().split('\t')
        variant = line[0]
        out = []
        if variant in variants:
            for sv,sv_type in variants[variant]:
                corr = "r"
                if sv_type == "anova":
                    corr = "V"
                r = variants[variant][(sv,sv_type)]
                # 7-mCNV (r = 0.382)
                out.append("%s (%s = %.3f)" % (sv,corr,r))
        if len(out) == 0:
            out = [variant,"."]
        else:
            out = [variant,", ".join(out)]
        print("\t".join(out))
