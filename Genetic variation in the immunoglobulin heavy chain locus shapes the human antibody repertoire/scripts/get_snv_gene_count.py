#!/bin/env python
import sys

variants = {}

with open(sys.argv[1],'r') as fh:
    for line in fh:
        line = line.rstrip().split('\t')
        if line[20] != "True":
            continue
        variant = line[0]
        gene = line[1]
        if variant not in variants:
            variants[variant] = {}
        if gene not in variants[variant]:
            variants[variant][gene] = 0
        variants[variant][gene] += 1

with open(sys.argv[2],'w') as outfh:
    for variant in variants:
        out = [variant,len(variants[variant])]
        outfh.write("%s\n" % "\t".join(map(str,out)))
