#!/bin/env python
import sys

variant_genes = {}
with open(sys.argv[2],'r') as fh:
    for line in fh:
        line = line.rstrip().split('\t')
        if line[0] not in variant_genes:
            variant_genes[line[0]] = []
        variant_genes[line[0]].append(line[1])
        
fn1 = sys.argv[1]
with open(fn1,'r') as fh:
    for line in fh:
        line = line.rstrip().split('\t')
        if line[0] in variant_genes:
            if line[1] in variant_genes[line[0]]:
                print("\t".join(line))

