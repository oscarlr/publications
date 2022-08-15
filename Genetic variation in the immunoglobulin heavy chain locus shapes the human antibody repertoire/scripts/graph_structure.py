#!/bin/env python
import sys
import networkx as nx
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import math

matplotlib.rcParams['figure.figsize'] = 12, 12

fn = "%s/%s" % (sys.argv[2],sys.argv[1])
edges_fn = "%s/results/edges.txt" % sys.argv[2]

snps = {}

with open(fn,'r') as fh:
    for line in fh:
        line = line.rstrip().split('\t')
        snp = line[0]
        gene = line[1]

        if snp not in snps:
            snps[snp] = set()

        snps[snp].add(gene)

G = nx.Graph()

for snp in snps:
    genes = snps[snp]
    added_pairs = []
    for g1 in genes:
        for g2 in genes:
            if g1 == g2:
                continue
            if sorted([g1,g2]) in added_pairs:
                continue
            if G.has_edge(g1,g2):
                G[g1][g2]['weight'] += 1
            else:
                G.add_edge(g1,g2,weight=1)
                #G[g1][g2]['weight'] = 0
            added_pairs.append(sorted([g1,g2]))

with open(edges_fn,'w') as outfh:
    for (u, v, wt) in G.edges.data('weight'):
        outfh.write("%s\t%s\t%s\n" % (u,v,wt))        


G1 = nx.Graph([(u,v,d) for (u,v,d) in  G.edges(data=True) if d['weight'] >= 3])

nx.write_gml(G1,"%s/figures/complete_graph.gml" % sys.argv[2])

#sys.exit()

S = [G1.subgraph(c).copy() for c in nx.connected_components(G1)]

for i,s1 in enumerate(S):    
    pos=nx.spring_layout(s1)   #nx.spectral_layout(s1)  #nx.spring_layout(s1,k=20) # pos = nx.nx_agraph.graphviz_layout(G)
    #pos = nx.circular_layout(s1)
    nx.draw_networkx(s1,pos)
    labels = nx.get_edge_attributes(s1,'weight')
    nx.draw_networkx_edge_labels(s1,pos,edge_labels=labels)
    plt.tight_layout()
    plt.savefig("%s/figures/complete_graph_c_%i.png" % (sys.argv[2],i))
    plt.clf()

#sys.exit()

for i,nodes in enumerate(list(nx.find_cliques(G1))):
    c = G1.subgraph(nodes)
    pos=nx.spring_layout(c) # pos = nx.nx_agraph.graphviz_layout(G)
    nx.draw_networkx(c,pos)
    labels = nx.get_edge_attributes(c,'weight')
    nx.draw_networkx_edge_labels(c,pos,edge_labels=labels)
    plt.tight_layout()
    plt.savefig("%s/figures/clique_%s.png" % (sys.argv[2],i))
    plt.clf()

sys.exit()

d_genes = set()
v_genes = set()

for j in [3]: #[3,5,10,50]:
    G_weight = nx.Graph([(u,v,d) for (u,v,d) in  G.edges(data=True) if d['weight'] >= j])
    for i,nodes in enumerate(list(nx.find_cliques(G_weight))):
        c = G.subgraph(nodes)
        for n in nodes:
            if "IGHV" in n:
                v_genes.add(n)
            if "IGHD" in n:
                d_genes.add(n)
        pos=nx.spring_layout(c) # pos = nx.nx_agraph.graphviz_layout(G)
        nx.draw_networkx(c,pos)
        labels = nx.get_edge_attributes(c,'weight')
        nx.draw_networkx_edge_labels(c,pos,edge_labels=labels)
        
        plt.savefig("graphs_min_%s/complete_graph_%s.png" % (j,i))
        plt.clf()

print(len(v_genes))
print(len(d_genes))
