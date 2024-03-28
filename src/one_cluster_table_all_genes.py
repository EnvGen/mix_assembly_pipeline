#!/usr/bin/env python3

import os
import argparse
import gzip
import io


usage = 'python one_cluster_table_all_genes.py -i -r -o'
description = 'This program combines mix_assebmly clusters and Ind assembly clusters into one table. Column 1 rep mix assembly gene, column 2 genes in the cluster'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument('-i', dest='i', help='file from clustering ind-assembly',required=True)
parser.add_argument('-m', dest='m', help='file from clustering mix-assembly',required=True)
parser.add_argument('-o', dest='o', help='output table',required=True)
args = parser.parse_args()


def read_cluster_file(filein):
    dic_rep={}
    with open(filein,"r") as fin:
        for line in fin:
            line=line.rstrip()
            line=line.split("\t")
            rep=line[0]
            dic_rep[rep]=set(line[1].split(","))
    return dic_rep

print("Reading Ind assembly clusters")
ind_rep_clusters=read_cluster_file(args.i)
print("Reading Mix assembly clusters")
mix_rep_clusters=read_cluster_file(args.m)

print("Printing deconvoluted Mix assembly clusters")
with open(args.o, "w") as fout:
    co_ptrn="CO::"
    for r,l in mix_rep_clusters.items():
        deconv_genes=set()
        co_genes=set([ c for c in l if c.startswith(co_ptrn)])
        ind_genes=set([ i for i in l if not i.startswith(co_ptrn)])
        if len(co_genes) > 0:
            deconv_genes.update(co_genes)
        if len(ind_genes) > 0:
            for j in ind_genes:
                deconv_genes.update(ind_rep_clusters[j])


        print("{}\t{}".format(r,",".join(deconv_genes)), file=fout)

