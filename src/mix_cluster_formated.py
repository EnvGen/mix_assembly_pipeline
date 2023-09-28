#!/usr/bin/env python3

import os
import argparse


usage = 'python mix_cluster_formated.py -i -r -o'
description = 'This program creates a table of representative genes and genes in each cluster'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument('-i', dest='i', help='file from clustering',required=True) #"Cluster/mix_cluster.tsv"
parser.add_argument('-r', dest='r', help='representative genes.faa',required=True)
parser.add_argument('-o', dest='o', help='output table rfam',required=True)
args = parser.parse_args()


rep_genes=set()
with open(args.r, "r") as finr:
    for line in finr:
        line=line.rstrip()
        if line.startswith(">"):
            id=line.split()[0][1:]
            rep_genes.add(id)


dic_rep={}
with open(args.i,"r") as fin:
    for line in fin:
        line=line.rstrip()
        line=line.split("\t")
        rep=line[0]
        if rep in rep_genes:
            if rep in dic_rep:
                dic_rep[rep].append(line[1])
            else:
                dic_rep[rep]=[line[1]]

#print(rep_genes)
#print(dic_rep)
with open(args.o, "w") as fout:
    print("#Rep_genes\tGenes_in_cluster", file=fout)
    for r,l in dic_rep.items():
        print("{}\t{}".format(r,",".join(l)), file=fout)
