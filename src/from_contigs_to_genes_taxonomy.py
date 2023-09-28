# /usr/bin/env python3
import os
import argparse
import sys


usage= 'usage: python from_contigs_to_genes_taxonomy.py [-p] [-t] [-o]'
description = 'This program retreives the contigs taxonomy affiliation and transfered to corresponding representative genes'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument ('-p', metavar="proteins mix_assembly", dest= 'p', help='path to representative protein file - mix_assembly', required = True)
parser.add_argument ('-t', metavar="contigs taxonomy table", dest= 't', help='path to contigs taxonomy table', required = True)
parser.add_argument ('-o', metavar="representative genes taxonomy output file name", dest= 'o', help='output file, representative genes taxonomy affiliation table', required = True)
parser.add_argument ('-k', metavar="rep_genes_taxonomy_krona.txt", dest= 'k', help='output file, representative genes taxonomy affiliation krona file', required = True)

args = parser.parse_args()

def contigs_taxonomy(input_file):
    contig_taxo = {}
    with open(input_file) as fin:
        for line in fin:
            line=line.rstrip()
            LS=line.split("\t")
            if LS[3] != "unclassified":
                contig_taxo[LS[0]]=LS[8] #0 contig name, 8 taxonomy affiliation
            else:
                contig_taxo[LS[0]]=LS[3] # 3 unclassified
    return(contig_taxo)



C_taxo=contigs_taxonomy(args.t)

with open(args.p, "r") as fin, open(args.o, "w") as fout, open(args.k, "w") as foutk:
    for line in fin:
        line=line.rstrip()
        if line.startswith(">"):
            gene=line.split()[0][1:]
            contig="_".join(gene.split("_")[:-1])
            taxo=C_taxo[contig]
            taxo_list=taxo.split(";")
            taxo_tab="\t".join([t for t in taxo_list if t != "-_cellular organisms"])
            print("{}\t{}".format(gene, taxo), file=fout)
            print("{}\t{}".format(1, taxo_tab), file=foutk)
