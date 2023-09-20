#!/usr/bin/env python3

import os
import argparse


usage = 'python clean_mix_assembly_genes.py -g -p -m -f -o -u'
description = 'This program removes false postive genes predicted by prodigal with RFAM anotation related to rRNA or tRNA; and the genes related to t.thermophilus'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument('-g', dest='g', help='gene sequences fasta file', required = True)
parser.add_argument('-p', dest='p', help='protein sequences fasta file', required = True)
parser.add_argument('-m', dest='m', help='gene names from special list of gene to be removed, e.g T thermophilus', default="")
parser.add_argument('-f', dest='f', help='gene names false positives rfam', default='list_rRNA_genes_in_mix.txt')
parser.add_argument('-o', dest='o', help='output file, genes', required = True)
parser.add_argument('-u', dest='u', help='output file, proteins', required = True)
args = parser.parse_args()

def get_false_positives(file, removable):
    with open(file, "r") as fin:
        for line in fin:
            line=line.rstrip()
            if "::" in line:
                line=linne.split("::")[1]
            removable.add(line)
    return removable


def print_file(filein, fileout, exclude):
    copy = False
    with open(filein, "r") as flin, open(fileout, "w") as flout:
        for line in flin:
            line = line.rstrip()
            if line.startswith(">"):
                #>SRR3727505::SRR3727505_k119_2_1 # 3 # 65 # -1 # ID=1_1;partial=10;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.571
                header=line
                id=line.split()[0][1:]
                id=id.split("::")[1]
                if not id in exclude:
                    copy=True
            else:
                if copy:
                    print("{}\n{}".format(header, line), file=flout)
                    copy = False

set_of_genes_to_remove=set()
if os.path.isfile(args.m):
    print("Removing specific genes")
    set_of_genes_to_remove=get_false_positives(args.m, set_of_genes_to_remove)
set_of_genes_to_remove=get_false_positives(args.f, set_of_genes_to_remove)
print_file(args.g, args.o, set_of_genes_to_remove)
print_file(args.p, args.u, set_of_genes_to_remove)

