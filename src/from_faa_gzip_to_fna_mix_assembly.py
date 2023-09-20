# /usr/bin/env python3
import os
import argparse
import sys
import gzip
import io


usage= 'usage: from_faa_to_fna_mix_assembly.py [-id] [-gd] [-od] [-o]'
description = 'This program converts the protein representative assembled genes to nucleotide representative assembled genes'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument ('-i', metavar="ind_assembly fna.gz files", dest= 'i', help='Path to representative gene sequence file - ind_assembly, .fna.gz', required = True)
parser.add_argument ('-c', metavar="co_assembly fna.gz files", dest= 'c', help='Path to gene sequence file - co_assembly, .fna.gz', required = True)
parser.add_argument ('-m', metavar="proteins mix_assembly", dest= 'm', help='path to representative protein file - mix_assembly, e.g., DB_clu_rep.fasta', required = True)
parser.add_argument ('-o', metavar="output file name", dest= 'o', help='Path to output file, containing the representative gene sequences - mix_assembly, e.g., mix_assembly.genes.fna (not interleaved)', required = True)
args = parser.parse_args()


# opening file with list of representative proteins from mix_assembly
with open(args.m, "r") as genes:
    list_genes = set()
    for line in genes:
        line = line.rstrip()
        if line.startswith(">"):
            gen = line.split()[0][1:]
            list_genes.add(gen) #storing in a set of representative genes mix_assembly
#print(list_genes)

with open(args.o, "a") as fout:
  with gzip.open(args.i, 'rb') as input_file:
    with io.TextIOWrapper(input_file, encoding='utf-8') as ind_assembly_seqs:
        #extracting the representative nucleotide sequences from ind_assembly
            copy = False
            first_time = True
            add_last_line = False
            for line1 in ind_assembly_seqs:
                line1 = line1.rstrip()
                if line1.startswith(">"):
                    gene_id = line1.split()[0][1:]
                    if gene_id in list_genes: #if it is a representative gene
                        add_last_line = True
                        if first_time:
                            print(line1, file=fout) #printing >gene ID
                            first_time = False
                        else:
                            print('\n', end="", file=fout)
                            print(line1, file=fout) #printing >gene ID
                        copy = True
                    else:
                        copy = False
                else:
                    if copy:
                        print(line1, end="", file=fout) #printing nucleotide sequence (interleaved or not)

            #last line before adding new sequences
            if add_last_line:
                   print('\n', end="", file=fout)

  with gzip.open(args.c, 'rb') as co_assembly_input_file:
     with io.TextIOWrapper(co_assembly_input_file, encoding='utf-8') as co_assembly_seqs:
        #extracting the representative nucleotide sequences from co_assembly
            add_last_line = False
            copy = False
            first_time = True
            for line2 in co_assembly_seqs:
                line2 = line2.rstrip()
                if line2.startswith(">"):
                    gene_id = line2.split()[0][1:]
                    if gene_id in list_genes: #if it is a representative gene
                        add_last_line = True
                        if first_time:
                            print(line2, file=fout) #printing >gene ID
                            first_time = False
                        else:
                            print('\n', end="", file=fout)
                            print(line2, file=fout) #printing >gene ID
                        copy = True
                    else:
                        copy = False
                else:
                    if copy:
                        print(line2, end="", file=fout) #printing nucleotide sequence (interleaved or not)

            #last line before adding new sequences
            if add_last_line:
                     print('\n', end="", file=fout)
