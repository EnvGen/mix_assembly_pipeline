# /usr/bin/env python3
import os
import argparse
import sys
import gzip
import io

usage= 'usage: from_faa_to_fna.py [-id] [-r] [-o]'
description = 'This program converts the protein representative assembled genes to nucleotide representative assembled genes'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument ('-id', metavar="fna.gz files", dest= 'id', help='Path to directory with nucleotide sequenes, .fna.gz', required = True)
parser.add_argument ('-r', metavar="rep proteins.faa.gz", dest= 'r', help='path to representative proteins of the samples', required = True)
parser.add_argument ('-o', metavar="output file name", dest= 'o', help='file name output, containing the representative nucleotide sequences, e.g ind_assembly.final.fna (not interleaved)', required = True)

args = parser.parse_args()

def get_rep_proteins(faain):
    ids=set()
    with gzip.open(faain, 'rb') as input_file:
        with io.TextIOWrapper(input_file, encoding='utf-8') as fin:
            for line in fin:
                line=line.rstrip()
                if line.startswith(">"):
                #>P6071_510_k141_33_1 # 2 # 631 # -1 # ID=32_1;partial=11;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.297
                    ids.add(line.split()[0][1:]) #P6071_510_k141_33_1
    return(ids)

list_genes = get_rep_proteins(args.r) #list of representative genes
files=[ d for d in os.listdir(args.id) if d.endswith(".fna.gz") ]

for name in files: #ok in this case, because all samples have representative genes

#extracting the representative nucleotide sequences from the sample

    with gzip.open(os.path.join(args.id,name), 'rb') as input_file:
        with io.TextIOWrapper(input_file, encoding='utf-8') as seqs, open(args.o, "a") as fout:

            copy = False
            first_time = True
            for line1 in seqs:
                line1 = line1.rstrip()
                if line1.startswith(">"):
                    gene_id = line1.split()[0][1:]
                    if gene_id in list_genes: #if it is a representative gene
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
            print('\n', end="", file=fout)
