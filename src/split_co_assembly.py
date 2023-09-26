# /usr/bin/env python3
import os
import argparse
import sys
import io
import gzip
import math

usage= 'usage: python split_co_assembly.py [-i] [-c] [-o] [-n] [-p]'
description = 'This program splits a fasta (interleaved or not) file (.fa or .fa.gz) into N small files'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument ('-i', metavar="input file", dest= 'i', required = True)
parser.add_argument ('-o', metavar="output directory name", dest= 'o', required = True)
parser.add_argument ('-c', metavar="Number of chuncks", dest= 'c', required = True)
parser.add_argument ('-n', metavar="Number of contigs", dest= 'n', required = True)
parser.add_argument ('-p', metavar="suffix", dest= 'p', default = "_co_assembly_contig.fasta")
args = parser.parse_args()

def list_of_seqs_in_files():
    seqs_file_x=math.ceil(int(args.n)/int(args.c))
    seqs_file_y=math.floor(int(args.n)/int(args.c))
    counter_x=int(args.n) - seqs_file_y*int(args.c)
    counter_y=int(args.c) - counter_x
    myList1 = [ seqs_file_y  for i in range(counter_y) ]
    myList2 = [ seqs_file_x  for i in range(counter_x) ]
    myList1.extend(myList2)
    return(myList1)

if not os.path.exists(args.o):
    os.makedirs(args.o)

sequence = ''
seqs_file=list_of_seqs_in_files()
counter=0
Nseqs=0
chunk_files={}

if args.i.endswith("gz"):
  with gzip.open(args.i, 'rb') as input_file:
    with io.TextIOWrapper(input_file, encoding='utf-8') as fin:
        for line in fin:
            line = line.rstrip()
            if line.startswith(">"):
                if sequence != '':
                   chunk_files[id]=sequence
                   sequence = ''
                id=line
                Nseqs += 1
            else:
               sequence += line
            if Nseqs > seqs_file[counter] :
               with open(os.path.join(args.o,str(counter)+args.p), "w") as fout:
                   for ids,cont in chunk_files.items():
                       print("{}\n{}".format(ids,cont), file=fout)
               counter += 1
               Nseqs = 1
               chunk_files={}
        #printig info from the last sequence
        chunk_files[id]=sequence
        with open(os.path.join(args.o,str(counter)+args.p), "w") as fout:
            for ids,cont in chunk_files.items():
                print("{}\n{}".format(ids,cont), file=fout)

else:
    with open(args.i, 'r') as fin:
        for line in fin:
            line = line.rstrip()
            if line.startswith(">"):
                if sequence != '':
                   chunk_files[id]=sequence
                   sequence = ''
                id=line
                Nseqs += 1
            else:
               sequence += line
            if Nseqs > seqs_file[counter] :
               with open(os.path.join(args.o,str(counter)+args.p), "w") as fout:
                   for ids,cont in chunk_files.items():
                       print("{}\n{}".format(ids,cont), file=fout)
               counter += 1
               Nseqs = 1
               chunk_files={}
        #printig info from the last sequence
        chunk_files[id]=sequence
        with open(os.path.join(args.o,str(counter)+args.p), "w") as fout:
            for ids,cont in chunk_files.items():
                print("{}\n{}".format(ids,cont), file=fout)
