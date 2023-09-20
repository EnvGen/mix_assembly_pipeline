#!/usr/bin/env python3

import os
import argparse
import re
import gzip
import io

usage = 'create_GTDB_fasta_db.py -d -b -a -o'
description = 'This program creates a fasta DB under the required mmseqs2 formart for GTDB'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument('-d', dest='d', help='directory GTDB proteins', default="protein_faa_reps")
parser.add_argument('-b', dest='b', help='bacteria tax id file .tsv', default="bac120_taxonomy.tsv")
parser.add_argument('-a', dest='a', help='archea tax id file .tsv', default="ar122_taxonomy.tsv")
parser.add_argument('-o', dest='o', help='output non-interleaved fasta DB', default="TESTproteins.faa")
args = parser.parse_args()

taxo_arch={}
taxo_bact={}
with open(args.b, "r") as finb:
    for line in finb:
        line=line.rstrip()
        #RS_GCF_014075335.1      d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia flexneri
        line=line.split()
        taxo_bact[line[0]]=line[1]

with open(args.a, "r") as fina:
    for linea in fina:
        linea=linea.rstrip()
        linea=linea.split()
        taxo_arch[linea[0]]=linea[1]

files_arch=[f for f in os.listdir(os.path.join(args.d, "archaea")) if f.endswith(".faa.gz")]
files_bact=[f for f in os.listdir(os.path.join(args.d, "bacteria")) if f.endswith(".faa.gz")]

with open(args.o, "w") as fout:
    for file in files_arch:
        id=file.split("_protein")[0]
        if id in taxo_arch:
            with gzip.open(os.path.join(args.d,"archaea", file), "r") as input_file:
                with io.TextIOWrapper(input_file, encoding='utf-8') as fin:
                    copy=False
                    for line1 in fin:
                        line1=line1.rstrip()
                        if line1.startswith(">"):
                            if copy:
                                print(header, file=fout)
                                print(seq, file=fout)
                                copy=True

                            header=line1.split()[0][1:]
                            header=">"+id+"~"+header+" "+taxo_arch[id]
                            seq=""
                        else:
                            seq+=line1
    #printing out last sequence
                    print(header, file=fout)
                    print(seq, file=fout)

    for fileb in files_bact:
        idb=fileb.split("_protein")[0] #GB_GCA_003266145.1_protein.faa
        if idb in taxo_bact:
            with gzip.open(os.path.join(args.d,"bacteria", fileb), "r") as input_file2:
                with io.TextIOWrapper(input_file2, encoding='utf-8') as fin2:                   
                    copy=False
                    for lineb in fin2:
                        lineb=lineb.rstrip()
                        if lineb.startswith(">"):
                            if copy:
                                print(headerb, file=fout)
                                print(seqb, file=fout)
                                copy=True

                            headerb=lineb.split()[0][1:]
                            headerb=">"+idb+"~"+headerb+" "+taxo_bact[idb]
                            seqb=""
                        else:
                            seqb+=lineb
    #printing out last sequence
                    print(headerb, file=fout)
                    print(seqb, file=fout)
