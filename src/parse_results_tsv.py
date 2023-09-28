#!/usr/bin/env python3

import os
import argparse
import sys

usage = 'parse_results_tsv.py -i -o'
description = 'This program selects the best hit from several database and prints them out in a tsv file'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument('-i', dest='i', help='list of path to .tsv files', required=True)
parser.add_argument('-o', dest='o', help='output file name', default="TEST_taxonomy.tsv")
args = parser.parse_args()


files = [f for f in args.i.split(",") ]

query_hits={}
for file in files:
    print("extracting hits from {}".format(file))
    with open(file, "r") as f1:
#        counter = 0
        for line in f1:
            line=line.rstrip()
            #P1994_103_k141_7697     414004  strain  Cenarchaeum symbiosum A 4       1       1       1.000   -_cellular organisms;d_Archaea;-_TACK group;p_Thaumarchaeota;o_Cenarchaeales;f_Cenarchaeaceae;g_Cenarchaeum;s_Cenarchaeum symbiosum;-_Cenarchaeum symbiosum A
            LS=line.split("\t")
            query=LS[0]
            support_received=float(LS[7]) #the support received
            Nreads=float(LS[6])
            if query in query_hits:
                if support_received > float(query_hits[query].split("\t")[7]) and Nreads > float(query_hits[query].split("\t")[6]):
                    query_hits[query]=line
            else:
                query_hits[query]=line



with open(args.o, "w") as fu:
    for l in query_hits.values():
                print(l, file=fu)
