#!/usr/bin/env python3

import os
import argparse
import sys


usage = 'krona_clean_results_tsv.py -i -o'
description = 'This program creates input file for krona'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument('-i', dest='i', help='path to .tsv file', required=True)
parser.add_argument('-o', dest='o', help='output file name', default="TEST_taxonomy_krona.txt")
args = parser.parse_args()

with open(args.i, "r") as f1, open(args.o, "w") as fout:
	for line in f1:
		line=line.rstrip()
		LS=line.split("\t")
		taxo="unclassified"

		if float(LS[7]) > 0.0:
			taxo1=LS[8].split(";")
			taxo="\t".join([t for t in taxo1 if t != "-_cellular organisms"])
		print("{}\t{}".format(1, taxo), file=fout)
