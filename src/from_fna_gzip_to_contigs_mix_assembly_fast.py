# /usr/bin/env python3
import os
import argparse
import sys
import gzip
import io



usage= 'usage: from_fna_gzip_to_contigs_mix_assembly.py [-id] [-gd] [-od] [-o]'
description = 'This program retreives the nucleotide representative assembled contigs from the representative assembled genes ID'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument ('-l', metavar="List of ind assembly files", dest= 'l', required = True)
parser.add_argument ('-c', metavar="co_assembly file", dest= 'c', help='Path to assembly file - co_assembly, co_assembly.large.normalized.contigs.renamed.fa.gz', required = True)
parser.add_argument ('-m', metavar="proteins mix_assembly", dest= 'm', help='path to representative protein file - mix_assembly, e.g., DB_clu_rep.fasta', required = True)
parser.add_argument ('-o', metavar="output file name", dest= 'o', help='output file, containing the representative contigs sequences - mix_assembly, e.g., mix_assembly.contigs.fna (not interleaved)', required = True)

args = parser.parse_args()


def contigs_ind_assemblies(input_file):
    ind_files = set()
    with open(input_file) as fin:
        for line in fin:
            line=line.rstrip()
            if not line.startswith("#"):
                line=line.split(",")     
                ind_files.add(line[1])
    return(ind_files)



# opening file with list of representative proteins from mix_assembly
with open(args.m, "r") as genes:
    list_contigs_ind = set()
    list_contigs_co = set()
    gene_full_name={}
    for line in genes:
        line = line.rstrip()
        if line.startswith(">"):
            gen = line.split()[0][1:]
            contig = "_".join(map(str,gen.split("_")[:-1]))
            contign=contig.split("::")[1]
            gene_full_name[contign]=contig
            if contig.startswith("CO::"):
                
                list_contigs_co.add(contign) #storing in a set of representative contig mix_assembly
            else:
                
                list_contigs_ind.add(contign) #storing in a set of representative contig mix_assembly

    samples = contigs_ind_assemblies(args.l)

for name in samples:

#extracting the representative conting sequences from the sample
    with open(args.o, "a") as fout:

        with gzip.open(name, 'rb') as input_file:
            with io.TextIOWrapper(input_file, encoding='utf-8') as seqs:

                add_last_line = False
                copy = False
                first_time = True
                for line1 in seqs:
                    line1 = line1.rstrip()
                    if line1.startswith(">"):
                        gene_id = line1.split()[0][1:]
                        if gene_id in list_contigs_ind: #if it is a representative gene
                            list_contigs_ind.remove(gene_id) # Remove contig from list to avoid being printing out more than once.
                            add_last_line = True
                            if first_time:
                                New_header= ">"+str(gene_full_name[gene_id])+" "+" ".join(line1.split()[1:]) 
                                print(New_header, file=fout) #printing >gene ID
                                first_time = False
                            else:
                                print('\n', end="", file=fout)
                                New_header= ">"+str(gene_full_name[gene_id])+" "+" ".join(line1.split()[1:]) 
                                print(New_header, file=fout) #printing >gene ID
                            copy = True
                        else:
                            copy = False
                    else:
                        if copy:
                            print(line1, end="", file=fout) #printing nucleotide sequence (interleaved or not)

                #last line before adding new sequences
                if add_last_line:
                    print('\n', end="", file=fout)


#with open(args.o, "a") as fout:
        with gzip.open(args.c, 'rb') as co_assembly_input_file:
             with io.TextIOWrapper(co_assembly_input_file, encoding='utf-8') as co_assembly_seqs:
                #extracting the representative nucleotide sequences from co_assembly
                    add_last_line = False
                    copy = False
                    first_time = True
                    for line2 in co_assembly_seqs:
                        line2 = line2.rstrip()
                        if line2.startswith(">"):
                            gene_id2 = line2.split()[0][1:]
                            
                            if gene_id2 in list_contigs_co: #if it is a representative gene
                                list_contigs_co.remove(gene_id2) # Remove contig from list to avoid being printing out more than once.
                                add_last_line = True
                                if first_time:
                                    New_header2= ">"+str(gene_full_name[gene_id2])+" "+" ".join(line2.split()[1:]) 
                                    print(New_header2, file=fout) #printing >gene ID
                                    first_time = False
                                else:
                                    print('\n', end="", file=fout)
                                    New_header2= ">"+str(gene_full_name[gene_id2])+" "+" ".join(line2.split()[1:]) 
                                    print(New_header2, file=fout) #printing >gene ID
                                copy = True
                            else:
                                copy = False
                        else:
                            if copy:
                                print(line2, end="", file=fout) #printing nucleotide sequence (interleaved or not)

                    #last line before adding new sequences
                    if add_last_line:
                             print('\n', end="", file=fout)


