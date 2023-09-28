#!/usr/bin/env python3

import os
import argparse


usage = 'python scr/merge_mix_assembly_genes_annotation_tables.py -d -e -t'
description = 'This program creates a table resuming all annotation from representative genes'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument('-d', dest='d', help='file dbCAN annotation table',required=True)
parser.add_argument('-p', dest='p', help='file Pfam annotation table',required=True)
parser.add_argument('-e', dest='e', help='file EggNOG annotation', required=True)
parser.add_argument('-o', dest='o', help='output all annotation table',required=True)

args = parser.parse_args()


def clean_eggnog(filein, incl):
# 0 - 3     query_name, seed_eggNOG_ortholog, seed_ortholog_evalue, seed_ortholog_score,
# 4- 7      eggNOG_OGs, max_annot_lvl, COG_cat, Description
# 8 - 12    Preferred_name, GOs, EC, KEGG_ko, KEGG_Pathway,
# 13 - 16   KEGG_Module, KEGG_Reaction, KEGG_rclass, BRITE,
# 17 - 20     KEGG_TC, CAZy, BiGG_Reaction, PFAMs

    with open(filein, "r") as flin:
        for line in flin:
            line = line.rstrip()
            if not line.startswith("#"):
                #eggNOG_OGs\tCOG_cat\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tBiGG_Reaction\tDescription
                list_line=line.split("\t")
                id=list_line[0]
                select_items=[4,6,8,9,10,11,12,13,14,15,16,17,19,7]
                sel_descp=[list_line[i] for i in select_items]
                incl[id]=sel_descp

    return incl

def get_dbcan(file, dic, includ):
    with open(file, "r") as fin:
        '''
        #                                                                 --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
        # target name        accession  query name             accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
        #------------------- ----------   -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
        GT81.hmm             -          co_assembly_k127_985_2 -            6.6e-17   55.8   0.0   8.4e-17   55.4   0.0   1.1   1   0   0   1   1   1   1 -
        GT2.hmm              -          co_assembly_k127_985_2 -            1.9e-16   55.0   0.1   2.5e-16   54.6   0.1   1.2   1   0   0   1   1   1   1 -
        '''
        for line in fin:
            line=line.rstrip()
            if not line.startswith("#"):
                query=line.split()[0].rstrip() #2
                score=line.split()[5]
                if query in dic:
                    check=dic[query].split()[5]
                    if check > score:
                        dic[query]=line
                        includ[query]=line.split()[2].split(".")[0] #0
                else:
                    dic[query]=line
                    includ[query]=line.split()[2].split(".")[0]

    return dic, includ

def get_pfam(file, dic, include):
    with open(file, "r") as fin:
        '''
        #        	            	           	                     	           	--- full 	sequence	 ---- 	--- best 	1	 domain	 ---- 	--- domain	 number 	estimation 	----
        # target name        	accession  	query name           	accession  	  E-value	  score 	 bias 	  E-value	 	 score 	 bias 	  exp reg 	clu  ov 	env dom rep	 inc 	description of 	target
        #------------------- 	---------- 	-------------------- 	---------- 	---------	 ------ 	----- 	---------	 	------ 	----- 	  --- --- 	--- --- 	--- --- ---	 --- 	---------------	------
        TPR_19   	            PF14559.10 	LAVU02000112.1_1     	-          	1×10^-18	   67.6 	  1.2 	  7.6e-07	 	  29.6 	  0.0 	  3.4   2 	  1   1 	  3   3   1	1	Tetratricopepti	de repeat
        '''
        for line in fin:
            line=line.rstrip()
            if not line.startswith("#"):
                target=line.split()[0].rstrip() #2
                score=line.split()[5]

                if target in dic:
                    check=dic[target].split()[5]
                    if check > score:
                        dic[target]=line
                        include[target]=line.split()[1:3] #:2
                else:
                    dic[target]=line
                    include[target]=line.split()[1:3] #:2

    return dic, include

def print_annotation(fileout, dic):
    with open(fileout, "w") as fout:
        print("#                                                                 --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----", file=fout)
        print("# target name        accession  query name             accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target", file=fout)
        print("#------------------- ----------   -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------", file=fout)
        for l in dic.values():
            print(l, file=fout)



dbcan_dict={}
incl_dbcan={}
dbcan_dict, incl_dbcan=get_dbcan(args.d, dbcan_dict, incl_dbcan)

#print_annotation(args.o, dbcan_dict)

print("Total genes with dbcan annotation: {}".format(len(dbcan_dict.keys())))

pfam_dict={}
incl_pfam={}
pfam_dict, incl_pfam=get_pfam(args.p, pfam_dict, incl_pfam)

print("Total genes with pfam annotation: {}".format(len(pfam_dict.keys())))

incl_egg={}
incl_egg=clean_eggnog(args.e, incl_egg)


print("Total genes with EggNOG annotation: {}".format(len(incl_egg.keys())))

genes=set(incl_pfam.keys())
genes.update(set(incl_dbcan.keys()))
genes.update(set(incl_egg.keys()))

print("Total genes with at least one annotation: {}".format(len(genes)))

with open(args.o, "w") as fout:
#    print("#GeneID\tdbCAN_family\tPFAM_family\tPFAM_accession\tEggNOG_annotation", file=fout)
    print("#GeneID\tdbCAN_family\tPFAM_family\tPFAM_accession\teggNOG_OGs\tCOG_cat\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tBiGG_Reaction\tDescription", file=fout)
    for gene in genes:
            texto=str(gene)
            if gene in incl_dbcan:
                    texto+="\t"+str(incl_dbcan[gene])
            else:
                    texto+="\t-"
            if gene in incl_pfam:
                    texto+="\t"+str("\t".join(incl_pfam[gene]))
            else:
                    texto+="\t-"*2  #len(inc_pfam[gene])
            if gene in incl_egg:
                    texto+="\t"+str("\t".join(incl_egg[gene]))
            else:
                    texto+="\t-"*14 #len(incl_egg[gene])

            print(texto, file=fout)
