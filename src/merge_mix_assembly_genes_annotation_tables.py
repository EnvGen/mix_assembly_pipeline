#!/usr/bin/env python3

import os
import argparse


usage = 'python scr/merge_mix_assembly_genes_annotation_tables.py -d -e -t'
description = 'This program creates a table resuming all annotation from representative genes'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument('-d', dest='d', help='file dbCAN annotation table',required=True)
parser.add_argument('-p', dest='p', help='file Pfam annotation table',required=True)
parser.add_argument('-r', dest='r', help='file Rfam annotation table',required=True)
parser.add_argument('-f', dest='f', help='gene names false positives rfam', default='list_rRNA_genes_in_mix.txt')
parser.add_argument('-e', dest='e', help='file EggNOG annotation', required=True)
parser.add_argument('-o', dest='o', help='output all annotation table',required=True)

args = parser.parse_args()


def clean_eggnog(filein, incl):
# 0 - 3     query_name, seed_eggNOG_ortholog, seed_ortholog_evalue, seed_ortholog_score,
# 4- 7      eggNOG_OGs, max_annot_lvl, COG_cat, Description
# 8 - 12    Preferred_name, GOs, EC, KEGG_ko, KEGG_Pathway,
# 13 - 16   KEGG_Module, KEGG_Reaction, KEGG_rclass, BRITE,
# 17 - 20     KEGG_TC, CAZy, BiGG_Reaction, PFAMs

#1 query_name
#2 seed_eggNOG_ortholog
#3 seed_ortholog_evalue
#4 seed_ortholog_score
#5 best_tax_level
#6 Preferred_name
#7 GOs
#8 EC
#9 KEGG_ko
#10 KEGG_Pathway
#11 KEGG_Module
#12 KEGG_Reaction
#13 KEGG_rclass
#14 BRITE
#15 KEGG_TC
#16 CAZy
#17 BiGG_Reaction
#18 taxonomic scope
#19 eggNOG OGs
#20 best eggNOG OG
#21 COG Functional cat.
#22 eggNOG free text desc.


    with open(filein, "r") as flin:
        for line in flin:
            line = line.rstrip()
            if not line.startswith("#"):
                #eggNOG_OGs\tCOG_cat\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tBiGG_Reaction\tDescription\tbest_tax_level\tCAZy
                list_line=line.split("\t")
                id=list_line[0]
                #select_items=[4,6,8,9,10,11,12,13,14,15,16,17,19,7]
                select_items=[18, 20, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 21, 4, 15]
                sel_descp=[list_line[i] if i < len(list_line) else "-" for i in select_items]
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


def get_false_positives(file):
    removable=set()
    with open(file, "r") as fin:
        for line in fin:
            line=line.rstrip()
            if "::" in line:
                removable.add(line)
    return removable



def get_rfam(file, dic, include, removable):


    '''
    # target name           accession  query name                    accession  hmmfrom hmm to alifrom  ali to envfrom  env to  modlen strand   E-value  score  bias  description of target
    #   ------------------- ----------          -------------------- ---------- ------- ------- ------- ------- ------- ------- ------- ------ --------- ------ ----- ---------------------
    #SSU_rRNA_bacteria      RF00177    P1994_119::P1994_119_k141_6488_2 -              748     866       1     113       1     114    1533    +     7.3e-25   82.1   3.1  Bacterial small subunit ribosomal RNA

    '''

    with open(file, "r") as fin:

        for line in fin:
            line=line.rstrip()
            if not line.startswith("#"):
                target=line.split()[2].rstrip() #2
                score=line.split()[13]
                if not target in removable:
                    if target in dic:
                        check=dic[target].split()[13]
                        if check > score:
                            dic[target]=line
                            include[target]=line.split()[:2] #:2
                    else:
                        dic[target]=line
                        include[target]=line.split()[:2] #:2

    return dic, include




def get_pfam(file, dic, include):
    with open(file, "r") as fin:
        '''
#                                                                               --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# target name                        accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
#                ------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
CO::co_assembly_k127_15632429_1      -          1-cysPrx_C           PF10417.13   1.1e-13   58.1   0.3   1.9e-13   57.3   0.3   1.4   1   0   0   1   1   1   1 # 69 # 593 # -1 # ID=739885_1;partial=00;start_type=TTG;rbs_motif=None;rbs_spacer=None;gc_cont=0.461
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
                        include[target]=line.split()[2:4] #:2
                else:
                    dic[target]=line
                    include[target]=line.split()[2:4] #:2

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

rfam_dict={}
incl_rfam={}
set_of_genes_to_remove=get_false_positives(args.f)
rfam_dict, incl_rfam=get_rfam(args.r, rfam_dict, incl_rfam, set_of_genes_to_remove)

print("Total genes with rfam annotation: {}".format(len(rfam_dict.keys())))

incl_egg={}
incl_egg=clean_eggnog(args.e, incl_egg)


print("Total genes with EggNOG annotation: {}".format(len(incl_egg.keys())))

genes=set(incl_pfam.keys())
genes.update(set(incl_dbcan.keys()))
genes.update(set(incl_egg.keys()))

print("Total genes with at least one annotation: {}".format(len(genes)))

with open(args.o, "w") as fout:
#    print("#GeneID\tdbCAN_family\tPFAM_family\tPFAM_accession\tEggNOG_annotation", file=fout)
    print("#GeneID\tdbCAN_family\tRFAM_description\tRFAM_accession\tPFAM_family\tPFAM_accession\teggNOG_OGs\tCOG_cat\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tBiGG_Reaction\tDescription\tbest_tax_level\tCAZy", file=fout)
    for gene in genes:
            texto=str(gene)
            if gene in incl_dbcan:
                    texto+="\t"+str(incl_dbcan[gene])
            else:
                    texto+="\t-"
            if gene in incl_rfam:
                    texto+="\t"+str("\t".join(incl_rfam[gene]))
            else:
                    texto+="\t-"*2  #len(inc_rfam[gene])
            if gene in incl_pfam:
                    texto+="\t"+str("\t".join(incl_pfam[gene]))
            else:
                    texto+="\t-"*2  #len(inc_pfam[gene])
            if gene in incl_egg:
                    texto+="\t"+str("\t".join(incl_egg[gene]))
            else:
                    texto+="\t-"*16 #len(incl_egg[gene])

            print(texto, file=fout)

