import os

configfile: "./mix_gc_config.json"

workdir: config["work_dir"]

def contigs_ind_assemblies(input_file):
  ind_samples = []
  ind_files = []
  with open(input_file) as fin:
    for line in fin:
      line=line.rstrip()
      if not line.startswith("#"):
        line=line.split(",")
        ind_samples.append(line[0])
        ind_files.append(line[1])
  return(ind_samples, ind_files)


names, path_files = contigs_ind_assemblies(config["path_to_ind_assembly_list"])


rule all:
    input:  "Gene_catalog/rep_proteins.faa", "Gene_catalog/rep_genes.fna",
            "Gene_catalog/rep_annotation.tsv", "Gene_catalog/rep_contigs.fasta.gz",
            "Gene_catalog/rep_genes_taxonomy.tsv", "Gene_catalog/rep_clusters.tsv",
            "Gene_catalog/rep_genes_taxonomy_krona.html"

rule unzip_and_check_ind_contigs_files:
    input:  config["path_to_ind_assembly_list"]
    output: temp(expand("Ind_assembly_contigs_temp/{s}.contigs.fna", s=names))
    message: "Reading Ind-assembly contigs"
    shell: """
              cat {input} | while read line
                do
                    s=$(echo $line | cut -d "," -f1)
                    f=$(echo $line | cut -d "," -f2)
                    if [[ -s "$f" ]]; then
                        if [[ "$f" =~ \.gz$ ]]; then
                            gunzip -cd $f | sed s/">"/">$s::"/ > Ind_assembly_contigs_temp/$s.contigs.fna
                        else
                            cat $f | sed s/">"/">$s::"/ > Ind_assembly_contigs_temp/$s.contigs.fna
                        fi
                    else
                        echo "ERROR: $f doesn't exists or is empty"
                        exit 1
                    fi
                done
           """

rule ind_proteins:
    input:  "Ind_assembly_contigs_temp/{s}.contigs.fna"
    output: p="Ind_assembly_dir/{s}.faa.gz",
            g="Ind_assembly_dir/{s}.fna.gz",
            f="Ind_assembly_dir/{s}.gff.gz"
    message: "Gene calling on Ind-assembly contigs"
    params: pm=config["prodigal_params"],
            p="Ind_assembly_dir/{s}.faa",
            g="Ind_assembly_dir/{s}.fna",
            f="Ind_assembly_dir/{s}.gff"
    conda: "mix_gc_env"
    threads: 1
    shell:  """
                prodigal -i {input} -a {params.p} -d {params.g} -f gff -o {params.f} {params.pm}
                gzip {params.p} {params.g} {params.f}
            """

rule co_in_chunks:
    input: config["path_to_co_assembly_contigs"]
    output: temp(dynamic("Co_assembly_dir_temp/{cunk}_co_assembly_contig.fasta"))
    message: "Co-assembly contigs in chuncks"
    params: int(config["threads"])
    threads: 1
    shell:  """
                if [[ {input} =~ \.gz$ ]]; then
                    n=$(gunzip -cd {input} | grep -c ">")
                else
                    n=$(grep -c ">" {input})
                fi
                python src/split_co_assembly.py -i {input} -c {params} -n $n -o Co_assembly_dir_temp
            """


rule co_proteins_chunck:
    input: "Co_assembly_dir_temp/{cunk}_co_assembly_contig.fasta"
    output: p=temp("Co_assembly_dir_temp/{cunk}_co_assembly.faa.gz"),
            g=temp("Co_assembly_dir_temp/{cunk}_co_assembly.fna.gz"),
            f=temp("Co_assembly_dir_temp/{cunk}_co_assembly.gff.gz")
    message: "Gene calling on Co-assembly contigs chunk"
    params: pm=config["prodigal_params"], t="Co_assembly_dir_temp/{cunk}_input_prodigal_tmp.fasta",
            p=temp("Co_assembly_dir_temp/{cunk}_co_assembly.faa"),
            g=temp("Co_assembly_dir_temp/{cunk}_co_assembly.fna"),
            f=temp("Co_assembly_dir_temp/{cunk}_co_assembly.gff")
    conda: "mix_gc_env"
    threads: 1
    shell:  """
               cat {input} | sed s/">"/">CO::"/ > {params.t}
               prodigal -i {params.t} -a {params.p} -d {params.g} -f gff -o {params.f} {params.pm}
               gzip {params.p} {params.g} {params.f}
               rm {params.t}
            """

rule merge_gene_calling_co_assembly:
    input:  p=dynamic("Co_assembly_dir_temp/{cunk}_co_assembly.faa.gz"),
            g=dynamic("Co_assembly_dir_temp/{cunk}_co_assembly.fna.gz"),
            f=dynamic("Co_assembly_dir_temp/{cunk}_co_assembly.gff.gz")
    output: p="Co_assembly_dir/co_assembly.faa.gz",
            g="Co_assembly_dir/co_assembly.fna.gz",
            f="Co_assembly_dir/co_assembly.gff.gz"
    message: "Concatenating gene calling co-assembly files"
    threads: 1
    shell:  """
               cat {input.p} > {output.p}
               cat {input.g} > {output.g}
               cat {input.f} > {output.f}
            """


rule cluster_ind:
    input: expand("Ind_assembly_dir/{s}.faa.gz", s=names)
    output: "Cluster/ind_assembly_rep_seq.fasta.gz"
    message: "Clustering Ind-assembly proteins"
    params: tmp=config["tmp_dir"]+"/ind",
            o="Cluster/ind_assembly",
            uz="Cluster/ind_assembly_rep_seq.fasta",
            pm=config["cluster_params"]
    conda: "mix_gc_env"
    threads: config["threads"]
    shell:  """
              mkdir -p {params.tmp}
              mmseqs easy-cluster {input} {params.o} {params.tmp} {params.pm} --threads {threads} -v 0
              gzip {params.uz}
            """


rule cluster_mix:
    input:  "Cluster/ind_assembly_rep_seq.fasta.gz","Co_assembly_dir/co_assembly.faa.gz"
    output: "Cluster/mix_rep_seq.fasta", "Cluster/mix_cluster.tsv"
    message: "Clustering Co-assembly and representative Ind-assembly proteins"
    params: tmp=config["tmp_dir"]+"/mix",
            o="Cluster/mix",
            pm=config["cluster_params"]
    conda: "mix_gc_env"
    threads: config["threads"]
    shell:  """
              mkdir -p {params.tmp}
              mmseqs easy-cluster {input} {params.o} {params.tmp} {params.pm} --threads {threads} -v 0
            """

rule ind_rep_genes:
  input: i=expand("Ind_assembly_dir/{s}.fna.gz", s=names), r="Cluster/ind_assembly_rep_seq.fasta.gz"
  output: "Cluster/ind_assembly_rep_seq.fna.gz"
  message: "Printing out Ind-assembly representative genes"
  params: i="Ind_assembly_dir"
  shell:  """
              python src/from_faa_gzip_to_fna.py -id {params.i} -r {input.r} -o Cluster/ind_assembly_rep_seq.fna
              gzip Cluster/ind_assembly_rep_seq.fna
          """

rule Mix_genes:
    input:  i="Cluster/ind_assembly_rep_seq.fna.gz",
            c="Co_assembly_dir/co_assembly.fna.gz",
            m="Cluster/mix_rep_seq.fasta"
    output: "Cluster/mix_rep_genes.fna"
    message: "Printing out Mix-assembly representative genes"
    threads: 1
    shell:  """
              python src/from_faa_gzip_to_fna_mix_assembly.py -i {input.i} -c {input.c} -m {input.m} -o {output}
            """

rule rfam_sto:
  output: "Rfam_db/Rfam.hmm"
  params: "Rfam_db/RFAM.sto"
  threads: config["threads"]
  message: "Downloading Rfam database"
  conda: "mix_gc_env"
  shell:  """
            if [[ ! -s Rfam_db/Rfam.seed ]]; then
              wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.seed.gz -P Rfam_db/
              gunzip Rfam_db/Rfam.seed.gz
            fi
            grep -v '#=GC' Rfam_db/Rfam.seed | grep -v 'Binary file Rfam.seed matches' > {params}
            hmmbuild  --dna --cpu {threads} {output} {params} > /dev/null
          """

rule rfam:
    input: g="Cluster/mix_rep_genes.fna", r="Rfam_db/Rfam.hmm"
    output: "annotation/rfam/mix_RFAM.tblout"
    threads: config["threads"]
    conda: "mix_gc_env"
    message: "rfam on Mix-assembly representative genes"
    shell: """
                if [[ ! -s {input.r}.h3i ]]; then
                    hmmpress {input.r}
                fi
                    nhmmscan --noali --cut_ga --cpu {threads} --tblout {output} {input.r} {input.g} > /dev/null
            """

rule remove_false_positive:
  input: r="annotation/rfam/mix_RFAM.tblout",
         p="Cluster/mix_rep_seq.fasta",
         g="Cluster/mix_rep_genes.fna"
  output: p="Gene_catalog/rep_proteins.faa",
          g="Gene_catalog/rep_genes.fna"
  params: m=config["list_specific_genes_to_remove"]
  message: "Removing false positive predicted genes"
  threads: 1
  shell:  """
                if [ ! -s list_rRNA_genes_in_mix.txt ]; then
                    grep 'rRNA' {input.r} | awk '{{print $3}}' > list_rRNA_genes_in_mix.txt
                fi
                python src/clean_mix_assembly_genes.py -g {input.g} -p {input.p} -m {params.m} -o {output.g} -u {output.p}
          """

rule pfam:
    input: "Gene_catalog/rep_proteins.faa"
    output: "annotation/pfam/mix_assembly_pfam.tsv"
    params: i=config["path_pfam_db"]
    threads: config["threads"]
    conda: "mix_gc_env"
    message: "Pfam annotation"
    shell: """
                if [[ ! -s {params.i} ]]; then
                    wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
                    gunzip Pfam-A.hmm.gz
                    mv Pfam-A.hmm {params.i}
                fi
                if [[ ! -s {params.i}.h3i ]]; then
                    hmmpress {params.i}
                fi
                hmmscan --noali --cut_ga --cpu {threads} --tblout {output} {params.i} {input} > /dev/null
           """

rule dbcan:
    input: "Gene_catalog/rep_proteins.faa"
    output: "annotation/dbcan/mix_assembly_dbcan.tsv"
    params: i=config["path_dbcan_db"]
    threads: config["threads"]
    conda: "mix_gc_env"
    message: "dbcan annotation"
    shell:  """
                if [[ ! -s {params.i} ]]; then
                        wget https://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-old@UGA/dbCAN-fam-HMMs.txt
                        mv dbCAN-fam-HMMs.txt {params.i}
                fi
                if [[ ! -s {params.i}.h3i ]]; then
                    hmmpress {params.i}
                fi
                hmmscan --noali --cpu {threads} --tblout {output} {params.i} {input} > /dev/null
            """

rule eggnog:
    input: "Gene_catalog/rep_proteins.faa"
    output: "annotation/eggnog/rep_proteins.emapper.annotations"
    params: db=config["path_to_eggnog_db"]
    conda: "eggnog_mapper_env"
    threads: config["threads"]
    message: "Eggnog annotation"
    shell:  """
                filehit=annotation/eggnog/rep_proteins.emapper.hits
                if [[ -s $filehit ]]; then
                    echo "INFO: $filehit is present, annotation will be resumed"
                    emapper.py -i {input} -o rep_proteins --output_dir annotation/eggnog --cpu {threads} --data_dir {params.db} -m diamond --resume
                else
                    emapper.py -i {input} -o rep_proteins --output_dir annotation/eggnog --cpu {threads} --data_dir {params.db} -m diamond
                fi
            """

rule annotation_table:
    input:  d="annotation/dbcan/mix_assembly_dbcan.tsv", p="annotation/pfam/mix_assembly_pfam.tsv",
            e="annotation/eggnog/rep_proteins.emapper.annotations"
    output: "Gene_catalog/rep_annotation.tsv"
    message: "Printing out Mix_assembly representative genes annotations"
    threads: 1
    shell: "python src/merge_mix_assembly_genes_annotation_tables.py -d {input.d} -p {input.p} -e {input.e} -o {output}"

rule mix_contigs:
    input:  m="Gene_catalog/rep_proteins.faa",
            i=path_files,
            c=config["path_to_co_assembly_contigs"]
    output: "Gene_catalog/rep_contigs.fasta.gz"
    message: "Printing out Mix_assembly contigs"
    params: l=config["path_to_ind_assembly_list"],
            o="Gene_catalog/rep_contigs.fasta"
    shell:  """
              python src/from_fna_gzip_to_contigs_mix_assembly_fast.py -l {params.l} -c {input.c} -m {input.m} -o {params.o}
              gzip {params.o}
            """

#taxonomy

if config["taxonomy_DB"] == "gtdb":

    rule download_gtdb_ref_proteins:
      output: d=directory(config["GTDB"]["GTDB_dir"]), b=config["GTDB"]["bact_tsv"], a=config["GTDB"]["arch_tsv"]
      message: "Downloading GTDB proteins database"
      shell:  """
                    mkdir -p GTDB_aa_db
                    cd GTDB_aa_db/
                    wget https://data.gtdb.ecogenomic.org/releases/latest/genomic_files_reps/gtdb_proteins_aa_reps.tar.gz
                    tar -xf gtdb_proteins_aa_reps.tar.gz
                    wget -A _taxonomy.tsv -r -l 1 -nd https://data.gtdb.ecogenomic.org/releases/latest/
                    cd ..
              """

    rule DB_fasta:
      input: d=config["GTDB"]["GTDB_dir"], b=config["GTDB"]["bact_tsv"], a=config["GTDB"]["arch_tsv"]
      output: "mmseqs2DBs/GTDB_mmseqs2_format.fasta"
      message: "Creating GTDB mmseqs2DBs (step 1/4)"
      shell:"python src/create_GTDB_fasta_db.py -d {input.d} -b {input.b} -a {input.a} -o {output}"

    rule transformation:
      input: "mmseqs2DBs/GTDB_mmseqs2_format.fasta"
      output: expand("{temp}/taxonomy/{files}.dmp", temp=config["tmp_dir"], files=["delnodes","merged", "names", "nodes"]), config["tmp_dir"]+"/taxonomy/mapping"
      params: w=config["work_dir"], t=config["tmp_dir"]
      message: "Creating GTDB mmseqs2DBs seq_mapping and seq_taxonomy (step 2/4)"
      shell: "bash src/build_name_dmp_node_dmp.sh {input} {params.w} {params.t}"

    rule createdb:
      input: i="mmseqs2DBs/GTDB_mmseqs2_format.fasta",
             a=expand("{temp}/taxonomy/{files}.dmp", temp=config["tmp_dir"], files=["delnodes","merged", "names", "nodes"]),
             b=config["tmp_dir"]+"/taxonomy/mapping"
      output: db="mmseqs2DBs/GTDB/SeqDB", h="mmseqs2DBs/GTDB/SeqDB_h"
      message: "Creating GTDB mmseqs2DBs (step 3/4)"
      conda: "mix_gc_env"
      shell: "mmseqs createdb {input.i} {output.db}"

    rule taxdb:
      input:  "mmseqs2DBs/GTDB/SeqDB"
      output: "mmseqs2DBs/GTDB/SeqDB_mapping", b="mmseqs2DBs/GTDB/SeqDB_taxonomy"
      params: t1=config["tmp_dir"]+"/taxo", t2=config["tmp_dir"]+"/taxonomy", t3=config["tmp_dir"]+"/taxonomy/mapping"
      message: "Creating GTDB mmseqs2DBs (step 4/4)"
      conda: "mix_gc_env"
      shell:  """
                mkdir -p {params.t1}
                mmseqs createtaxdb {input} {params.t1} --ncbi-tax-dump {params.t2} --tax-mapping-file {params.t3}
              """

    rule taxonomy_gtdb:
      input:  c="Gene_catalog/rep_contigs.fasta.gz",
              db="mmseqs2DBs/GTDB/SeqDB",
              a="mmseqs2DBs/GTDB/SeqDB_mapping",
              b="mmseqs2DBs/GTDB/SeqDB_taxonomy"
      output: "Taxonomy_gtdb/Mix_easy_tax_lca.tsv"
      message: "Representative genes taxonomy assignment - GTDB"
      conda: "mix_gc_env"
      threads: config["threads"]
      params: p=config["mmseqs_taxonomy_params"], t=config["tmp_dir"]+"/taxo_gtdb", o="Taxonomy_gtdb/Mix_easy_tax"
      shell:  """
                mkdir -p {params.t}
                mmseqs easy-taxonomy {input.c} {input.db} {params.o} {params.t} {params.p} --threads {threads}
              """


rule download_create_uniprotDB:
  output: db="mmseqs2DBs/Uniprot/SeqDB",
          h="mmseqs2DBs/Uniprot/SeqDB_h",
          a="mmseqs2DBs/Uniprot/SeqDB_mapping",
          b="mmseqs2DBs/Uniprot/SeqDB_taxonomy"
  threads: config["threads"]
  params: t=config["tmp_dir"]+"/taxo", type=config["uniprot"]["uniprot_db_type"]
  message: "Downloading Uniprot database"
  conda: "mix_gc_env"
  shell:  """
            mkdir -p {params.t}
            mmseqs databases {params.type} {output.db} {params.t} --threads {threads}
          """

if config["taxonomy_DB"] == "uniprot" and config["uniprot"]["by_chuncks"] == False :

   rule taxonomy_uniprot:
        input:  c="Gene_catalog/rep_contigs.fasta.gz",
                db="mmseqs2DBs/Uniprot/SeqDB",
                a="mmseqs2DBs/Uniprot/SeqDB_mapping",
                b="mmseqs2DBs/Uniprot/SeqDB_taxonomy"
        output: "Taxonomy/Mix_easy_tax_lca.tsv", r="Taxonomy/Mix_easy_tax_report"
        threads: config["threads"]
        message: "Representative genes taxonomy assignment - Uniprot"
        conda: "mix_gc_env"
        params: p=config["mmseqs_taxonomy_params"], t=config["tmp_dir"]+"/taxo", o="Taxonomy/Mix_easy_tax"
        shell:  """
                        mkdir -p {params.t}
                        mmseqs easy-taxonomy {input.c} {input.db} {params.o} {params.t} {params.p} --threads {threads}
                """

   rule taxonomy_uniport_tsv:
        input: "Taxonomy/Mix_easy_tax_lca.tsv"
        output: "Gene_catalog/rep_genes_taxonomy.tsv"
        threads: 1
        message: "Parsing taxonomy assignments"
        shell: "python src/parse_results_tsv.py -i {input} -o {output}"

if config["taxonomy_DB"] == "gtdb" or config["taxonomy_DB"] == "uniprot" and config["uniprot"]["by_chuncks"] == True :

  rule filtre_Virus:
    input:  db="mmseqs2DBs/Uniprot/SeqDB",
            h="mmseqs2DBs/Uniprot/SeqDB_h",
            a="mmseqs2DBs/Uniprot/SeqDB_mapping",
            b="mmseqs2DBs/Uniprot/SeqDB_taxonomy"
    output: db="mmseqs2DBs/Uniprot_Virus/SeqDB",
            h="mmseqs2DBs/Uniprot_Virus/SeqDB_h",
            a="mmseqs2DBs/Uniprot_Virus/SeqDB_mapping",
            b="mmseqs2DBs/Uniprot_Virus/SeqDB_taxonomy"
    threads: config["threads"]
    message: "Filtering Uniprot database - Virus"
    conda: "mix_gc_env"
    params: p="--taxon-list 10239",
            w=config["work_dir"],
            s=config["work_dir"]+"/mmseqs2DBs/",
            ai="Uniprot/SeqDB_mapping", ao="Uniprot_Virus/SeqDB_mapping",
            bi="Uniprot/SeqDB_taxonomy", bo="Uniprot_Virus/SeqDB_taxonomy"
    shell: '''
              mmseqs filtertaxseqdb {input.db} {output.db} {params.p} --threads {threads} -v 0
              cd {params.s}
              cp {params.ai} tempo_map & mv tempo_map {params.ao}
              cp {params.bi} tempo_tax & mv tempo_tax {params.bo}
              cd {params.w}
            '''

  rule taxonomy_virus:
    input:  c="Gene_catalog/rep_contigs.fasta.gz",
            db="mmseqs2DBs/Uniprot_Virus/SeqDB",
            a="mmseqs2DBs/Uniprot_Virus/SeqDB_mapping",
            b="mmseqs2DBs/Uniprot_Virus/SeqDB_taxonomy"
    output: "Taxonomy_Virus/Mix_easy_tax_lca.tsv"
    threads: config["threads"]
    message: "Taxonomy assignment - Virus"
    conda: "mix_gc_env"
    params: p=config["mmseqs_taxonomy_params"], t=config["tmp_dir"]+"/taxo_virus", o="Taxonomy_Virus/Mix_easy_tax"
    shell:  """
              mkdir -p {params.t}
              mmseqs easy-taxonomy {input.c} {input.db} {params.o} {params.t} {params.p} --threads {threads}
            """

  rule filtre_Euka:
    input:  db="mmseqs2DBs/Uniprot/SeqDB",
            h="mmseqs2DBs/Uniprot/SeqDB_h",
            a="mmseqs2DBs/Uniprot/SeqDB_mapping",
            b="mmseqs2DBs/Uniprot/SeqDB_taxonomy"
    output: db="mmseqs2DBs/Uniprot_Euka/SeqDB",
            h="mmseqs2DBs/Uniprot_Euka/SeqDB_h",
            a="mmseqs2DBs/Uniprot_Euka/SeqDB_mapping",
            b="mmseqs2DBs/Uniprot_Euka/SeqDB_taxonomy"
    threads: config["threads"]
    message: "Filtering Uniprot database - Eukaryotes"
    conda: "mix_gc_env"
    params: p="--taxon-list 2759",
            w=config["work_dir"],
            s=config["work_dir"]+"/mmseqs2DBs/",
            ai="Uniprot/SeqDB_mapping", ao="Uniprot_Euka/SeqDB_mapping",
            bi="Uniprot/SeqDB_taxonomy", bo="Uniprot_Euka/SeqDB_taxonomy"
    shell: '''
              mmseqs filtertaxseqdb {input.db} {output.db} {params.p} --threads {threads} -v 0
              cd {params.s}
              cp {params.ai} tempo_map & mv tempo_map {params.ao}
              cp {params.bi} tempo_tax & mv tempo_tax {params.bo}
              cd {params.w}
            '''


  rule taxonomy_Euka:
    input:  c="Gene_catalog/rep_contigs.fasta.gz",
            db="mmseqs2DBs/Uniprot_Euka/SeqDB",
            a="mmseqs2DBs/Uniprot_Euka/SeqDB_mapping",
            b="mmseqs2DBs/Uniprot_Euka/SeqDB_taxonomy"
    output: "Taxonomy_Euka/Mix_easy_tax_lca.tsv"
    threads: config["threads"]
    message: "Taxonomy assignment - Eukaryotes"
    conda: "mix_gc_env"
    params: p=config["mmseqs_taxonomy_params"], t=config["tmp_dir"]+"/taxo_euka", o="Taxonomy_Euka/Mix_easy_tax"
    shell:  """
              mkdir -p {params.t}
              mmseqs easy-taxonomy {input.c} {input.db} {params.o} {params.t} {params.p} --threads {threads}
            """

if config["taxonomy_DB"] == "uniprot" and config["uniprot"]["by_chuncks"] == True :
  rule filtre_arch:
    input:  db="mmseqs2DBs/Uniprot/SeqDB",
            h="mmseqs2DBs/Uniprot/SeqDB_h",
            a="mmseqs2DBs/Uniprot/SeqDB_mapping",
            b="mmseqs2DBs/Uniprot/SeqDB_taxonomy"
    output: db="mmseqs2DBs/Uniprot_arch/SeqDB",
            h="mmseqs2DBs/Uniprot_arch/SeqDB_h",
            a="mmseqs2DBs/Uniprot_arch/SeqDB_mapping",
            b="mmseqs2DBs/Uniprot_arch/SeqDB_taxonomy"
    threads: config["threads"]
    message: "Filtering Uniprot database - Archeae"
    conda: "mix_gc_env"
    params: p="--taxon-list 2157",
            w=config["work_dir"],
            s=config["work_dir"]+"/mmseqs2DBs/",
            ai="Uniprot/SeqDB_mapping", ao="Uniprot_arch/SeqDB_mapping",
            bi="Uniprot/SeqDB_taxonomy", bo="Uniprot_arch/SeqDB_taxonomy"
    shell: '''
              mmseqs filtertaxseqdb {input.db} {output.db} {params.p} --threads {threads} -v 0
              cd {params.s}
              cp {params.ai} tempo_map & mv tempo_map {params.ao}
              cp {params.bi} tempo_tax & mv tempo_tax {params.bo}
              cd {params.w}
            '''

  rule taxonomy_arch:
    input:  c="Gene_catalog/rep_contigs.fasta.gz",
            db="mmseqs2DBs/Uniprot_arch/SeqDB",
            a="mmseqs2DBs/Uniprot_arch/SeqDB_mapping",
            b="mmseqs2DBs/Uniprot_arch/SeqDB_taxonomy"
    output: "Taxonomy_arch/Mix_easy_tax_lca.tsv"
    message: "Taxonomy assignment - Archeae"
    threads: config["threads"]
    conda: "mix_gc_env"
    params: p=config["mmseqs_taxonomy_params"], t=config["tmp_dir"]+"/taxo_arch", o="Taxonomy_arch/Mix_easy_tax"
    shell:  """
              mkdir -p {params.t}
              mmseqs easy-taxonomy {input.c} {input.db} {params.o} {params.t} {params.p} --threads {threads}
            """

  rule filtre_bact:
    input:  db="mmseqs2DBs/Uniprot/SeqDB",
            h="mmseqs2DBs/Uniprot/SeqDB_h",
            a="mmseqs2DBs/Uniprot/SeqDB_mapping",
            b="mmseqs2DBs/Uniprot/SeqDB_taxonomy"
    output: db="mmseqs2DBs/Uniprot_Bact/SeqDB",
            h="mmseqs2DBs/Uniprot_Bact/SeqDB_h",
            a="mmseqs2DBs/Uniprot_Bact/SeqDB_mapping",
            b="mmseqs2DBs/Uniprot_Bact/SeqDB_taxonomy"
    threads: config["threads"]
    message: "Filtering Uniprot database - Bacteria"
    conda: "mix_gc_env"
    params: p="--taxon-list 2",
            w=config["work_dir"],
            s=config["work_dir"]+"/mmseqs2DBs/",
            ai="Uniprot/SeqDB_mapping", ao="Uniprot_Bact/SeqDB_mapping",
            bi="Uniprot/SeqDB_taxonomy", bo="Uniprot_Bact/SeqDB_taxonomy"
    shell: '''
              mmseqs filtertaxseqdb {input.db} {output.db} {params.p} --threads {threads}
              cd {params.s}
              cp {params.ai} tempo_map & mv tempo_map {params.ao}
              cp {params.bi} tempo_tax & mv tempo_tax {params.bo}
              cd {params.w}
            '''

  rule taxonomy_Bact:  
    input:  c="Gene_catalog/rep_contigs.fasta.gz",
            db="mmseqs2DBs/Uniprot_Bact/SeqDB",
            a="mmseqs2DBs/Uniprot_Bact/SeqDB_mapping",
            b="mmseqs2DBs/Uniprot_Bact/SeqDB_taxonomy"
    output: "Taxonomy_Bact/Mix_easy_tax_lca.tsv"
    message: "Taxonomy assignment - Bacteria"
    threads: config["threads"]
    conda: "mix_gc_env"
    params: p=config["mmseqs_taxonomy_params"], t=config["tmp_dir"]+"/taxo_Bact", o="Taxonomy_Bact/Mix_easy_tax"
    shell:  """
              mkdir -p {params.t}
              mmseqs easy-taxonomy {input.c} {input.db} {params.o} {params.t} {params.p} --threads {threads}
            """


if config["taxonomy_DB"] == "gtdb":
  rule taxonomy_combined:
    input:  "Taxonomy_gtdb/Mix_easy_tax_lca.tsv",
            "Taxonomy_Virus/Mix_easy_tax_lca.tsv",
            "Taxonomy_Euka/Mix_easy_tax_lca.tsv"
    output: "Gene_catalog/rep_genes_taxonomy.tsv"
    params: i="Taxonomy_gtdb/Mix_easy_tax_lca.tsv,Taxonomy_Virus/Mix_easy_tax_lca.tsv,Taxonomy_Euka/Mix_easy_tax_lca.tsv"
    threads: 1
    message: "Merging Taxonomy assignment - GTDB/Uniprot"
    shell: "python src/parse_results_tsv.py -i {params.i} -o {output}"

elif config["taxonomy_DB"] == "uniprot" and config["uniprot"]["by_chuncks"] == True :
  rule taxonomy_uniprot_combined:
    input:  "Taxonomy_Bact/Mix_easy_tax_lca.tsv",
            "Taxonomy_arch/Mix_easy_tax_lca.tsv",
            "Taxonomy_Virus/Mix_easy_tax_lca.tsv",
            "Taxonomy_Euka/Mix_easy_tax_lca.tsv"
    output: "Gene_catalog/rep_genes_taxonomy.tsv"
    params: i="Taxonomy_Bact/Mix_easy_tax_lca.tsv,Taxonomy_arch/Mix_easy_tax_lca.tsv,Taxonomy_Virus/Mix_easy_tax_lca.tsv,Taxonomy_Euka/Mix_easy_tax_lca.tsv"
    threads: 1
    message: "Merging Taxonomy assignment - Uniprot"
    shell: "python src/parse_results_tsv.py -i {params.i} -o {output}"


rule krona:
   input: "Gene_catalog/rep_genes_taxonomy.tsv"
   output: "Gene_catalog/rep_genes_taxonomy_krona.txt"
   shell: "python src/krona_clean_results_tsv.py -i {input} -o {output}"

rule html:
   input: "Gene_catalog/rep_genes_taxonomy_krona.txt"
   output: "Gene_catalog/rep_genes_taxonomy_krona.html"
   conda: "mix_gc_env"
   shell: "ktImportText {input} -o {output}"

rule table_rep_clusters:
    input: i="Cluster/mix_cluster.tsv", r="Gene_catalog/rep_proteins.faa"
    output:"Gene_catalog/rep_clusters.tsv"
    threads: 1
    shell: "python src/mix_cluster_formated.py -i {input.i} -r {input.r} -o {output}"
