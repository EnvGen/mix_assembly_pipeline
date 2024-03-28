# Installation

	1. Clone the repository:

		git clone https://github.com/EnvGen/mix_assembly_pipeline

	2. Download Eggnog mapper and EggNOG database:

			conda create -n eggnog_mapper_env -c bioconda -c conda-forge eggnog-mapper=2.0.1 -y
			conda activate eggnog_mapper_env
			download_eggnog_data.py --data_dir <your option. Same path to be set in the mix_gc_config.json - "Path_to_eggnog_db" >
			conda deactivate

	3. If you want to use CAT (Contig Annotation Tool) with GTDB (Genome Taxonomy Database):

			conda create -n CAT_env -c bioconda -c conda-forge cat=5.3 -y

	4. Create the conda environment mix_gc_env:

		conda create -n mix_gc_env -c conda-forge -c bioconda prodigal=2.6.3 mmseqs2=13.45111 hmmer=3.3.2 krona=2.7 -y

	5. Install snakemake:

		conda create -n snakemake_env -c bioconda -c conda-forge snakemake=7.25.0 -y

# Usage

	1. Prepare your Ind_contigs_list.csv file: one line per sample following this format: <sample name>,<path to assembly file .fasta or fasta.gz>
 		Files can be compressed (only .gz). Make sure samples names nor contigs header use "::", e.g., >This_is_not_allowed::in_contigs_headers, No_::_in_sample_name


	2. Modify config file:


			nano mix_gc_config.json

			Directories:
			"work_dir": "/abs/path/to/mix_gc_pipeline",  #working directory
			"tmp_dir":"TMP_dir",  #Temporary directory

			input:
			"path_to_ind_assembly_list": "Ind_contigs_list.csv",
			"path_to_co_assembly_contigs": "/abs/path/to/co_assembly_contigs.fa.gz",

			Gene calling parameters:
			"prodigal_params": "-p meta -q",

			Clustering parameters:
			"cluster_params":"-c 0.95 --min-seq-id 0.95 --cov-mode 1 --cluster-mode 2",

			Functional annotation:
			"list_specific_genes_to_remove": "None", #in case you have a list of genes that you specifically want to remove. File with one gene name per line
			"path_to_eggnog_db": "abs/path/to/eggnog_db/data",
			"path_dbcan_db":"dbcan_db/dbCAN-fam-HMMs.txt", #It not available the pipeline will download it
			"path_pfam_db":"Pfam_db/Pfam-A.hmm", ##It not available the pipeline will download it

			Taxonomy assignment:

			Using CAT:
			"include_CAT_GTDB_taxonomy": True, # Using CAT contigs taxonomy with GTDB as reference database
			"path_to_CAT_gtdb_dir": "/abs/path/to/CAT_GTDB_database",

			Using Mmseqs2:
			"taxonomy_DB": "uniprot", #options gtdb, uniprot - The pipeline will download the databases. GTDB works only on prokaryotes, you need to include a Uniprot database here after for the taxonomy affiliation of other Kindoms.  
			"GTDB": {
			        "GTDB_dir":"GTDB_aa_db/protein_faa_reps",
			        "bact_tsv": "GTDB_aa_db/bac120_taxonomy.tsv",
			        "arch_tsv": "GTDB_aa_db/ar53_taxonomy.tsv"
			        },
			"uniprot": {
			            "uniprot_db_type":"UniRef90", #. Always provide an Uniprot database. See mmseqs2 manual for more Uniprot databases (https://mmseqs.com/latest/userguide.pdf)
			            "by_chuncks": True,  #Options: True, False. If True, it splits the taxonomy database into Kindoms to reduce the RAM requirements
			            },
			"mmseqs_taxonomy_params":"--tax-lineage 1 -v 1 --report-mode 1",
			"extra_mmseqs_taxonomy_params_Virus": "--orf-filter 0",

			To easy the functional annotation and taxonomy affilitaion steps, the process can be run by chunks
			"n_chuncks_annotations": 20 - Number of chuncks.

			"n_chuncks_taxonomy": 50


			Performance:

			"Big_jobs": {"Threads": 20,  #number of cpus to be used                   
			             "memory_per_cpu": 6400 # (MB)- Required when using HPC (high-performance computing) and several nodes
			                              },
			"Small_jobs": {"Threads": 1,
			               "memory_per_cpu": 6400 # (MB) - Required when using HPC (high-performance computing) and several nodes
			                                 },
			"slurm_extra":"", #Optional argument added to "big_jobs" when using HPC (high-performance computing) and several nodes. For instance, on Uppmax "-C fat" to get access to the fat nodes (12800 mb per cpu)


			Save modifications (CTRL +x, y).

	3. Run the pipeline:

			Using one node:

			conda activate snakemake_env
			snakemake -s mix_gc.smk --use-conda -j <number of cpus to be used>


			If you want to run the pipeline using several nodes in a HPC:

			a. Download cookiecutter:

			conda create -n cookiecutter_env -c conda-forge cookiecutter -y
			conda activate  cookiecutter_env

			b. Create the profile directory and answer the questions:

			profile_dir="/abs/path/to/mix_gc_pipeline/.config/snakemake"
			mkdir -p "$profile_dir"
			template="gh:Snakemake-Profiles/slurm"
			cookiecutter --output-dir "$profile_dir" "$template"
			conda deactivate

			c. Run the pipeline:

			conda activate snakemake_env
			snakemake --profile /abs/path/to/mix_gc_pipeline/.config/snakemake/<name_of_your_slurm_profile_file> -s mix_gc.smk


# Output
Gene_catalog is the main output folder, containing the following files:

 		Contigs and genes names(headers) will have the prefix sample name:: (if the gene comes from a individual assembly contig) or co:: (if the gene comes from a co_assembly contig)

		Gene catalog:
		        * rep_genes.fna.gz -- gene sequences in FASTA format    
		        * rep_proteins.faa.gz - - protein sequences in FASTA format  

		Functional annotations:					
		        * rep_annotations.tsv.gz - - This table regroups dbCAN, PFAM, RFAM (only the best hit, based on the highest score, is included) and EggNOG annotations per gene.

		Taxonomic annotations:
		        * Mmseqs2_rep_genes_taxonomy.tsv -- Taxonomy affiliation using mmseq2 and Uniref90
     			* CAT_rep_genes_taxonomy.tsv -- Taxonomy affiliation using CAT and GTDB

		Miscellaneous files:
			* Mmseqs2_rep_genes_taxonomy_krona.html - - Krona charts of representative gene taxonomic annotations (using mmseq2 and Uniref90)
		        * CAT_rep_genes_taxonomy_krona.html - - Krona charts of representative gene taxonomic annotations (using CAT and GTDB)   
			* rep_clusters_all.tsv - - Mix-assembly clusters with representative genes in the first column and cluster members in the second column (both ind and co-assembly)

	  		* rep_contigs.fasta.gz -- Contigs from which representative Mix-assembly genes were predicted.
	  		* Mmseqs2_rep_contigs_taxonomy_krona.html - - Krona charts of contigs (rep_contigs.fasta.gz) taxonomic annotations (using mmseq2 and Uniref90)
			* Mmseqs2_rep_contigs_taxonomy.tsv -- Tontigs (rep_contigs.fasta.gz taxonomy affiliation table using mmseq2 and Uniref90
