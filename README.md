# Installation

	1. Clone the repository:

		git clone https://github.com/lfdelzam/mix_gc_pipeline

	2. Download Eggnog mapper and EggNOG database:

			conda create -n eggnog_mapper_env -c bioconda -c conda-forge eggnog-mapper=2.0.1 -y
			conda activate eggnog_mapper_env
			download_eggnog_data.py --data_dir <your option. Same path to be set in the mix_gc_config.json - "Path_to_eggnog_db" >
			conda deactivate

	3. Create the conda environment mix_gc_env:

		conda create -n mix_gc_env -c conda-forge -c bioconda prodigal=2.6.3 mmseqs2=13.45111 hmmer=3.3.2 krona=2.7 -y

	4. Install snakemake:

		conda create -n snakemake_env -c bioconda -c conda-forge snakemake=7.25.0 -y

# Usage

	1. Prepare your Ind_contigs_list.csv file: one line per sample following this format: <sample name>,<path to assembly file .fasta or fasta.gz>


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
			"path_pfam_db":"/Users/luisdelgado/Documents/Mix_assembly_pipeline/Pfam_db/Pfam-A.hmm", ##It not available the pipeline will download it

			Taxonomy assignment:
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
			"mmseqs_taxonomy_params":"--lca-mode 4 --orf-filter 0 --tax-lineage 1 --filter-hits 1 -v 0 --min-seq-id 0.4 -s 1 --report-mode 1",

			Threads:
			"threads":20 - number of cpus to be used

			To easy the functional annotation steps, the process can be run by chunks
			"n_chuncks_annotations": 20 - Number of chuncks. If higher than number of threads, it will be set to max number of cpus available.



			Save modifications.

	3. Run the pipeline:

			conda activate snakemake_env
			snakemake -s mix_gc.smk --use-conda -j 20

	4.	Gene_catalog is the main output folder, containing the following files:


				rep_annotation.tsv
				rep_clusters.tsv
				rep_contigs_taxonomy_krona.html
				rep_contigs_taxonomy_krona.txt
				rep_contigs_taxonomy.tsv
				rep_contigs.fasta.gz
				rep_genes_taxonomy_krona.html
				rep_genes_taxonomy_krona.txt
				rep_genes_taxonomy.tsv
				rep_genes.fna
				rep_proteins.faa
