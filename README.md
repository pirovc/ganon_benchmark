# ganon_benchmark

benchmarking pipeline used to evaluate ganon as presented in https://www.biorxiv.org/content/10.1101/406017v3

## description:

- The pipeline is divided into build and classify steps.
- The folder config/ contains examples of configuration files for building and classifying
- Configuration files used for generating the results presented in the manuscript are in files/config/. Files should be edited in a way the paths are correctly set. 

## dependencies:

- snakemake>=5.2.0
- AMBER==2.0.21-beta (https://github.com/CAMI-challenge/AMBER/archive/658a268bc86622e6f41304d38c1856b727c995ed.tar.gz)
* it's necessary to update the "tools_path" in the configuration files with path to AMBER. All dependencies are covered by conda.

## running analysis:

### 0) Prepare files

	# unpack the aux files
	gzip -d files/*refseq*/*acc_len_taxid_assembly.txt.gz

	# unpack assembly_summary.txt
	gzip -d files/downloads/*/assembly_summary.txt.gz

	# unpack taxdump files
	gzip -d files/taxdump/*.dmp.gz

	# unpack gt files
	gzip -d files/reads/cami/*/*readid_assembly_taxid.txt.gz

### 1) Download reference sequences

	# Recover reference sequences (.fna and .faa) for RefSeqOLD
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Bacteria/all.fna.tar.gz
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Bacteria/all.faa.tar.gz

	wget https://raw.githubusercontent.com/pirovc/genome_updater/a301b927362aced47ed88a20b3f396df0bf5258d/genome_updater.sh
	chmod +x genome_updater.sh
	
	# Recover reference sequences (.fna and .faa) for RefSeqCG with 32 threads
	# 39246 Files, 22GB
	genome_updater.sh -e files/downloads/abfv_refseq_cg/assembly_summary.txt -o files/downloads/abfv_refseq_cg/ -i -m -f "genomic.fna.gz,protein.faa.gz" -t 32 -v recovered
	
	# Recover reference sequences (fna and faa) for RefSeqALL with 32 threads
	# 295426 files, 270GB
	genome_updater.sh -e files/downloads/abfv_refseq_all/assembly_summary.txt -o files/downloads/abfv_refseq_all/ -i -m -f "genomic.fna.gz,protein.faa.gz" -t 32 -v recovered

### 2) Generate reference files

#### RefSeqOLD

	# .fna
	mkdir files/downloads/refseq_old_fna
	tar xf files/downloads/all.fna.tar.gz -C files/downloads/refseq_old_fna
	cat files/downloads/refseq_old_fna/*/*.fna | awk '{if(substr($0,0,1)==">"){split($0,sep,"|"); print ">" sep[4] sep[5];}else{print $0}}' > files/bacteria_refseq_old/all.fna
	
	# .faa
	mkdir files/downloads/refseq_old_faa
	tar xf files/downloads/all.faa.tar.gz -C files/downloads/refseq_old_faa
	cat files/downloads/refseq_old_faa/*/*.faa | awk '{if(substr($0,0,1)==">"){split($0,sep,"|"); print ">" sep[4] sep[5];}else{print $0}}' > files/bacteria_refseq_old/all.faa

#### RefSeq-CG

	# .fna
	zcat files/downloads/abfv_refseq_cg/recovered/files/*_genomic.fna.gz > files/abfv_refseq_cg/20181219_abfv_refseq_cg.fna
	
	# .faa
	zcat files/downloads/abfv_refseq_cg/recovered/files/*_protein.faa.gz > files/abfv_refseq_cg/20181219_abfv_refseq_cg.faa

#### RefSeq-CG-top-3

	# .fna
	cut -f 4 files/abfv_refseq_cg_t3a/20181219_abfv_refseq_cg_t3a_acc_len_taxid_assembly.txt | sort | uniq | while read f; do zcat files/downloads/abfv_refseq_cg/recovered/files/${f}*_genomic.fna.gz >> files/abfv_refseq_cg_t3a/20181219_abfv_refseq_cg_t3a.fna 2> files/abfv_refseq_cg_t3a/20181219_abfv_refseq_cg_t3a.fna.log; done
	
	# .faa
	cut -f 4 files/abfv_refseq_cg_t3a/20181219_abfv_refseq_cg_t3a_acc_len_taxid_assembly.txt | sort | uniq | while read f; do zcat files/downloads/abfv_refseq_cg/recovered/files/${f}*_protein.faa.gz >> files/abfv_refseq_cg_t3a/20181219_abfv_refseq_cg_t3a.faa 2> files/abfv_refseq_cg_t3a/20181219_abfv_refseq_cg_t3a.fna.log; done

#### RefSeq-ALL

	# .fna
	zcat files/downloads/abfv_refseq_all/recovered/files/*_genomic.fna.gz > files/abfv_refseq_all/20181219_abfv_refseq_all.fna
	
	# .faa
	zcat files/downloads/abfv_refseq_all/recovered/files/*_protein.faa.gz > files/abfv_refseq_all/20181219_abfv_refseq_all.faa

#### RefSeq-ALL-top-3

	# .fna
	cut -f 4 files/abfv_refseq_all_t3a/20181219_abfv_refseq_all_t3a_acc_len_taxid_assembly.txt | sort | uniq | while read f; do zcat files/downloads/abfv_refseq_all/recovered/files/${f}*_genomic.fna.gz >> files/abfv_refseq_all_t3a/20181219_abfv_refseq_all_t3a.fna 2> files/abfv_refseq_all_t3a/20181219_abfv_refseq_all_t3a.fna.log; done
	
	# .faa
	cut -f 4 files/abfv_refseq_all_t3a/20181219_abfv_refseq_all_t3a_acc_len_taxid_assembly.txt | sort | uniq | while read f; do zcat downloads/abfv_refseq_all/recovered/files/${f}*_protein.faa.gz >> //abfv_refseq_all_t3a/20181219_abfv_refseq_all_t3a.faa 2> files/abfv_refseq_all_t3a/20181219_abfv_refseq_all_t3a.fna.log; done

### 3) Obtaining reads

A sub-set of a random 1 million reads from the complete set used to evaluate the tools is provided in this repository under `files/reads/`. 

To get the whole read sets used, please follow the instructions at: https://data.cami-challenge.org/participate

### 4) Running the pipeline

To test the commands before running, use the addition `-npr` parameter.

	# Build all indices
	snakemake --snakefile build/Snakefile --configfile files/config/build_complete_eval.yaml --use-conda --cores 48 -k > build.out 2>&1 

	# Classify simulated reads at taxonomic level
	snakemake --snakefile classify/Snakefile --configfile files/config/classify_toy_taxid.yaml --use-conda --cores 48 -k > classify_sim_taxid.out 2>&1 
	
	# Classify simulated reads at assembly level
	snakemake --snakefile classify/Snakefile --configfile files/config/classify_toy_assembly.yaml --use-conda --cores 48 -k > classify_sim_assembly.out 2>&1 

	# Classify real reads at taxonomic level
	snakemake --snakefile classify/Snakefile --configfile files/config/classify_challenge_taxid.yaml --use-conda --cores 48 -k > classify_real_taxid.out 2>&1 

	# Plots
	scripts/plots.sh
	
## build config:

The main input files are:
 - fasta: .fasta/.fna file with sequences, header should be in the NCBI standards (>accessions.version)
 - faa: .faa file with sequences (for diamond)
 - nodes: nodes.dmp from taxonomy matching fasta/faa files
 - merged: merged.dmp from taxonomy matching fasta/faa files
 - names: names.dmp from taxonomy matching fasta/faa files
 - acc_len_taxid_assembly: aux. file with the fields: accession.version <tab> sequence len. <tab> taxid [<tab> assembly id]

### Folder structure and output files:

	workdir/
	  db_name/
	    tool_db/
	      db_config/
	         prefix.\*: database/index files
	         prefix.check: size of all essential database files
	         prefix.log: stdout and stderr from the tool run
	         prefix.time: output from "/usr/bin/time -v"
	         prefix.stats: tool <tab> database <tab> parameters <tab> time elapsed <tab> time seconds <tab> peak memory bytes <tab> index size bytes

## classify config:

The main input files are:
 - databases generated from the build step
 - fq1: .fq.gz/.fastq.gz reads
 - gt: ground truth file with the fields: accession.version <tab> assembly id <tab> taxonomic id

### Folder structure and output files:

	workdir/
	   db_name/
	     sample_name/
	       amber/: amber output folder
	         index.html: main amber results
	       amber.gt.bioboxes: ground truth for amber
	       amber.log: stdout and stderr from amber run
	       eval_merged.dmp: merged.dmp used for the evaluations
	       eval_names.dmp: names.dmp used for the evaluations
	       eval_nodes.dmp: nodes.dmp used for the evaluations
	       tool-run_config-db_config.bioboxes.gz: read id <tab> assembly id <tab> taxonomic id
	       tool-run_config-db_config.cumu.{npz,tsv}: cumulative-based evaluation results
	       tool-run_config-db_config.rank.{npz,tsv}: rank-based evaluation results
	       tool-run_config-db_config.eval.log: stdout and stderr from the run evaluation run
	       tool-run_config-db_config.log: stdout and stderr from the tool run
	       tool-run_config-db_config.time: output from "/usr/bin/time -v"
	       tool-run_config-db_config.stats: tool <tab> sample_name <tab> db_name <tab> db_config <tab> run_config <tab> time elapsed <tab> time elapsed seconds <tab> Mbp/m (or equivalent run time in \*seconds\*) <tab> peak memory bytes
	       plot_cumu.png: plot of the cumulative-based evaluation for this sample
	       plot_rank.png: plot of the rank-based evaluation for this sample

### evaluation files:

	\*.cumu.tsv = cumulative-based results
		gt: number of reads with entries on the ground truth
		db: number of reads with entries on the ground refence set used
		class: number of reads classified
		tp: number of true positives
		fp: number of false positives
		sens_max_db: sensitivity (based on db)
		sens: sensitivity
		prec: precision
		f1s: f1-score
	

	\*.rank.tsv = rank-based results
		gt: number of reads with entries on the ground truth
		db: number of reads with entries on the ground refence set used
		class: number of reads classified
		tp: number of true positives (tp_direct + tp_indirect)
		fp: number of false positives (fp_lower + fp_higher)
		cs_db: cumulative sum of db
		cs_class: cumulative sum of class
		cs_tp: cumulative sum of tp
		cs_fp: cumulative sum of fp
		tp_direct: number of direct true positives 
		tp_indirect: number of indirect true positives 
		fp_lower: number of lower false positives 
		fp_higher: number of higher false positives 
		sens_max_db: sensitivity (based on db)
		sens: sensitivity
		prec: precision
		f1s: f1-score
