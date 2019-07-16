#!/usr/bin/env python3

import sys, gzip, pickle, argparse
from collections import defaultdict

def main():

	parser = argparse.ArgumentParser(prog='ganon benchmark evaluation', conflict_handler="resolve", add_help=True)
	parser.add_argument('-k', metavar='<ranks>', 					required=False, dest="ranks", 		 			type=str, nargs="*", default="", help="Evaluated ranks. Default: 'root' 'superkingdom' 'phylum' 'class' 'order' 'family' 'genus' 'species' 'species+'")
	parser.add_argument('-r', metavar='<input_results>', 			required=True,  dest="input_results",  			type=str, help="readid <tab> assembly id (or 0) <tab> taxid")
	parser.add_argument('-g', metavar='<input_ground_truth>', 		required=True,  dest="input_ground_truth",  	type=str, help="readid <tab> taxid <tab> assembly id (or 0)")
	parser.add_argument('-d', metavar='<input_database_profile>',	required=True,  dest="input_database_profile", 	type=str, help="accession id <tab> seq. length <tab> taxid [<tab> assembly id]")
	parser.add_argument('-n', metavar='<input_nodes>', 				required=True,  dest="input_nodes",  			type=str, help="nodes.dmp - all input taxids will be normalized by this version")
	parser.add_argument('-m', metavar='<input_merged>', 			required=True,  dest="input_merged",  			type=str, help="merged.dmp - all input taxids will be normalized by this version")
	parser.add_argument('-o', metavar='<output_eval>', 				required=False, dest="output_eval",  			type=str, help="Tabular evaluation output (taxonomic ranks)")
	parser.add_argument('-z', metavar='<output_npz>', 				required=False, dest="output_npz",  			type=str, help="Pre-cluster sequences into rank/taxid/specialization, so they won't be splitted among bins [none,specialization name,taxid,species,genus,...] Default: none")
	args = parser.parse_args()

	if not args.ranks:
		fixed_ranks = ['root','superkingdom','phylum','class','order','family','genus','species','species+']
	else:
		fixed_ranks = args.ranks
	
	nodes, ranks = read_nodes(args.input_nodes)
	merged = read_merged(args.input_merged)
	gt_file = args.input_ground_truth
	db_file = args.input_database_profile
	res_file = args.input_results

	eval_file = open(args.output_eval,'w') if args.output_eval else sys.stdout
	eval_npz_file = open(args.output_npz, 'wb') if args.output_npz else None

	# datastructure to save lineages (save computation)
	taxid_lineage = defaultdict(lambda: ['']*len(fixed_ranks))

	# Ground truth: readid <tab> taxid <tab> assembly
	gt = defaultdict(tuple)
	for line in gzip.open(gt_file, 'rt') if gt_file.endswith(".gz") else open(gt_file,'r'):
		fields = line.rstrip().split("\t") 
		taxid = fields[1]
		if taxid not in nodes: 
			if taxid in merged:
				taxid = merged[taxid]
			else:
				continue #skip taxid not found
		readid = fields[0]
		assembly = fields[2]
		gt[readid] = (assembly, taxid)
		
	# Database profile: accession <tab> len <tab> taxid <tab> assembly
	db_assembly = set() # store all assemblies db
	db_leaf_taxids = set()
	for line in gzip.open(db_file, 'rt') if db_file.endswith(".gz") else open(db_file,'r'):
		fields = line.rstrip().split("\t")
		taxid = fields[2]
		if taxid not in nodes: 
			if taxid in merged:
				taxid = merged[taxid]
			else:
				continue #skip taxid not found
		db_leaf_taxids.add(taxid)
		if len(fields)>3: db_assembly.add(fields[3])
	db_taxids = set() # store all taxid lineage on db to check gt later
	for leaf_dbtaxid in db_leaf_taxids:
		t = leaf_dbtaxid
		while t!="0":
			db_taxids.add(t)
			t = nodes[t]

	# Results: readid <tab> assembly <tab> taxid
	# if no assembly readid <tab> 0 <tab> taxid
	res = defaultdict(tuple)
	for line in gzip.open(res_file, 'rt') if res_file.endswith(".gz") else open(res_file,'r'):
		fields = line.rstrip().split("\t")
		taxid = fields[2]
		if taxid not in nodes: 
			if taxid in merged:
				taxid = merged[taxid]
			else:
				continue #skip taxid not found
		readid = fields[0].split("/")[0]
		assembly = fields[1]
		res[readid] = (assembly, taxid)
	
	stats = {'classified': 0, 'unclassified': 0, 'tp': 0, 'fp': 0}
	classified_ranks = defaultdict(int)
	classified_ranks_assembly = 0
	db_ranks = defaultdict(int)
	db_ranks_assembly = 0
	gt_ranks = defaultdict(int)
	gt_ranks_assembly = 0
	tp_ranks = defaultdict(int)
	tp_ranks_assembly = 0
	fp_ranks = defaultdict(int) 
	fp_ranks_assembly = 0
	
	for readid, (gt_assembly, gt_taxid) in gt.items():

		# make lineage if not made yet
		if not taxid_lineage.get(gt_taxid):
			taxid_lineage[gt_taxid] = get_lineage(nodes,ranks,gt_taxid,fixed_ranks)


		# Check if there's assembly id on ground truth and account for it
		if gt_assembly!="0":
			gt_ranks_assembly+=1
			if gt_assembly in db_assembly: #if assembly is present in the database (=could be classified)
				db_ranks_assembly+=1

		for idx,fr in enumerate(fixed_ranks):
			if taxid_lineage[gt_taxid][idx] in db_taxids:
				db_ranks[fr]+=1

		gt_leaf_rank, _ = rank_up_to(gt_taxid, nodes, ranks, fixed_ranks)
		# account for taxonomic level available in gt
		gt_ranks[gt_leaf_rank]+=1
			
		if readid in res.keys(): #if read is classified
			res_assembly = res[readid][0]
			res_taxid = res[readid][1]

			# has a unique assembly classification
			if res_assembly!="0": 
				classified_ranks_assembly+=1
				if res_assembly == gt_assembly: #it is correct
					tp_ranks_assembly+=1
				else:
					fp_ranks_assembly+=1


			# make lineage if not made yet
			if not taxid_lineage.get(res_taxid):
				taxid_lineage[res_taxid] = get_lineage(nodes,ranks,res_taxid,fixed_ranks)

			res_leaf_rank, _ = rank_up_to(res_taxid, nodes, ranks, fixed_ranks)

			# compare from classification rank up
			for idx,fr in enumerate(fixed_ranks[:fixed_ranks.index(res_leaf_rank)+1]):
				classified_ranks[fr]+=1
				if taxid_lineage[gt_taxid][idx]==taxid_lineage[res_taxid][idx]:
					tp_ranks[fr]+=1
				else:
					fp_ranks[fr]+=1



		else:
			stats['unclassified']+=1

	total_reads_gt = len(gt)
	stats['classified'] = total_reads_gt - stats['unclassified']
	
	final_stats = defaultdict(dict)
	header = ["total_reads_gt:"+str(total_reads_gt), "gt","db","class","tp","fp", "sens_max_db", "sens", "prec", "f1s"] 
	
	print("\t".join(header), file=eval_file)

	sens_assembly = tp_ranks_assembly/total_reads_gt
	sens_max_assembly = tp_ranks_assembly/float(db_ranks_assembly) if db_ranks_assembly>0 else 0
	prec_assembly = tp_ranks_assembly/classified_ranks_assembly if classified_ranks_assembly>0 else 0
	f1s_assembly = (2*sens_assembly*prec_assembly)/float(sens_assembly+prec_assembly) if sens_assembly+prec_assembly>0 else 0
	print("assembly", gt_ranks_assembly, db_ranks_assembly, classified_ranks_assembly, tp_ranks_assembly, fp_ranks_assembly, "%.5f" % sens_max_assembly, "%.5f" % sens_assembly, "%.5f" % prec_assembly, "%.5f" % f1s_assembly, sep="\t", file=eval_file)

	for fr in fixed_ranks[::-1]:
		tp = tp_ranks[fr]
		fp = fp_ranks[fr]
		sens = tp/total_reads_gt
		sens_max = tp/float(db_ranks[fr]) if db_ranks[fr]>0 else 0
		prec = tp/classified_ranks[fr] if classified_ranks[fr]>0 else 0
		f1s = (2*sens*prec)/float(sens+prec) if sens+prec>0 else 0
		
		print(fr, gt_ranks[fr], db_ranks[fr], classified_ranks[fr], tp, fp, "%.5f" % sens_max, "%.5f" % sens, "%.5f" % prec, "%.5f" % f1s, sep="\t", file=eval_file)
		
		if eval_npz_file:
			final_stats['db'][fr] = db_ranks[fr]
			final_stats['gt'][fr] = gt_ranks[fr]
			final_stats['classified'][fr] = classified_ranks[fr]
			final_stats['tp'][fr] = tp
			final_stats['fp'][fr] = fp
			final_stats['sensitivity'][fr] = sens
			final_stats['precision'][fr] = prec
			final_stats['f1_score'][fr] = f1s
	
	if eval_npz_file: 
		pickle.dump(final_stats, eval_npz_file)
		eval_npz_file.close()

	if eval_file: eval_file.close()

# find the closest fixed rank
def rank_up_to(taxid, nodes, ranks, fixed_ranks):
	original_rank = ranks[taxid]
	original_taxid = taxid
	while taxid!="0":
		if(ranks[taxid] in fixed_ranks):
			#everything below species (not being assembly) is counted as species+
			if original_rank!="species" and ranks[taxid]=="species":
				return "species+", original_taxid
			else:
				return ranks[taxid], taxid
		taxid = nodes[taxid]
	return "root", "1" #no standard rank identified
	
def read_nodes(nodes_file):
	# READ nodes -> fields (1:TAXID 2:PARENT_TAXID)
	nodes = {}
	ranks = {} # Only collect ranks in case pre_cluster is required
	with open(nodes_file,'r') as fnodes:
		for line in fnodes:
			taxid, parent_taxid, rank, _ = line.split('\t|\t',3)
			ranks[taxid] = rank
			nodes[taxid] = parent_taxid
	nodes["1"] = "0" #Change parent taxid of the root node to 0 (it's usually 1 and causes infinite loop later)
	return nodes, ranks

def read_merged(merged_file):
	# READ nodes -> fields (1:OLD TAXID 2:NEW TAXID)
	merged = {}
	if merged_file:
		with open(merged_file,'r') as fmerged:
			for line in fmerged:
				old_taxid, new_taxid, _ = line.rstrip().split('\t|',2)
				merged[old_taxid] = new_taxid
	return merged

def get_lineage(nodes, ranks, taxid, fixed_ranks):
    lin = [""]*len(fixed_ranks)
    t = taxid
    while t!="1": 
    	r,tx = rank_up_to(t, nodes, ranks, fixed_ranks)
    	lin[fixed_ranks.index(r)]=tx
    	if tx=="1": break
    	t = nodes[tx]
    return lin

if __name__ == "__main__":
	main()
