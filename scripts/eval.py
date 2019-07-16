#!/usr/bin/env python3

import sys, gzip, pickle, argparse
from LCA import LCA
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
	parser.add_argument('-a', metavar='<output_eval_assembly>', 	required=False, dest="output_eval_assembly",  	type=str, help="Tabular evaluation output (assembly)")
	parser.add_argument('-z', metavar='<output_npz>', 				required=False, dest="output_npz",  			type=str, help="Pre-cluster sequences into rank/taxid/specialization, so they won't be splitted among bins [none,specialization name,taxid,species,genus,...] Default: none")
	args = parser.parse_args()

	if not args.ranks:
		fixed_ranks = ['root','superkingdom','phylum','class','order','family','genus','species','species+']
	else:
		fixed_ranks = args.ranks
	
	nodes, ranks = read_nodes(args.input_nodes)
	merged = read_merged(args.input_merged)
	db_file = args.input_database_profile
	gt_file = args.input_ground_truth
	res_file = args.input_results

	eval_file = open(args.output_eval,'w') if args.output_eval else sys.stdout
	eval_assembly_file = open(args.output_eval_assembly,'w') if args.output_eval_assembly else None
	eval_npz_file = open(args.output_npz, 'wb') if args.output_npz else None

	# Ground truth: readid <tab> taxid <tab> assembly
	gt = defaultdict(tuple)
	gt_leaf_taxids = set()
	for line in gzip.open(gt_file, 'rt') if gt_file.endswith(".gz") else open(gt_file,'r'):
		fields = line.rstrip().split("\t") 
		taxid = int(fields[1])
		if taxid not in nodes: 
			if taxid in merged:
				taxid = merged[taxid]
			else:
				continue #skip taxid not found
		readid = fields[0]
		assembly = fields[2]
		gt_leaf_taxids.add(taxid)
		gt[readid] = (assembly, taxid)

	# Database profile: accession <tab> len <tab> taxid <tab> assembly
	db_assembly = set() # store all assemblies db
	db_leaf_taxids = set()
	for line in gzip.open(db_file, 'rt') if db_file.endswith(".gz") else open(db_file,'r'):
		fields = line.rstrip().split("\t")
		taxid = int(fields[2])
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
		while t!=0:
			db_taxids.add(t)
			t = nodes[t]
			
	# Results: readid <tab> assembly <tab> taxid
	# if no assembly readid <tab> 0 <tab> taxid
	res = defaultdict(tuple)
	res_leaf_taxids = set()
	for line in gzip.open(res_file, 'rt') if res_file.endswith(".gz") else open(res_file,'r'):
		fields = line.rstrip().split("\t")
		taxid = int(fields[2])
		if taxid not in nodes: 
			if taxid in merged:
				taxid = merged[taxid]
			else:
				continue #skip taxid not found
		readid = fields[0].split("/")[0]
		assembly = fields[1]
		res_leaf_taxids.add(taxid)
		res[readid] = (assembly, taxid)
	
	# Filter nodes for used taxids
	filtered_nodes = {}

	for leaf_taxid in gt_leaf_taxids|res_leaf_taxids: #union gt and res leaf taxids, make lineage for LCA
		t = leaf_taxid
		while t!=0:
			filtered_nodes[t] = nodes[t]
			t = nodes[t]
	# pre-calculate LCA 
	L = LCA(filtered_nodes)

	# check lineage of the gt taxids to check at which rank it could be classified
	rank_gttaxid = {}
	for leaf_gttaxid in gt_leaf_taxids:
		t = leaf_gttaxid
		while t!=0:
			if t in db_taxids:
				rank_gttaxid[leaf_gttaxid] = rank_up_to(t, nodes, ranks, fixed_ranks)
				break
			t = nodes[t]


	# 1 - check if there is an accession (uniq assignment) and if it's correct
	# 2 - check if taxid matches
	#	if taxid tool = taxid gt -> correct identification at taxonomic level
	#	if lca = tool results -> TP, meaning it got it right in a lower taxonomic level, read was more specific
	#	if lca = gt -> FP, meaning the classification was too specific
	
	stats = {'classified': 0, 'unclassified': 0, 'tp': 0, 'fp': 0}
	classified_ranks = defaultdict(int)
	classified_ranks_assembly = 0
	db_ranks = defaultdict(int)
	db_ranks_assembly = 0
	gt_ranks = defaultdict(int)
	gt_ranks_assembly = 0
	tp_direct_ranks = defaultdict(int)
	tp_direct_ranks_assembly = 0
	fp_direct_ranks = defaultdict(int) 
	fp_direct_ranks_asssembly = 0
	# taxonomic exclusive stats
	lca_lower_ranks = defaultdict(int)
	lca_higher_ranks = defaultdict(int)
	
	for readid, (gt_assembly, gt_taxid) in gt.items():

		# Check if there's assembly id on ground truth and account for it
		if gt_assembly!="0":
			gt_ranks_assembly+=1
			if gt_assembly in db_assembly: #if assembly is present in the database (=could be classified)
				db_ranks_assembly+=1
		
		# account for taxonomic level available in gt
		gt_ranks[rank_up_to(gt_taxid, nodes, ranks, fixed_ranks)]+=1

		#if rank level taxid is present in the database (=could be classified)
		if gt_taxid in rank_gttaxid: 
			db_ranks[rank_gttaxid[gt_taxid]]+=1
			
		if readid in res.keys(): #if read is classified
			res_assembly = res[readid][0]
			res_taxid = res[readid][1]
			
			# has a unique assembly classification
			if res_assembly!="0": 
				classified_ranks_assembly+=1
				if res_assembly == gt_assembly: #it is correct
					tp_direct_ranks_assembly+=1
				else:
					fp_direct_ranks_asssembly+=1
			
			# Always accounts the taxid clasification (even with assembly)
			r = rank_up_to(res_taxid, nodes, ranks, fixed_ranks)
			classified_ranks[r]+=1
			if r=="root": # root classification is equal to false
				fp_direct_ranks[r]+=1
			else:
				if res_taxid == gt_taxid: #tp -> perfect classification
					tp_direct_ranks[r]+=1
				else:
					lca = L(gt_taxid,res_taxid)
					if lca==res_taxid: # tp -> conservative classification (gt is lower on tree)
						lca_lower_ranks[r]+=1
					elif lca==gt_taxid: # fp -> classification to specific (gt is higher on tree)
						lca_higher_ranks[r]+=1
					else: # fp -> lca is higher than gt and res
						fp_direct_ranks[r]+=1
		else:
			stats['unclassified']+=1

	stats['classified'] = len(gt) - stats['unclassified']
	stats['tp'] = sum(tp_direct_ranks.values()) + sum(lca_lower_ranks.values())
	stats['fp'] = stats['classified'] - stats['tp']
	
	final_stats = defaultdict(dict)
	header = ["-","db","gt","classified","tp","fp","tp_direct","tp_lca_lower", "fp_direct", "fp_lca_higher", "sensitivity_max_db", "sensitivity", "precision",  "f1_score"] 
	
	# assembly stats
	if eval_assembly_file:
		print("\t".join(header), file=eval_assembly_file)
		tp = tp_direct_ranks_assembly
		fp = fp_direct_ranks_asssembly
		sens = tp/len(gt)
		sens_max = tp/float(db_ranks_assembly) if db_ranks_assembly>0 else 0
		prec = tp/classified_ranks_assembly if classified_ranks_assembly>0 else 0
		f1s = (2*sens*prec)/float(sens+prec) if sens+prec>0 else 0
		print("assemly", db_ranks_assembly, gt_ranks_assembly, classified_ranks_assembly, tp, fp, tp, 0, fp, 0, sens_max, sens, prec, f1s, sep="\t", file=eval_assembly_file)
		eval_assembly_file.close()
		if eval_npz_file:
			final_stats['db']['assembly'] = db_ranks_assembly
			final_stats['gt']['assembly'] = gt_ranks_assembly
			final_stats['classified']['assembly'] = classified_ranks_assembly
			final_stats['tp']['assembly'] = tp
			final_stats['fp']['assembly'] = fp
			final_stats['tp_direct']['assembly'] = tp
			final_stats['tp_lca_lower']['assembly'] = 0
			final_stats['fp_direct']['assembly'] = fp
			final_stats['fp_lca_higher']['assembly'] = 0
			final_stats['sensitivity_max_db']['assembly'] = sens_max
			final_stats['sensitivity']['assembly'] = sens
			final_stats['precision']['assembly'] = prec
			final_stats['f1_score']['assembly'] = f1s

	# taxonomic stats
	print("\t".join(header), file=eval_file)
	print("-","-",len(gt), stats['classified'], stats['tp'], stats['fp'], sum(tp_direct_ranks.values()), sum(lca_lower_ranks.values()), sum(fp_direct_ranks.values()), sum(lca_higher_ranks.values()), sep="\t", file=eval_file)
	cs_tp = 0
	cs_fp = 0
	cs_class = 0
	cs_gt = 0
	cs_db = 0
	for fr in fixed_ranks[::-1]:
		tp = tp_direct_ranks[fr] + lca_lower_ranks[fr]
		fp = fp_direct_ranks[fr] + lca_higher_ranks[fr]
		cs_tp+=tp
		cs_fp+=fp
		cs_class+=classified_ranks[fr]
		cs_gt+=gt_ranks[fr]

		if fr=="root": #if root, all available
			db_ranks[fr]=len(gt)
		else:
			db_ranks[fr]+=cs_db
			cs_db=db_ranks[fr] # make it cumulative
			
		sens = cs_tp/len(gt)
		sens_max = cs_tp/float(db_ranks[fr]) if db_ranks[fr]>0 else 0
		prec = cs_tp/float(cs_class) if cs_class>0 else 0
		f1s = (2*sens*prec)/float(sens+prec) if sens+prec>0 else 0
		
		print(fr, db_ranks[fr], gt_ranks[fr], classified_ranks[fr], tp, fp, tp_direct_ranks[fr], lca_lower_ranks[fr], fp_direct_ranks[fr], lca_higher_ranks[fr], sens_max, sens, prec, f1s, sep="\t", file=eval_file)
		
		if eval_npz_file:
			final_stats['db'][fr] = db_ranks[fr]
			final_stats['gt'][fr] = gt_ranks[fr]
			final_stats['classified'][fr] = classified_ranks[fr]
			final_stats['tp'][fr] = tp
			final_stats['fp'][fr] = fp
			final_stats['tp_direct'][fr] = tp_direct_ranks[fr]
			final_stats['tp_lca_lower'][fr] = lca_lower_ranks[fr]
			final_stats['fp_direct'][fr] = fp_direct_ranks[fr]
			final_stats['fp_lca_higher'][fr] = lca_higher_ranks[fr]
			final_stats['sensitivity_max_db'][fr] = sens_max
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
	while taxid!=0:
		if(ranks[taxid] in fixed_ranks):
			#everything below species (not being assembly) is counted as species+
			if original_rank!="species" and ranks[taxid]=="species":
				return "species+"
			else:
				return ranks[taxid]
		taxid = nodes[taxid]
	return "root" #no standard rank identified
	
def read_nodes(nodes_file):
	# READ nodes -> fields (1:TAXID 2:PARENT_TAXID)
	nodes = {}
	ranks = {} # Only collect ranks in case pre_cluster is required
	with open(nodes_file,'r') as fnodes:
		for line in fnodes:
			taxid, parent_taxid, rank, _ = line.split('\t|\t',3)
			ranks[int(taxid)] = rank
			nodes[int(taxid)] = int(parent_taxid)
	nodes[1] = 0 #Change parent taxid of the root node to 0 (it's usually 1 and causes infinite loop later)
	return nodes, ranks

def read_merged(merged_file):
	# READ nodes -> fields (1:OLD TAXID 2:NEW TAXID)
	merged = {}
	if merged_file:
		with open(merged_file,'r') as fmerged:
			for line in fmerged:
				old_taxid, new_taxid, _ = line.rstrip().split('\t|',2)
				merged[int(old_taxid)] = int(new_taxid)
	return merged
	
if __name__ == "__main__":
	main()
