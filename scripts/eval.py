#!/usr/bin/env python3

import sys, gzip, pickle, argparse
from LCA import LCA
from collections import defaultdict

def main():

	parser = argparse.ArgumentParser(prog='ganon benchmark evaluation', conflict_handler="resolve", add_help=True)
	parser.add_argument('-k', metavar='<ranks>', 					required=False, dest="ranks", 		 			type=str, nargs="*", default="", help="Evaluated ranks. Default: 'superkingdom' 'phylum' 'class' 'order' 'family' 'genus' 'species' 'species+'")
	parser.add_argument('-i', metavar='<input_results>', 			required=True,  dest="input_results",  			type=str, help="readid <tab> assembly id (or 0) <tab> taxid")
	parser.add_argument('-g', metavar='<input_ground_truth>', 		required=True,  dest="input_ground_truth",  	type=str, help="readid <tab> assembly id (or 0) <tab> taxid")
	parser.add_argument('-d', metavar='<input_database_profile>',	required=True,  dest="input_database_profile", 	type=str, help="accession id <tab> seq. length <tab> taxid [<tab> assembly id]")
	parser.add_argument('-n', metavar='<input_nodes>', 				required=True,  dest="input_nodes",  			type=str, help="nodes.dmp - all input taxids will be normalized by this version")
	parser.add_argument('-m', metavar='<input_merged>', 			required=True,  dest="input_merged",  			type=str, help="merged.dmp - all input taxids will be normalized by this version")
	parser.add_argument('-c','--output-tab-cumulative', metavar='<output_tab_cumu>', 	required=False, dest="output_tab_cumu",  			type=str, help="Output file for evaluation in npz (cumulative mode). Use - for STDOUT")
	parser.add_argument('--output-npz-cumulative', 		metavar='<output_npz_cumu>', 	required=False, dest="output_npz_cumu",  			type=str, help="Output file for evaluation in npz (cumulative mode)")
	parser.add_argument('-r', '--output-tab-rank', 		metavar='<output_tab_rank>', 	required=False, dest="output_tab_rank",  			type=str, help="Output file for evaluation in npz (rank mode). Use - for STDOUT")
	parser.add_argument('--output-npz-rank', 			metavar='<output_npz_rank>', 	required=False, dest="output_npz_rank",  			type=str, help="Output file for evaluation in npz (rank mode)")
	args = parser.parse_args()

	if not args.ranks:
		fixed_ranks = ['root','superkingdom','phylum','class','order','family','genus','species','species+']
	else:
		fixed_ranks = ['root'] + args.ranks
	
	nodes, ranks = read_nodes(args.input_nodes)
	merged = read_merged(args.input_merged)
	gt_file = args.input_ground_truth
	db_file = args.input_database_profile
	res_file = args.input_results

	if args.output_tab_cumu=="-":
		output_tab_cumu = sys.stdout 
	elif args.output_tab_cumu:
		output_tab_cumu = open(args.output_tab_cumu,'w')
	else:
		output_tab_cumu = None
	output_npz_cumu = open(args.output_npz_cumu,'wb') if args.output_npz_cumu else None

	if args.output_tab_rank=="-":
		output_tab_rank = sys.stdout 
	elif args.output_tab_rank:
		output_tab_rank = open(args.output_tab_rank,'w')
	else:
		output_tab_rank = None
	output_npz_rank = open(args.output_npz_rank,'wb') if args.output_npz_rank else None

	# Ground truth: readid <tab> assembly <tab> taxid
	gt = defaultdict(tuple)
	gt_leaf_taxids = set()
	for line in gzip.open(gt_file, 'rt') if gt_file.endswith(".gz") else open(gt_file,'r'):
		if line[0]=="@": continue
		fields = line.rstrip().split("\t") 
		taxid = check_taxid(fields[2], gt_file, nodes, merged)
		if taxid is None: continue # taxid not found
		readid = fields[0]
		assembly = fields[1]
		gt_leaf_taxids.add(taxid)
		gt[readid] = (assembly, taxid)

	# Database profile: accession <tab> len <tab> taxid <tab> assembly
	db_assembly = set() # store all assemblies db
	db_leaf_taxids = set() # store all leaf taxids
	for line in gzip.open(db_file, 'rt') if db_file.endswith(".gz") else open(db_file,'r'):
		if line[0]=="@": continue
		fields = line.rstrip().split("\t")
		taxid = check_taxid(fields[2], db_file, nodes, merged)
		if taxid is None: continue # taxid not found
		db_leaf_taxids.add(taxid)
		if len(fields)>3: db_assembly.add(fields[3])
	# store all taxid lineage on db to check gt later
	db_taxids = set() 
	for leaf_dbtaxid in db_leaf_taxids:
		t = leaf_dbtaxid
		while t!="0":
			db_taxids.add(t)
			t = nodes[t]
			
	# Results: readid <tab> assembly <tab> taxid
	# if no assembly readid <tab> 0 <tab> taxid
	res = defaultdict(tuple)
	res_leaf_taxids = set()
	for line in gzip.open(res_file, 'rt') if res_file.endswith(".gz") else open(res_file,'r'):
		if line[0]=="@": continue
		fields = line.rstrip().split("\t")
		taxid = check_taxid(fields[2], res_file, nodes, merged)
		if taxid is None: continue # taxid not found
		readid = fields[0]
		assembly = fields[1]
		res_leaf_taxids.add(taxid)
		res[readid] = (assembly, taxid)


	# Cumulative evaluation
	if output_tab_cumu or output_npz_cumu:

		# Filter nodes for used taxids to generate LCA
		filtered_nodes = {}
		for leaf_taxid in gt_leaf_taxids|res_leaf_taxids: #union gt and res leaf taxids, make lineage for LCA
			t = leaf_taxid
			while t!="0":
				filtered_nodes[t] = nodes[t]
				t = nodes[t]
		# pre-calculate LCA 
		L = LCA(filtered_nodes)

		# check lineage of the gt taxids to check at which rank it could be classified given the database
		rank_gttaxid = {}
		for leaf_gttaxid in gt_leaf_taxids:
			t = leaf_gttaxid
			if not rank_gttaxid.get(leaf_gttaxid): # if not yet calculated
				while t!="0":
					if t in db_taxids:
						rank_gttaxid[leaf_gttaxid], _ = rank_up_to(t, nodes, ranks, fixed_ranks)
						break
					t = nodes[t]

		# 1 - check if there is an accession (uniq assignment) and if it's correct
		# 2 - check if taxid matches
		#	if taxid tool = taxid gt -> correct identification at taxonomic level
		#	if lca = tool results -> TP, meaning it got it right in a lower taxonomic level, read was more specific
		#	if lca = gt -> FP, meaning the classification was too specific
		cumulative_eval(res, gt, nodes, ranks, fixed_ranks, L, rank_gttaxid, output_tab_cumu, output_npz_cumu)

	# Rank evaluation
	if output_tab_rank or output_npz_rank:
		rank_eval(res, gt, nodes, ranks, fixed_ranks, db_assembly, db_taxids, output_tab_rank, output_npz_rank)

	
	if output_tab_cumu: output_tab_cumu.close()
	if output_npz_cumu: output_npz_cumu.close()
	if output_tab_rank: output_tab_rank.close()
	if output_npz_rank: output_npz_rank.close()


# find the closest fixed rank
def rank_up_to(taxid, nodes, ranks, fixed_ranks):
	original_rank = ranks[taxid]
	original_taxid = taxid
	while taxid!="0":
		if(ranks[taxid] in fixed_ranks):
			#everything below species (not being assembly) is counted as species+
			if original_rank!="species" and ranks[taxid]=="species" and "species+" in fixed_ranks:
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
				old_taxid, new_taxid, _ = line.split('\t|',2)
				merged[old_taxid] = new_taxid.rstrip('\t')
	return merged

def print_log(text):
	sys.stderr.write(text+"\n")

def check_taxid(taxid, file, nodes, merged):
	if taxid not in nodes: 
		if taxid in merged:
			print_log("taxid changed [" + taxid + " -> " + merged[taxid] + ", " + file + "]")
			taxid = merged[taxid]
		else:
			print_log("taxid not found, entry ignored [" + taxid + ", " + file + "]")
			return None
	return taxid

def get_lineage(nodes, ranks, taxid, fixed_ranks):
    lin = [""]*len(fixed_ranks)
    t = taxid
    while t!="1": 
    	r,tx = rank_up_to(t, nodes, ranks, fixed_ranks)
    	lin[fixed_ranks.index(r)]=tx
    	if tx=="1": break
    	t = nodes[tx]
    return lin

def cumulative_eval(res,gt,nodes, ranks, fixed_ranks, L, rank_gttaxid, output_tab_cumu, output_npz_cumu):

	stats = {'classified': 0, 'unclassified': 0, 'tp': 0, 'fp': 0}
	classified_ranks = defaultdict(int)
	db_ranks = defaultdict(int)
	gt_ranks = defaultdict(int)
	tp_direct_ranks = defaultdict(int)
	fp_direct_ranks = defaultdict(int) 
	lca_lower_ranks = defaultdict(int)
	lca_higher_ranks = defaultdict(int)
	
	for readid, (_, gt_taxid) in gt.items():

		# account for taxonomic level available in gt
		leaf_rank_gt, _ = rank_up_to(gt_taxid, nodes, ranks, fixed_ranks)
		gt_ranks[leaf_rank_gt]+=1

		#if rank level taxid is present in the database (=could be classified)
		# define at which level a certain taxid could be classified given this db
		if gt_taxid in rank_gttaxid: 
			db_ranks[rank_gttaxid[gt_taxid]]+=1
			
		if readid in res.keys(): #if read is classified
			res_taxid = res[readid][1]
					
			# taxonomic clasification
			r, _ = rank_up_to(res_taxid, nodes, ranks, fixed_ranks)
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

	total_reads_gt = len(gt)
	stats['classified'] = total_reads_gt - stats['unclassified']
	stats['tp'] = sum(tp_direct_ranks.values()) + sum(lca_lower_ranks.values())
	stats['fp'] = stats['classified'] - stats['tp']
	
	final_stats = defaultdict(dict)
	header = ["-","gt","db","class","tp","fp","cs_db","cs_class","cs_tp","cs_fp","tp_direct","tp_lca_lower", "fp_direct", "fp_lca_higher", "sens_max_db", "sens", "prec",  "f1s"] 

	# taxonomic stats
	print("\t".join(header), file=output_tab_cumu)
	print("SUM", total_reads_gt, sum(db_ranks.values()), stats['classified'], stats['tp'], stats['fp'],"-","-","-","-", sum(tp_direct_ranks.values()), sum(lca_lower_ranks.values()), sum(fp_direct_ranks.values()), sum(lca_higher_ranks.values()),"-","-","-","-", sep="\t", file=output_tab_cumu)
	cs_db = 0
	cs_class = 0
	cs_tp = 0
	cs_fp = 0

	for fr in fixed_ranks[::-1]:
		tp = tp_direct_ranks[fr] + lca_lower_ranks[fr]
		fp = fp_direct_ranks[fr] + lca_higher_ranks[fr]
		
		cs_db+=db_ranks[fr] # make it cumulative
		cs_class+=classified_ranks[fr]
		cs_tp+=tp
		cs_fp+=fp
		
		#if root, all available
		if fr=="root": cs_db=total_reads_gt		
			
		sens_max_db = cs_tp/float(cs_db) if cs_db>0 else 0
		sens = cs_tp/total_reads_gt
		prec = cs_tp/float(cs_class) if cs_class>0 else 0
		f1s = (2*sens*prec)/float(sens+prec) if sens+prec>0 else 0
		
		print(fr, gt_ranks[fr], db_ranks[fr], classified_ranks[fr], tp, fp, cs_db, cs_class, cs_tp, cs_fp, tp_direct_ranks[fr], lca_lower_ranks[fr], fp_direct_ranks[fr], lca_higher_ranks[fr], "%.5f" % sens_max_db, "%.5f" % sens, "%.5f" % prec, "%.5f" % f1s, sep="\t", file=output_tab_cumu)
		
		if output_npz_cumu:
			final_stats['gt'][fr] = gt_ranks[fr]
			final_stats['db'][fr] = db_ranks[fr]
			final_stats['classified'][fr] = classified_ranks[fr]
			final_stats['tp'][fr] = tp
			final_stats['fp'][fr] = fp
			final_stats['tp_direct'][fr] = tp_direct_ranks[fr]
			final_stats['tp_lca_lower'][fr] = lca_lower_ranks[fr]
			final_stats['fp_direct'][fr] = fp_direct_ranks[fr]
			final_stats['fp_lca_higher'][fr] = lca_higher_ranks[fr]
			final_stats['sensitivity_max_db'][fr] = sens_max_db
			final_stats['sensitivity'][fr] = sens
			final_stats['precision'][fr] = prec
			final_stats['f1_score'][fr] = f1s
	
	if output_npz_cumu: pickle.dump(final_stats, output_npz_cumu)

def rank_eval(res, gt, nodes, ranks, fixed_ranks, db_assembly, db_taxids, output_tab_rank, output_npz_rank):

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
	
	# datastructure to save lineages already computed (saves time)
	taxid_lineage = defaultdict(lambda: ['']*len(fixed_ranks))

	for readid, (gt_assembly, gt_taxid) in gt.items():

		# Check if there's assembly id on ground truth and account for it
		if gt_assembly!="0":
			gt_ranks_assembly+=1
			if gt_assembly in db_assembly: #if assembly is present in the database (=could be classified)
				db_ranks_assembly+=1

		# make lineage if not made yet
		if not taxid_lineage.get(gt_taxid):
			taxid_lineage[gt_taxid] = get_lineage(nodes,ranks,gt_taxid,fixed_ranks)

		# Check if the taxid is on ground truth and account for it
		for idx,fr in enumerate(fixed_ranks):
			if taxid_lineage[gt_taxid][idx]!="": gt_ranks[fr]+=1
			if taxid_lineage[gt_taxid][idx] in db_taxids: # if the is present in the database (=could be classified)
				db_ranks[fr]+=1

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
	header = ["-", "gt","db","class","tp","fp", "sens_max_db", "sens", "prec", "f1s"] 
	
	print("\t".join(header), file=output_tab_rank)

	sens_assembly = tp_ranks_assembly/total_reads_gt
	sens_max_assembly = tp_ranks_assembly/float(db_ranks_assembly) if db_ranks_assembly>0 else 0
	prec_assembly = tp_ranks_assembly/classified_ranks_assembly if classified_ranks_assembly>0 else 0
	f1s_assembly = (2*sens_assembly*prec_assembly)/float(sens_assembly+prec_assembly) if sens_assembly+prec_assembly>0 else 0
	print("assembly", gt_ranks_assembly, db_ranks_assembly, classified_ranks_assembly, tp_ranks_assembly, fp_ranks_assembly, "%.5f" % sens_max_assembly, "%.5f" % sens_assembly, "%.5f" % prec_assembly, "%.5f" % f1s_assembly, sep="\t", file=output_tab_rank)

	for fr in fixed_ranks[::-1]:
		tp = tp_ranks[fr]
		fp = fp_ranks[fr]
		sens = tp/total_reads_gt
		sens_max = tp/float(db_ranks[fr]) if db_ranks[fr]>0 else 0
		prec = tp/classified_ranks[fr] if classified_ranks[fr]>0 else 0
		f1s = (2*sens*prec)/float(sens+prec) if sens+prec>0 else 0
		
		print(fr, gt_ranks[fr], db_ranks[fr], classified_ranks[fr], tp, fp, "%.5f" % sens_max, "%.5f" % sens, "%.5f" % prec, "%.5f" % f1s, sep="\t", file=output_tab_rank)
		
		if output_npz_rank:
			final_stats['db'][fr] = db_ranks[fr]
			final_stats['gt'][fr] = gt_ranks[fr]
			final_stats['classified'][fr] = classified_ranks[fr]
			final_stats['tp'][fr] = tp
			final_stats['fp'][fr] = fp
			final_stats['sensitivity'][fr] = sens
			final_stats['precision'][fr] = prec
			final_stats['f1_score'][fr] = f1s
	
	if output_npz_rank:  pickle.dump(final_stats, output_npz_rank)

if __name__ == "__main__":
	main()
