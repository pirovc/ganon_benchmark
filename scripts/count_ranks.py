#!/usr/bin/env python3

import sys
from collections import defaultdict

def main():
	input_nodes = sys.argv[1]
	acc_len_taxid = sys.argv[2]

	fixed_ranks = ['root','superkingdom','phylum','class','order','family','genus','species','species+']
	nodes, ranks = read_nodes(input_nodes)

	rank_count = defaultdict(int)
	with open(acc_len_taxid,'r') as file:
		for line in file:
			acc, length, taxid, assembly = line.split('\t',3)
			rank, changed_taxid = rank_up_to(taxid, nodes, ranks, fixed_ranks)
			rank_count[rank]+=1

	print(rank_count)

# find the closest fixed rank
def rank_up_to(taxid, nodes, ranks, fixed_ranks):
	if taxid not in ranks: return "NOT-FOUND", "0"
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


if __name__ == "__main__":
	main()