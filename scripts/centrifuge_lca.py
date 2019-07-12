import sys
from LCA import LCA
from collections import defaultdict

def main():
	unique_taxids = set()
	rmatches = defaultdict(list)
	for line in open(sys.argv[1]):
		#readid, taxid, kmercount = line.rstrip().split("\t")
		readid, assemblyid, taxid = line.rstrip().split("\t")
		taxid = int(taxid)
		unique_taxids.add(taxid)
		rmatches[readid].append((assemblyid,taxid))

	nodes = read_nodes(sys.argv[2])
	# Filter nodes for used taxids
	filtered_nodes = {}
	for leaf_taxid in unique_taxids:
		t = leaf_taxid
		while t!=0:
			try:
				filtered_nodes[t] = nodes[t]
				t = nodes[t]
			except KeyError: #if taxid is missing, link to root
				filtered_nodes[t] = 1
				t = 0
			
	L = LCA(filtered_nodes)
	
	for readid, matches in rmatches.items():
		if len(matches)==1:
			print(readid, matches[0][0], matches[0][1], sep="\t")
		else:
			tmp_lca = L(matches[0][1],matches[1][1])
			for i in range(2, len(matches)):
				tmp_lca = L(tmp_lca,matches[i][1])
			print(readid, "0", tmp_lca, sep="\t")		

def read_nodes(nodes_file):
	# READ nodes -> fields (1:TAXID 2:PARENT_TAXID)
	nodes = {}
	with open(nodes_file,'r') as fnodes:
		for line in fnodes:
			taxid, parent_taxid, rank, _ = line.split('\t|\t',3)
			nodes[int(taxid)] = int(parent_taxid)
	nodes[1] = 0 #Change parent taxid of the root node to 0 (it's usually 1 and causes infinite loop later)
	return nodes
	
if __name__ == "__main__":
	main()
