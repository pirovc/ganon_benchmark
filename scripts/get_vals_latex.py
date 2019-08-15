#!/usr/bin/env python3

import pickle, argparse, sys
from collections import OrderedDict

def main():
	parser = argparse.ArgumentParser(description='evaluate_profiles')
	parser.add_argument('-i', metavar='<input_files>', 	required=True, 	dest="input_files", type=str, nargs="*", help="Input files (.npz)")
	parser.add_argument('-c','--check', default=False, action='store_true', help='Verbose mode for ganon')
	parser.add_argument('-k', metavar='<ranks>', 		required=False, dest="ranks", 		type=str, nargs="*", default="", help="Evaluated ranks. Default: 'superkingdom' 'phylum' 'class' 'order' 'family' 'genus' 'species' 'species+' 'assembly'")
	args = parser.parse_args()

	if not args.ranks:
		selected_ranks = ['root','superkingdom','phylum','class','order','family','genus','species','species+']
	else:
		selected_ranks = args.ranks

	stats=OrderedDict()
	for input_file in sorted(args.input_files):
		s = pickle.load(open(input_file, "rb"))
		tool_name = input_file.split("/")[-1].split(".")[0].split("-")[0]
		stats[tool_name] = s

	for tool_name, s in stats.items():
		if args.check: 
			print(tool_name)
		else:
			print(tool_name, end="\t& ")
		for r in selected_ranks:
			if args.check: 
				print('precision', r, s['precision'][r], sep="\t")
			else:
				print("%.2f\\%%" % (s['precision'][r]*100), end="\t& ")
		for r in selected_ranks:
			if args.check: 
				print('sensitivity', r, s['sensitivity'][r], sep="\t")
			else:
				print("%.2f\\%%" % (s['sensitivity'][r]*100), end="\t& ")
		for r in selected_ranks:
			if args.check: 
				print('f1_score', r, s['f1_score'][r], sep="\t")
			else:
				print("%.2f\\%%" % (s['f1_score'][r]*100), end=" \\\\ " if selected_ranks.index(r)==len(selected_ranks)-1 else "\t& ")
		print("")

if __name__ == "__main__":
	main()
