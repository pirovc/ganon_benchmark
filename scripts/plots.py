#!/usr/bin/env python3

import argparse, pickle
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt

def main():

	parser = argparse.ArgumentParser(description='evaluate_profiles')
	parser.add_argument('-i', metavar='<input_files>', 	required=True, 	dest="input_files", type=str, nargs="*", help="Input files (.npz)")
	parser.add_argument('-l', metavar='<labels>', 		required=False, dest="labels", 		type=str, nargs="*", default="", help="Labels for the input files")
	parser.add_argument('-c', metavar='<colors>', 		required=False, dest="colors", 		type=str, nargs="*", default="", help="Colors for the input files")
	parser.add_argument('-n', metavar='<lines>', 		required=False, dest="lines", 		type=str, nargs="*", default="", help="Lines for the input files")
	parser.add_argument('-m', metavar='<markers>', 		required=False, dest="markers", 		type=str, nargs="*", default="", help="Markers for the input files")
	parser.add_argument('-k', metavar='<ranks>', 		required=False, dest="ranks", 		type=str, nargs="*", default="", help="Evaluated ranks. Default: 'superkingdom' 'phylum' 'class' 'order' 'family' 'genus' 'species' 'species+' 'assembly'")
	parser.add_argument('-o', metavar='<output_plot>', 		required=False, dest="output_plot", type=str, help="Outut figure file")
	args = parser.parse_args()
	fixed_ranks = ['root','superkingdom','phylum','class','order','family','genus','species','species+','assembly']	
	if not args.ranks:
		selected_ranks = fixed_ranks[1:]
	else:
		selected_ranks = [r for r in args.ranks if r in fixed_ranks]

	fontsize=14

	fixed_markers = ('P', 'v', 'o', '*', '.', 'd', 's', '<', 'X', '>')
	fixed_colors = ("#e6194b","#3cb44b","#ffe119","#0082c8","#f58231","#911eb4","#46f0f0","#f032e6","#d2f53c","#fabebe","#008080","#e6beff","#aa6e28","#fffac8","#800000","#aaffc3","#808000","#ffd8b1","#000080","#808080","black")

	fig, ax = plt.subplots(nrows=1, ncols=3, sharex=False, sharey=True)

	debug = False
	tool_names = []
	for tool_id, input_file in enumerate(args.input_files):
		stats = pickle.load(open(input_file, "rb"))

		tool_name = args.labels[tool_id] if args.labels else input_file.split("/")[-1].split(".")[0]
		color = args.colors[tool_id] if args.colors else fixed_colors[tool_id]
		linestyle = args.lines[tool_id] if args.lines else  "-"
		marker = args.markers[tool_id] if args.markers else "X"

		markersize = 10

		print(tool_id, tool_name, color, linestyle, marker, input_file)
		linewidth = 2
		alpha=0.8

		#if 'cs_db' in stats:
		#	ax[0].plot([stats['cs_db'][r]/stats['cs_db']['root'] for r in selected_ranks if r in stats['cs_db']], alpha=0.7, linewidth=2, linestyle=linestyle, marker=None, markersize=0, color="black", label="")
		#else:
		#	ax[0].plot([stats['db'][r]/stats['db']['root'] for r in selected_ranks if r in stats['db']], alpha=0.7, linewidth=2, linestyle=linestyle, marker=None, markersize=0, color="black", label="")

		ax[0].plot([stats['sensitivity'][r] for r in selected_ranks if r in stats['sensitivity'] and stats['sensitivity'][r]], alpha=alpha, linewidth=linewidth, linestyle=linestyle, marker=marker, markersize=markersize, color=color, label=tool_name)
		# Precision
		ax[1].plot([stats['precision'][r] for r in selected_ranks if r in stats['precision'] and stats['precision'][r]], alpha=alpha, linewidth=linewidth, linestyle=linestyle, marker=marker, markersize=markersize, color=color, label=tool_name)
		# F1 Score
		ax[2].plot([stats['f1_score'][r] for r in selected_ranks if r in stats['f1_score'] and stats['f1_score'][r]], alpha=alpha, linewidth=linewidth, linestyle=linestyle, marker=marker, markersize=markersize, color=color, label=tool_name)


	ax[0].set_ylim([-0.02,1.02])
	
	ax[0].tick_params(labelsize=fontsize)

	ax[0].set_title("Sensitivity", fontsize=fontsize)
	ax[0].set_xticks(list(range(len(selected_ranks))))
	ax[0].set_xticklabels(selected_ranks, rotation=45, fontsize=fontsize)


	ax[1].set_ylim([-0.02,1.02])
	ax[1].set_title("Precision", fontsize=fontsize)
	ax[1].set_xticks(list(range(len(selected_ranks))))
	ax[1].set_xticklabels(selected_ranks, rotation=45, fontsize=fontsize)

	ax[2].set_ylim([-0.02,1.02])
	ax[2].set_title("F1-Score", fontsize=fontsize)
	ax[2].set_xticks(list(range(len(selected_ranks))))
	ax[2].set_xticklabels(selected_ranks, rotation=45, fontsize=fontsize)
	handles,labels = ax[2].get_legend_handles_labels()
	labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
	ax[2].legend(handles,labels, bbox_to_anchor=(1, 1), ncol=1, fontsize=fontsize)

	plt.subplots_adjust(left=None, bottom=0.15, right=0.8, top=None, wspace=0.03, hspace=None)
	if args.output_plot:
		fig.set_size_inches(18.5, 8)
		#plt.savefig(args.output_plot, dpi=300, bbox_inches='tight')
		plt.savefig(args.output_plot, dpi=1200, bbox_inches='tight')
	else:
		figManager = plt.get_current_fig_manager()
		#figManager.window.showMaximized()
		plt.show()
	
if __name__ == "__main__":
	main()
