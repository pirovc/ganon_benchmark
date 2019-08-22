#!/bin/bash

blue_color="#abd9e9"
orange_color="#fdae61"
red_color="#d7191c"

centifuge_marker="v"
clark_marker="o"
diamond_marker="d"
ganon_marker="X"
kraken_marker="s"

evaluations=("cumu" "rank")
datasets=("challenge" "toy")
suffix="" #suffix="_1M"

for ev in ${evaluations[@]}; do
  for dt in ${datasets[@]}; do
  
    # Main results
    ganon_benchmark/scripts/plots.py \
    -i ganon_eval/new_results/20150602_bacteria_refseq_old/cami_${dt}_H01${suffix}/*.${ev}.npz \
       ganon_eval/new_results/20181219_abfv_refseq_cg_t3a/cami_${dt}_H01${suffix}/*.${ev}.npz \
       ganon_eval/new_results/20181219_abfv_refseq_all_t3a/cami_${dt}_H01${suffix}/*.${ev}.npz \
    -k 'superkingdom' 'phylum' 'class' 'order' 'family' 'genus' 'species' 'species+' \
    -l "centrifuge OLD" "clark OLD" "diamond OLD" "ganon OLD" "kraken OLD" \
       "centrifuge CG" "clark CG" "diamond CG" "ganon CG" "kraken CG" \
       "diamond ALL" "ganon ALL" \
    -c "${blue_color}" "${blue_color}" "${blue_color}" "${blue_color}" "${blue_color}" \
       "${orange_color}" "${orange_color}" "${orange_color}" "${orange_color}" "${orange_color}" \
       "${red_color}" "${red_color}" \
    -m "${centifuge_marker}" "${clark_marker}" "${diamond_marker}" "${ganon_marker}" "${kraken_marker}" \
       "${centifuge_marker}" "${clark_marker}" "${diamond_marker}" "${ganon_marker}" "${kraken_marker}" \
       "${diamond_marker}" "${ganon_marker}" \
    -o all_${dt}_${ev}${suffix}.png > all_${dt}_${ev}${suffix}.log
	
  done
done

