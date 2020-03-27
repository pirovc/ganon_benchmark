#!/bin/bash

blue_color="#abd9e9"
orange_color="#fdae61"
red_color="#d7191c"

centifuge_marker="v"
clark_marker="o"
diamond_marker="d"
ganon_marker="X"
kraken_marker="s"
kraken2_marker="p"

evaluations=("cumu" "rank")
datasets=("challenge" "toy")
suffix="_1M"

for ev in ${evaluations[@]}; do
  for dt in ${datasets[@]}; do
  
    # Main results
    ganon_benchmark/scripts/plots.py \
    -i ganon_eval/new_results/20150602_bacteria_refseq_old/cami_${dt}_H01${suffix}/*.${ev}.npz \
       ganon_eval/new_results/20181219_abfv_refseq_cg_t3a/cami_${dt}_H01${suffix}/*.${ev}.npz \
       ganon_eval/new_results/20181219_abfv_refseq_all_t3a/cami_${dt}_H01${suffix}/*.${ev}.npz \
    -k 'phylum' 'class' 'order' 'family' 'genus' 'species' 'species+' \
    -l "centrifuge OLD" "clark OLD" "diamond OLD" "ganon OLD" "kraken OLD" "kraken2 OLD" \
       "centrifuge CG-top-3" "clark CG-top-3" "diamond CG-top-3" "ganon CG-top-3" "kraken CG-top-3" "kraken2 CG-top-3" \
       "diamond ALL-top-3" "ganon ALL-top-3" "kraken2 ALL-top-3" \
    -c "${blue_color}" "${blue_color}" "${blue_color}" "${blue_color}" "${blue_color}" "${blue_color}" \
       "${orange_color}" "${orange_color}" "${orange_color}" "${orange_color}" "${orange_color}" "${orange_color}" \
       "${red_color}" "${red_color}" "${red_color}" \
    -m "${centifuge_marker}" "${clark_marker}" "${diamond_marker}" "${ganon_marker}" "${kraken_marker}" "${kraken2_marker}" \
       "${centifuge_marker}" "${clark_marker}" "${diamond_marker}" "${ganon_marker}" "${kraken_marker}" "${kraken2_marker}" \
       "${diamond_marker}" "${ganon_marker}" "${kraken2_marker}" \
    -o all_${dt}_${ev}${suffix}.pdf > all_${dt}_${ev}${suffix}.log

  done
done
