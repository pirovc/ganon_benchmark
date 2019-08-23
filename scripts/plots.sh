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
    -k 'phylum' 'class' 'order' 'family' 'genus' 'species' 'species+' \
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


    # compare top 3 and complete
    ganon_benchmark/scripts/plots.py \
    -i ganon_eval/new_results/20181219_abfv_refseq_cg/cami_${dt}_H01${suffix}/*.${ev}.npz \
       ganon_eval/new_results/20181219_abfv_refseq_cg_t3a/cami_${dt}_H01${suffix}/*.${ev}.npz \
    -k 'phylum' 'class' 'order' 'family' 'genus' 'species' 'species+' \
    -l "centrifuge CG" "clark CG" "diamond CG" "ganon CG" "kraken CG" \
       "centrifuge CG top 3" "clark CG top 3" "diamond CG top 3" "ganon CG top 3" "kraken CG top 3" \
    -c "${blue_color}" "${blue_color}" "${blue_color}" "${blue_color}" "${blue_color}" \
       "${orange_color}" "${orange_color}" "${orange_color}" "${orange_color}" "${orange_color}" \
    -m "${centifuge_marker}" "${clark_marker}" "${diamond_marker}" "${ganon_marker}" "${kraken_marker}" \
       "${centifuge_marker}" "${clark_marker}" "${diamond_marker}" "${ganon_marker}" "${kraken_marker}" \
    -o complete_vs_top_cg_${dt}_${ev}${suffix}.png > complete_vs_top_cg_${dt}_${ev}${suffix}.log
    ganon_benchmark/scripts/plots.py \
    -i ganon_eval/new_results/20181219_abfv_refseq_all/cami_${dt}_H01${suffix}/*.${ev}.npz \
       ganon_eval/new_results/20181219_abfv_refseq_all_t3a/cami_${dt}_H01${suffix}/*.${ev}.npz \
    -k 'phylum' 'class' 'order' 'family' 'genus' 'species' 'species+' \
    -l "diamond ALL" "ganon ALL" \
       "diamond ALL top 3" "ganon ALL top 3"  \
    -c "${blue_color}" "${blue_color}" \
       "${orange_color}" "${orange_color}" \
    -m "${diamond_marker}" "${ganon_marker}" \
       "${diamond_marker}" "${ganon_marker}" \
    -o complete_vs_top_all_${dt}_${ev}${suffix}.png > complete_vs_top_all_${dt}_${ev}${suffix}.log

  done
done

# compare assembly and taxid modes
ganon_benchmark/scripts/plots.py \
  -i ganon_eval/new_results/20150602_bacteria_refseq_old/cami_toy_H01-assembly/*.rank.npz \
     ganon_eval/new_results/20150602_bacteria_refseq_old/cami_toy_H01/{centrifuge,ganon,kraken}*.rank.npz \
     ganon_eval/new_results/20181219_abfv_refseq_cg/cami_toy_H01-assembly/*.rank.npz \
     ganon_eval/new_results/20181219_abfv_refseq_cg_t3a/cami_toy_H01/{centrifuge,ganon,kraken}*.rank.npz \
     ganon_eval/new_results/20181219_abfv_refseq_all/cami_toy_H01-assembly/*.rank.npz \
     ganon_eval/new_results/20181219_abfv_refseq_all_t3a/cami_toy_H01/ganon*.rank.npz \
  -k 'phylum' 'class' 'order' 'family' 'genus' 'species' 'species+' 'assembly' \
  -l "centrifuge OLD (assembly)" "ganon OLD (assembly)" "krakenuniq OLD (assembly)" \
     "centrifuge OLD (taxid)" "ganon OLD (taxid)" "kraken OLD (taxid)" \
     "centrifuge CG (assembly)" "ganon CG (assembly)" "krakenuniq CG (assembly)" \
     "centrifuge CG top 3 (taxid)" "ganon CG top 3 (taxid)" "kraken CG top 3 (taxid)" \
     "ganon ALL (assembly)" \
     "ganon ALL top 3 (taxid)" \
  -c "${blue_color}" "${blue_color}" "${blue_color}" \
     "${blue_color}" "${blue_color}" "${blue_color}" \
     "${orange_color}" "${orange_color}" "${orange_color}" \
     "${orange_color}" "${orange_color}" "${orange_color}" \
     "${red_color}" \
     "${red_color}" \
  -m "${centifuge_marker}" "${ganon_marker}" "${kraken_marker}" \
     "${centifuge_marker}" "${ganon_marker}" "${kraken_marker}" \
     "${centifuge_marker}" "${ganon_marker}" "${kraken_marker}" \
     "${centifuge_marker}" "${ganon_marker}" "${kraken_marker}" \
     "${ganon_marker}" \
     "${ganon_marker}" \
  -n "dotted" "dotted" "dotted" \
     "solid" "solid" "solid" \
     "dotted" "dotted" "dotted" \
     "solid" "solid" "solid" \
     "dotted" \
     "solid" \
  -o assembly_taxid_${dt}_${ev}${suffix}.png > assembly_taxid_${dt}_${ev}${suffix}.log
