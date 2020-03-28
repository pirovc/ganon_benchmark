#!/usr/bin/env python3

import os, sys, pickle, gzip
from collections import defaultdict

info = pickle.load(gzip.open(sys.argv[1], 'rb'))

print('rank', info['rank'])
print('kmer_size', info['kmer_size'])
print('hash_functions', info['hash_functions'])
print('number_of_bins', info['number_of_bins'])
print()
print('bin_length', info['bin_length'])
print('fragment_length', info['fragment_length'])
print('overlap_length', info['overlap_length'])


cnt = defaultdict(set)
for line in info['bins']:
    if info['rank']=="assembly":
        acc, length, taxid, target, binid = line.split("\t")
    else:
        acc, length, taxid, binid = line.split("\t")
        target = taxid
        
    cnt[target].add(binid)
  
split_cnt = defaultdict(int)
for target,binids in cnt.items():
    split_cnt[len(binids)]+=1

print(split_cnt)
print("1", split_cnt[1])
split_cnt[1]=0
print(">1", sum(split_cnt.values()))
