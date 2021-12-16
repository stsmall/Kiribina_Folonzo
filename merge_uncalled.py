#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob
import gzip
import sys
from collections import defaultdict

contig = sys.argv[1]
input_files = glob.glob(f"*{contig}*")
missingdict = defaultdict(list)

for files in (input_files):
    with gzip.open(files, 'rb') as missing:
        for line in missing:
            line = line.decode()
            sample, chrom, *pos = line.split()
            missingdict[f"{sample}:{chrom}"].extend(pos)

with gzip.open(f"{contig}.merged.mask.out.gz", 'wt') as out:
    for sample in missingdict.keys():
        indv, chrom = sample.split(":")
        pos = set(missingdict[sample])  # sorted?
        pos_sorted = sorted(pos, key=int)
        pos_line = "\t".join(pos_sorted)
        out.write(f"{indv}\t{chrom}\t{pos_line}\n")
