# -*- coding: utf-8 -*-
import argparse
import sys
import numpy as np
from copy import deepcopy

vcf_file = sys.argv[1]
err_rate = 0.0001
contig_len = 10e6

positions = []
with open(vcf_file) as vcf:
    for line in vcf:
        if not line.startswith("#"):
            x_lin = line.split()
            positions.append(int(x_lin[1]))

mono = ['0|0' for i in x_lin[9:]]
x_lin[9:] = mono
nsamples = len(x_lin[9:])
errs = np.random.randint(1, contig_len, int(err_rate*contig_len))

f = open(vcf_file, 'a')
for i in errs:
    while i in positions:
        i += 1
    positions.append(i)
    null_line_r = deepcopy(x_lin)
    err_pos = np.random.randint(0, nsamples)
    null_line_r[err_pos] = '0|1'
    null_line_r[1] = str(i)
    new_err = "\t".join(null_line_r)
    f.write(f"{new_err}\n")
f.close()
# vcfstreamsort IN > OUT
