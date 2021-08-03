# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 11:26:52 2021

@author: Scott T. Small

calculate the difference between two sets of dXY, FST statistics from pi_xy

Example
-------

    python diff_stats.py pixy_3_KFMoz.10kb_out_dxy.txt

"""
import sys
import numpy as np

infile = sys.argv[1]

with open(infile, 'r') as dxy:
    header = next(dxy)
    header = header.split()
    w1 = []
    w2 = []
    kf = []
    km = []
    fm = []
    for lin in dxy:
        x_lin = lin.split()
        p1 = x_lin[0]
        p2 = x_lin[1]
        w1.append(int(x_lin[3]))
        w2.append(int(x_lin[4]))
        try:
            d = float(x_lin[5])
        except ValueError:
            d = np.nan
        if p1 == "K" and (p2 == "F" or p2 == "Fs"):
            kf.append(d)
        elif p1 == "K" and p2 == "M":
            km.append(d)
        elif (p1 == "F" or p1 == "Fs") and p2 == "M":
            fm.append(d)
    w1 = np.sort(list(set(w1)))
    w2 = np.sort(list(set(w2)))


with open(f"{infile}.diff", 'w') as dxy:
    kf_arr = np.array(kf)
    km_arr = np.array(km)
    fm_arr = np.array(fm)
    dxy.write(header[:6])
    for i, w in enumerate(w1):
        dxy.write(f"KF\tKM\t3\t{w}\t{w2[i]}\t{kf[i]-km[i]}\n")
        dxy.write(f"KF\tFM\t3\t{w}\t{w2[i]}\t{kf[i]-fm[i]}\n")
