# -*- coding: utf-8 -*-
"""
Created on Fri Jul 30 11:22:15 2021

@author: Scott T. Small

"""

import msprime
import numpy as np
from hejase_stats import calc_cc10
from hejase_stats import calc_tmrcah

p_nodes1 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
p_nodes2 = [16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31]
p_nodes_cc = [p_nodes1, p_nodes2]
ts = msprime.sim_ancestry(16, random_seed=234)

def test_tmrcah(ts, p_nodes1):
    mid, tmrcah_rel, time_rel, time_rel2 = calc_tmrcah(ts, p_nodes1)
    assert np.around(mid, 2)[0] == 0.5
    assert np.around(tmrcah_rel, 2)[0] == 0.48
    assert np.around(time_rel, 2)[0] == 1.84
    assert np.around(time_rel2, 2)[0] == 1.35

def test_cc10(ts, p_nodes_cc):
    mid, cc10_ls, time_rel = calc_cc10(ts, p_nodes_cc)
    assert np.around(mid, 2)[0] == 0.5
    assert np.around(cc10_ls, 2) == [0.02, 0.06, 0.08, 0.09, 0.14, 0.17, 0.31, 0.31, 0.37, 0.51]
    assert np.around(time_rel, 2)[0] == 1.35