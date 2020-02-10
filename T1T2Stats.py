#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from Bio import Phylo as phy
import numpy as np
from scipy import stats


def distanceAB_T1T2(tree, AB_T1, AB_T2):
    """
    """
    distanceAB = tree.distance(sp2, sp3)  # AB
    distanceAC = tree.distance(sp2, sp4)  # (B),C
    distanceBC = tree.distance(sp3, sp4)  # (A),C
    AB_T2.append(distanceAB)
    AB_T1.append((distanceAC + distanceBC)/2)
    return(AB_T1, AB_T2)


def distanceBC_T1T2(tree, BC_T1, BC_T2):
    """
    """
    distanceBC = tree.distance(sp3, sp4)  # BC
    distanceBA = tree.distance(sp4, sp2)  # (C),A
    distanceCA = tree.distance(sp3, sp2)  # (B),A
    BC_T2.append(distanceBC)
    BC_T1.append((distanceCA + distanceBA/2))
    return(BC_T1, BC_T2)


def distanceAC_T1T2(tree, AC_T1, AC_T2):
    """
    """
    distanceAC = tree.distance(sp2, sp4)  # AC
    distanceAB = tree.distance(sp2, sp3)  # (A),B
    distanceCB = tree.distance(sp4, sp3)  # (C),B
    AC_T2.append(distanceAC)
    AC_T1.append((distanceAB + distanceCB)/2)
    return(AC_T1, AC_T2)


def t1t2_stats_out(AB_T1, AB_T2, BC_T1, BC_T2, AC_T1, AC_T2):
    """
    """
    # T2
    sp23T2 = np.mean(AB_T2)
    sp34T2 = np.mean(BC_T2)
    sp24T2 = np.mean(AC_T2)
    # MWU test
#    T2 = [sp23T2, sp34T2, sp24T2]
#    tix = T2.index(min(T2))
#    T2_list = [AB_T2, BC_T2, AC_T2]
#    xT2u, xT2p = stats.mannwhitneyu(AB_T2, BC_T2, alternative='greater')
#    yT2u, yT2p = stats.mannwhitneyu(AB_T2, AC_T2, alternative='greater')
#    zT2u, zT2p = stats.mannwhitneyu(BC_T2, AC_T2, alternative='greater')

    # T1
    sp234T1 = np.mean(AB_T1)
    sp342T1 = np.mean(BC_T1)
    sp423T1 = np.mean(AC_T1)
#    # MWU test
#    xT1u, xT1p = stats.mannwhitneyu(AB_T1, BC_T1, alternative='greater')
#    yT1u, yT1p = stats.mannwhitneyu(AB_T1, AC_T1, alternative='greater')
#    zT1u, zT1p = stats.mannwhitneyu(BC_T1, AC_T1, alternative='greater')

    # internal node length
    sp2C2 = sp234T1 - sp23T2
    sp3C2 = sp342T1 - sp34T2
    sp4C2 = sp423T1 - sp24T2

    print(f"({sp2},{sp3}),{sp4}  T2:{sp23T2} T1:{sp234T1} C2:{sp2C2}")
    print(f"({sp3},{sp4}),{sp2}  T2:{sp34T2} T1:{sp342T1} C2:{sp3C2}")
    print(f"({sp4},{sp2}),{sp3}  T2:{sp24T2} T1:{sp423T1} C2:{sp4C2}")

#    print(f"{sp2},{sp3} gt {sp3},{sp4} = pvalue {xT2p}")
#    print(f"{sp2},{sp3} gt {sp2},{sp4} = pvalue {yT2p}")
#    print(f"{sp3},{sp4} gt {sp2},{sp4} = pvalue {zT2p}")
#
#    print(f"({sp2},{sp3}),{sp4} gt ({sp3},{sp4}),{sp2} = pvalue {xT1p}")
#    print(f"({sp2},{sp3}),{sp4} gt ({sp2},{sp4}),{sp3} = pvalue {yT1p}")
#    print(f"({sp3},{sp4}),{sp2} gt ({sp2},{sp4}),{sp3} = pvalue {zT1p}")


if __name__ == "__main__":
    # =========================================================================
    #  Gather args
    # =========================================================================
    tree_file = sys.argv[1]
    sp2 = sys.argv[2]  # A
    sp3 = sys.argv[3]  # B
    sp4 = sys.argv[4]  # C
    # =========================================================================
    #  Main executions
    # =========================================================================
    BC_T2 = []
    BC_T1 = []
    c_BC = 0
    AC_T1 = []
    AC_T2 = []
    c_AC = 0
    AB_T1 = []
    AB_T2 = []
    c_AB = 0
    trees = phy.parse(tree_file, 'newick')
    for tree in trees:  # for each genealogy
        if (tree.distance(sp3, sp4) < tree.distance(sp2, sp3)) and (tree.distance(sp3, sp4) < tree.distance(sp2, sp4)):  # (3,4),2
            c_BC += 1
            BC_T1, BC_T2 = distanceBC_T1T2(tree, BC_T1, BC_T2)
        elif (tree.distance(sp2, sp3) < tree.distance(sp3, sp4)) and (tree.distance(sp2, sp3) < tree.distance(sp2, sp4)):  # (2,3),4
            c_AB += 1
            AB_T1, AB_T2 = distanceAB_T1T2(tree, AB_T1, AB_T2)
        elif (tree.distance(sp2, sp4) < tree.distance(sp2, sp3)) and (tree.distance(sp2, sp4) < tree.distance(sp3, sp4)):  # (2,4),3
            c_AC += 1
            AC_T1, AC_T2 = distanceAC_T1T2(tree, AC_T1, AC_T2)
    t1t2_stats_out(AB_T1, AB_T2, BC_T1, BC_T2, AC_T1, AC_T2)
