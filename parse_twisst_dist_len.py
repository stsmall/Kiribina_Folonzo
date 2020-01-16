#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 14:09:48 2020
@author: Scott T. Small

Finding the tree heights or coalescent times between species pairs is
best done using ETE3. This is simple for trees with 1 allele per species but
takes much more time for trees with >1 allele per species. The program twisst
can estimate distances and branch lengths when calculating weights.

twisst.py -t TREES_in -w WEIGHTS_out -D DISTANCE_out -L BRANCH_LEN_out
          --outputTopos TOPO_out -g FUN -g LIK -g VAN -g LON -g PAR -g RIV
          --groupsFile twisst.43.groups --method complete

This script will parse the Dist and Len files produced by twisst. Namely
specifying topologies or frequencies.


Example
-------
    python parse_twisst_dist_len.py -i 3L.blen.txt --blen -f 0.05 --bplot -c 3L
    python parse_twisst_dist_len.py -i X.blen.txt --blen -t topo3 topo10 topo11 --bplot
    python parse_twisst_dist_len.py -i 3R.dist.txt -p FUN-VAN -t topo3 topo10 --dist --chr 3R


"""
import sys
import argparse
import numpy as np
import pandas as pd
from os import path
from collections import defaultdict
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
mpl.rcParams['pdf.fonttype'] = 42


def boxplot_dist(chrm, div_df, pair, topo, topoplot, pairsplot):
    """plot the divergences
    """
    if pairsplot:
        label = pair
        if topo:
            df_list = div_df[topo].loc[pair].sort_index()
            x = [l for l in df_list]
            poplist = list(df_list.index.values)
        else:
            df_list = div_df.loc[pair].sort_index()
            x = [l for l in df_list]
            poplist = list(df_list.index.values)
    elif topoplot:
        label = topo
        x = [l for l in div_df[topo]]
        poplist = list(div_df[topo].index.values)
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.despine(ax=ax, offset=5)
    lw = 1.5
    box = ax.boxplot(
        x,
        labels=poplist,  patch_artist=True,
        medianprops={"color": "k", "linewidth": lw},
        whiskerprops={"color": "k"},
        capprops={"color": "k"},
        showfliers=False,
        flierprops={"c": "k", "markersize": 2})
    ax.set_ylabel(r'PwDiv', rotation=90, fontsize=12)
    ax.set_xticklabels(poplist, rotation=45, fontsize=8)
    # set list of colors
    colornames = list(mpl.colors.cnames.keys())[:len(poplist)]
    for patch, color in zip(box['boxes'], colornames):
        patch.set_facecolor(color)
        patch.set_linewidth(lw)
    fig.savefig("{}.{}.pdf".format(chrm, label), bbox_inches="tight")
    return(None)


def boxplot_branchlen(chrm, data, topos, toposfreq, minLen=1):
    """plots the nodeages
    """
    if topos:
        poplist = topos
    else:
        poplist = toposfreq
    topo_ix = [(int(x.strip("topo"))-1) for x in poplist]
    data_arr = np.array(data)[topo_ix]
    data_list = [a1[~np.isnan(a1)] for a1 in data_arr]
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.despine(ax=ax, offset=5)
    lw = 1.5
    box = ax.boxplot(
        data_list,
        labels=poplist,
        patch_artist=True,
        medianprops={"color": "k", "linewidth": lw},
        whiskerprops={"color": "k"},
        capprops={"color": "k"},
        showfliers=False,
        flierprops={"c": "k", "markersize": 2})
    ax.set_ylabel(r'NodeAge', rotation=90, fontsize=12)
    ax.set_xticklabels(poplist, rotation=45, fontsize=8)
    # set list of colors
    colornames = list(mpl.colors.cnames.keys())[:len(poplist)]
    for patch, color in zip(box['boxes'], colornames):
        patch.set_facecolor(color)
        patch.set_linewidth(lw)
    fig.savefig("{}.NodeDepth.pdf".format(chrm), bbox_inches="tight")
    return(None)


def parse_divergence(infile):
    """parses out distance file from twisst. returns a distribution for each
    species pair this should be applicable to other species

    Parameters
    ----------
    infile: str, file
        infile with distance data

    Returns
    -------
    div_df: obj
        Pandas DataFrame

    """
    div_dict = defaultdict(lambda: defaultdict(list))
    with open(infile, 'r') as f:
        header = next(f).split()
        for line in f:
            for i, div in enumerate(line.split()):
                topo, ind1, ind2 = header[i].split("_")
                pairname = "{}-{}".format(ind1, ind2)
                if np.isnan(div):
                    pass
                else:
                    div_dict[topo][pairname].append(float(div))
    div_df = pd.DataFrame(div_dict)

    return(div_df)


def calc_mean_distance(div_df, topo, pair):
    """Calculates mean and median from

    Parameters
    ----------
    div_df: obj
        pandas dataframe
    topo: str
        which topology
    pair:str
        which pair

    Returns
    -------
    pmean: flt
        numpy mean
    pmedian: flt
        numpy median
    quant_U: flt
        numpy upper quantile 97.5
    qunat_L: flt
        numpy lower quantile 2.5

    """
    pmean = np.nanmean(div_df[topo].loc[pair])
    pmedian = np.nanmedian(div_df[topo].loc[pair])
    quant_U = np.nanpercentile(div_df[topo].loc[pair], 97.5)
    quant_L = np.nanpercentile(div_df[topo].loc[pair], 2.5)

    return pmean, pmedian, quant_U, quant_L


def mean_distances(chrm, div_df, topos, pairs):
    """
    Parameters
    ----------
    chrm: str
        chromosome name
    div_df: dataframe
        dataframe of divergence data
    topos: List(str)
        list of topos for distance calculations
    pairs: List(str)
        list of pairs for distance calculations

    Returns
    -------
    NONE

    """
    # calculate mean distances
    if path.exists("{}.meandist.out".format(chrm)):
        f = open("{}.meandist.out".format(chrm), 'a')
    else:
        f = open("{}.meandist.out".format(chrm), 'w')
    if topos:  # specific topos
        for topo in topos:
            if pairs:  # specific pairs
                for pair in pairs:
                    pmean, pmedian, quant_U, quant_L = calc_mean_distance(div_df, topo, pair)
                    f.write(f"{topo}:{pair}:{pmean} {pmedian} [{quant_U}-{quant_L}]\n")
            else:
                for pair in list(div_df.index.values):
                    pmean, pmedian, quant_U, quant_L = calc_mean_distance(div_df, topo, pair)
                    f.write(f"{topo}:{pair}:{pmean} {pmedian} [{quant_U}-{quant_L}]\n")
    else:
        for topo in list(div_df.columns.values):
            if pairs:  # specific pairs
                for pair in pairs:
                    pmean, pmedian, quant_U, quant_L = calc_mean_distance(div_df, topo, pair)
                    f.write(f"{topo}:{pair}:{pmean} {pmedian} [{quant_U}-{quant_L}]\n")
            else:
                for piar in list(div_df.index.values):
                    pmean, pmedian, quant_U, quant_L = calc_mean_distance(div_df, topo, pair)
                    f.write(f"{topo}:{pair}:{pmean} {pmedian} [{quant_U}-{quant_L}]\n")
    f.close()


def calc_mrca(chrm, topos, blen_data, min_freq, tree_count):
    """Finds the MRCA for each topos

    Parameters
    ----------
    chrm: str
        chromosome or how to name the outfile
    blen_data: List(List)
        zipped list of MRCA from topologies
    min_freq: flt
        min frequency cut off to keep a topo for plotting

    Returns
    -------
    topos_freq: list
        topos passing the min_freq cutoff

    """
    topos_freq = []
    t = open(f"{chrm}.fulldata.out", "w")
    with open(f"{chrm}.nodedepth.out", 'w') as f:
        for i, mrca in enumerate(blen_data):
            topo = f"topo{i + 1}"
            nancount = sum(np.isnan(mrca))
            freq = (1 - (nancount/tree_count))
            if freq >= min_freq:
                topos_freq.append(f"{topo}")  # pass min_freq
                if (topos is None) or (topo in topos):
                    # print full data
                    for blen in mrca:
                        if np.isnan(blen):
                            pass
                        else:
                            t.write(f"{chrm}\t{topo}\t{blen}\n")
                    mean = np.nanmean(mrca)
                    median = np.nanmedian(mrca)
                    quant_low = np.nanpercentile(mrca, 2.5)
                    quant_up = np.nanpercentile(mrca, 97.5)
                    f.write(f"{topo} {freq:.3f}:{mean:.3f} {median:.3f} [{quant_low:.3f}-{quant_up:.3f}]\n")
    return topos_freq


def sum_branch_lengths(chrm, infile, min_freq, topos, step=10, outgroup_pos=2):
    """calculates the nodeage of the MRCA using a branchlength file produced by
    twisst

    to alter this for other species you must:
        1) change step size
        2) change outgroup blens

    The blen file is

    topo1_col1 topo1_col2 topo1_col3 topo2_col1 topo2_col2 topo2_col3
    val val val val val val val val val val val

    Parameters
    ----------
    outgroup_pos: int
        skip the first 2 as they are outgroup root lengths and we want ingroup
    step: int
        there are 8 node times and 2 outgroup time ... this is 10
    topos: List(str)
        list of topos for distance calculations

    Returns
    -------
    blen_data: List(List)
        zipped list of MRCA from topologies
    topos_freq: list
        topos passing the min_freq cutoff

    """
    blen_boxplot = []
    tree_count = 0  # total window count
    topodict = {}
    with open(infile, 'r') as blen_file:
        for line in blen_file:
            if line.startswith("#"):
                x = line.split()
                topo = x[0][1:]
                leaf_names = [leaf.split("--")[1] for leaf in x[3:]]
                topodict[topo] = leaf_names
            else:
                tree_count += 1
                start = 0
                stop = step
                blen = list(map(float, line.split()))
                blen_list = []
                while stop <= len(blen):
                    blen_vals = blen[start + outgroup_pos:stop]
                    if np.isnan(blen_vals[0]):
                        blen_list.append(np.nan)
                    else:
                        topo_key = start//step  # topos are column header
                        blen_topo = topodict[f"topo{topo_key+1}"]
                        leaf_ix = [blen_topo.index(bt) for bt in blen_topo if bt.count('_') == 0]
                        # should return 3 or more
                        blen_max = max([blen_vals[bl] for bl in leaf_ix])
                        # get all values for leaf index and max
                        blen_maxix = blen_topo[blen_vals.index(blen_max)]
                        #  here is the leaf name 'Fun'
                        leaf_dist = [blen_vals[blen_topo.index(x1)] for x1 in blen_topo if blen_maxix in x1]
                        # now return all vals with instances of 'Fun' in topo
                        leaf_sum = sum(leaf_dist)
                        blen_list.append(leaf_sum)
                    start += step
                    stop += step
                blen_boxplot.append(blen_list)
    # boxplot of node depth
    blen_data = list(zip(*blen_boxplot))
    topos_freq = calc_mrca(chrm, topos, blen_data, min_freq, tree_count)
    return blen_data, topos_freq


def parse_args(args_in):
    parser = argparse.ArgumentParser(prog="sys.argv[0].py",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', "--infile", type=str, required=True,
                        help="infile from distmatrix twisst")
    parser.add_argument('-t', "--topos", nargs='+', required=False,
                        help="list of topos e.g., topo82")
    parser.add_argument('-p', "--pairs", nargs='+', required=False,
                        help="list of interested pairs, e.g. Fun-Lik")
    parser.add_argument('-f', "--minFreq", type=float, default=0.0)
    parser.add_argument("--dist", action="store_true",
                        help="a divergence file")
    parser.add_argument("--blen", action="store_true",
                        help="a file of branch"
                        "lengths")
    parser.add_argument('-c', "--chr", type=str, default='chr',
                        help="prefix of outfile")
    parser.add_argument("--tplot", action="store_true",
                        help="  # for each topo returns all pairwise distances"
                        " on 1 plot")
    parser.add_argument("--pplot", action="store_true",
                        help="for each pair returns a boxplot of distances on"
                        " each topo")
    parser.add_argument("--bplot", action="store_true",
                        help="boxplot of nodeage for specified topos")
    return(parser.parse_args(args_in))


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    # =========================================================================
    #  Gather args
    # =========================================================================
    infile = args.infile
    topos = args.topos
    pairs = args.pairs
    dist = args.dist
    blen = args.blen
    chrom = args.chr
    min_freq = args.minFreq
    tplot = args.tplot
    pplot = args.pplot
    bplot = args.bplot
    # =========================================================================
    #  Main executions
    # =========================================================================
    if dist:
        div_df = parse_divergence(infile)
        # make boxplot
        if tplot:
            for topo in topos:
                boxplot_dist(chrom, div_df, pairs, topos)
        if pplot:
            for pair in pairs:
                boxplot_dist(chrom, div_df, pairs, topos)
    elif blen:
        blen_data, topos_freq = sum_branch_lengths(chrom, infile, min_freq, topos)
        if bplot:
            boxplot_branchlen(chrom, blen_data, topos, topos_freq)
