#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 17:01:00 2020
@author: Scott T. Small

This module is a direct implementation of:
    https://github.com/tskit-dev/tskit/issues/665

Example
-------

    $ python windowed_gnn.py TREE --ref 0 1 --foc 2

Notes
-----


"""
import sys
import argparse
import tskit
import json
import pandas as pd
import numpy as np
from os import path
# plotting
import matplotlib as mpl
import matplotlib.pyplot as plt
# pdf edits
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['pdf.fonttype'] = 42


COLOURS = ["#FDBF6F", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
           "#E31A1C", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928",
           "#d442f5", "#403f40", "#291dad"]

FOCAL_IND = "02-06207"  # sets a grey box highlight of hap in GNN plot


def plot_gnn_wg(gnndf, groups, focal_ind):
    """Plot gnn.

    Parameters
    ----------
    gnndf : TYPE
        DESCRIPTION.
    groups : TYPE
        DESCRIPTION.
    focal_inds : TYPE, optional
        DESCRIPTION. The default is "02-06207".

    Returns
    -------
    None.

    """
    # load df
    df = pd.read_csv("gnndf")

    # prepare data
    A = np.zeros((2, len(df)))
    for j, pop in enumerate(groups):
        A[j, :] = df[pop].values
    index = np.argsort(A[0])[::-1]
    inds = df.Isolate.values
    inds = list(inds[index])
    x1 = inds.index(focal_ind)
    x2 = inds.index(focal_ind, x1 + 1)
    A = A[:, index]

    # plotting
    colours = {g: COLOURS[i] for i, g in enumerate(groups)}
    fig, ax = plt.figure(figsize=(14, 4))
    x = np.arange(len(df))
    for j, region in enumerate(groups):
        ax.bar(x, A[j], bottom=np.sum(A[:j, :], axis=0), label=region, width=1,
               color=colours[region], align="edge")
    ax.set_xlim(0, len(df) - 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_ylabel("GNN Fraction")
    ax.set_title("Haplotypes sorted by GNN")

    # add focal outlier
    if focal_ind:
        for x in [x1, x2]:
            p = mpl.patches.Rectangle(
                (x, 0), width=1, height=1, fill=False, linestyle="--", color="grey")
            ax.add_patch(p)

    ax.legend(bbox_to_anchor=(1.02, 0.76))
    fig.savefig(f"{gnndf}.pdf", bbox_inches='tight')


def gnn_fx(outfile, ts, ref_samples, focal, groups, savedf=True):
    """Run mean GNN fx.

    Parameters
    ----------
    outfile : TYPE
        DESCRIPTION.
    ts : TYPE
        DESCRIPTION.
    ref_set : TYPE
        DESCRIPTION.
    focal : TYPE
        DESCRIPTION.
    groups : TYPE
        DESCRIPTION.
    savedf : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    None.

    """
    if not focal:
        focal = ts.samples()  # all samples

    # calc gnn
    gnn = ts.genealogical_nearest_neighbours(focal, ref_samples)

    # save results
    if savedf:
        sample_nodes = [ts.node(n) for n in ts.samples()]
        sample_ids = [n.id for n in sample_nodes]
        sample_names = [json.loads(ts.individual(n.individual).metadata)['Isolate'] for n in sample_nodes]
        sample_pops = [json.loads(ts.population(n.population).metadata)['Group'] for n in sample_nodes]
        gnn_table = pd.DataFrame(
            data=gnn,
            index=[pd.Index(sample_ids, name="Sample_node"),
                   pd.Index(sample_names, name="Isolate"),
                   pd.Index(sample_pops, name="Group")],
            columns=groups
            )

        gnn_table.to_csv(f"GNN.{outfile}.csv")


def parse_args(args_in):
    """Parse args."""
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("tree", type=str, help="tskit tree object")
    parser.add_argument("--ref", type=str, args="*", action="append",
                        help="reference sets to compare to focal node")
    parser.add_argument("--foc", type=str, args="*", action="append",
                        help="focal nodes")
    parser.add_argument("--gnn_windows", action="store_true",
                        help="run gnn in windows mode")
    parser.add_argument("--time_windows", action="store_true",
                        help="run gnn in time windows mode")
    return parser.parse_args(args_in)


def main():
    """Run main function."""
    args = parse_args(sys.argv[1:])
    # =========================================================================
    #  Gather args
    # =========================================================================
    tree = args.tree
    outfile = path.split(tree)[-1]
    ref_set = args.ref
    foc_set = args.foc
    tree_windows = args.gnn_windows
    time_windows = args.time_windows
    # =========================================================================
    #  Main executions
    # =========================================================================
    ts = tskit.load(tree)
    if not ref_set:
        ref_set = range(ts.num_populations)  # all populations
    groups = [json.loads(ts.population(i).metadata)["Group"] for i in ref_set]
    ref_samples = [ts.samples(population=pop_id) for pop_id in ref_set]

    if not tree_windows and not time_windows:
        gnn_fx(outfile, ts, ref_samples, foc_set, groups)
        plot_gnn_wg(f"GNN.{outfile}.csv", groups, FOCAL_IND)
    else:
        pass


if __name__ == "__main__":
    main()
