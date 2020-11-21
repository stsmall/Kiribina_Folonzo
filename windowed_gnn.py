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


def gnn_fx(outfile, ts, ref_samples, target_samples, ref_groups, savedf=True):
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
    # calc gnn
    gnn = ts.genealogical_nearest_neighbours(target_samples, ref_samples)

    # save results
    if savedf:
        #sample_nodes = [ts.node(n) for n in ts.samples()]
        sample_nodes = [ts.node(n) for n in target_samples]
        sample_ids = [n.id for n in sample_nodes]
        sample_names = [json.loads(ts.individual(n.individual).metadata)['Isolate'] for n in sample_nodes]
        sample_pops = [json.loads(ts.population(n.population).metadata)['Group'] for n in sample_nodes]
        gnn_table = pd.DataFrame(
            data=gnn,
            index=[pd.Index(sample_ids, name="Sample_node"),
                   pd.Index(sample_names, name="Isolate"),
                   pd.Index(sample_pops, name="Group")],
            columns=ref_groups
            )

        gnn_table.to_csv(f"GNN.{outfile}.csv")


def parse_time_windows(ts, time_windows):
    """Parse time windows.

    Parameters
    ----------
    ts : TYPE
        DESCRIPTION.
    time_windows : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    if time_windows is None:
        time_windows = [0.0, ts.max_root_time]
    return np.array(time_windows)


def windowed_gnn(ts,
                 focal,
                 reference_sets,
                 windows=None,
                 time_windows=None,
                 span_normalise=True,
                 time_normalise=True):
    """Run GNN in windows.

    Genealogical_nearest_neighbours with support for span- and time-based windows.

    Parameters
    ----------
    ts : Obj
        tskit tree object
    focal : list
        focal nodes
    reference_sets : list
        reference nodes
    windows : TYPE, optional
        DESCRIPTION. The default is None.
    time_windows : TYPE, optional
        DESCRIPTION. The default is None.
    span_normalise : TYPE, optional
        DESCRIPTION. The default is True.
    time_normalise : TYPE, optional
        DESCRIPTION. The default is True.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    A : TYPE
        DESCRIPTION.

    """
    reference_set_map = np.zeros(ts.num_nodes, dtype=int) - 1
    for k, reference_set in enumerate(reference_sets):
        for u in reference_set:
            if reference_set_map[u] != -1:
                raise ValueError("Duplicate value in reference sets")
            reference_set_map[u] = k

    windows_parsed = ts.parse_windows(windows)
    num_windows = windows_parsed.shape[0] - 1
    time_windows_parsed = parse_time_windows(ts, time_windows)
    num_time_windows = time_windows_parsed.shape[0] - 1
    A = np.zeros((num_windows, num_time_windows, len(focal), len(reference_sets)))
    K = len(reference_sets)
    parent = np.zeros(ts.num_nodes, dtype=int) - 1
    sample_count = np.zeros((ts.num_nodes, K), dtype=int)
    time = ts.tables.nodes.time
    norm = np.zeros((num_windows, num_time_windows, len(focal)))

    # Set the initial conditions.
    for j in range(K):
        sample_count[reference_sets[j], j] = 1

    window_index = 0
    for (t_left, t_right), edges_out, edges_in in ts.edge_diffs():
        for edge in edges_out:
            parent[edge.child] = -1
            v = edge.parent
            while v != -1:
                sample_count[v] -= sample_count[edge.child]
                v = parent[v]
        for edge in edges_in:
            parent[edge.child] = edge.parent
            v = edge.parent
            while v != -1:
                sample_count[v] += sample_count[edge.child]
                v = parent[v]

        # Update the windows
        assert window_index < num_windows
        while (windows_parsed[window_index] < t_right and window_index + 1 <= num_windows):
            w_left = windows_parsed[window_index]
            w_right = windows_parsed[window_index + 1]
            left = max(t_left, w_left)
            right = min(t_right, w_right)
            span = right - left
            # Process this tree.
            for j, u in enumerate(focal):
                focal_reference_set = reference_set_map[u]
                delta = int(focal_reference_set != -1)
                p = u
                while p != tskit.NULL:
                    total = np.sum(sample_count[p])
                    if total > delta:
                        break
                    p = parent[p]
                if p != tskit.NULL:
                    scale = span / (total - delta)
                    time_index = np.searchsorted(time_windows_parsed, time[p]) - 1
                    if time_index < num_time_windows and time_index >= 0:
                        for k, _reference_set in enumerate(reference_sets):
                            n = sample_count[p, k] - int(focal_reference_set == k)
                            A[window_index, time_index, j, k] += n * scale
                        norm[window_index, time_index, j] += span
            assert span > 0
            if w_right <= t_right:
                window_index += 1
            else:
                # This interval crosses a tree boundary, so we update it again
                # in the for the next tree
                break

    # Reshape norm depending on normalization selected
    # Return NaN when normalisation value is 0
    if span_normalise and time_normalise:
        # norm[norm == 0] = 1
        A /= norm.reshape((num_windows, num_time_windows, len(focal), 1))
    elif span_normalise and not time_normalise:
        norm = np.sum(norm, axis=1)
        # norm[norm == 0] = 1
        A /= norm.reshape((num_windows, 1, len(focal), 1))
    elif time_normalise and not span_normalise:
        norm = np.sum(norm, axis=0)
        # norm[norm == 0] = 1
        A /= norm.reshape((1, num_time_windows, len(focal), 1))

    A[np.all(A == 0, axis=3)] = np.nan

    # Remove dimension for windows and/or time_windows if parameter is None
    if windows is None and time_windows is not None:
        A = A.reshape((num_time_windows, len(focal), len(reference_sets)))
    elif time_windows is None and windows is not None:
        A = A.reshape((num_windows, len(focal), len(reference_sets)))
    elif time_windows is None and windows is None:
        A = A.reshape((len(focal), len(reference_sets)))
    return A


def gnn_windows_fx(outfile, ts, ref_samples, target_samples, ref_groups, foc,
                   gnn_win, gnn_time, median=True, savedf=True):
    """Calculate gnn in windows.

    Parameters
    ----------
    ts : TYPE
        DESCRIPTION.
    ref_samples : TYPE
        DESCRIPTION.
    foc_set : TYPE
        DESCRIPTION.
    groups : TYPE
        DESCRIPTION.

    Returns
    -------
    gnn_dict : TYPE
        DESCRIPTION.

    """
    windows = list(ts.breakpoints())  # all trees
    if gnn_win and gnn_time:
        gnn = windowed_gnn(ts, target_samples, ref_samples, windows=windows, time_windows=[0, 50])
    elif gnn_win and not gnn_time:
        gnn = windowed_gnn(ts, target_samples, ref_samples, windows=windows)
    elif gnn_time and not gnn_win:
        gnn = windowed_gnn(ts, target_samples, ref_samples, time_windows=[0, 50])
    if median:
        gnn_m = np.median(gnn, axis=1)
    else:
        gnn_m = np.mean(gnn, axis=1)

    if savedf:  # save to df
        left = list(ts.breakpoints())[:-1]
        right = list(ts.breakpoints())[1:]
        group = [foc] * len(left)
        gnn_table = pd.DataFrame(
            data=gnn_m,
            index=[pd.Index(left, name="left_bp"),
                   pd.Index(right, name="right_bp"),
                   pd.Index(group, name="Group")],
            columns=ref_groups
            )
        gnn_table.to_csv(f"GNN_windows.{outfile}.{foc}.csv")

    return gnn_m


def plot_gnn_windows_p1(outfile, ts, gnn_m, groups, foc, pfix=0.95, mark_outliers=True):
    """Plot output from gnn windows.

    Parameters
    ----------
    outfile : TYPE
        DESCRIPTION.
    ts : TYPE
        DESCRIPTION.
    gnn_dict : TYPE
        DESCRIPTION.
    groups : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    left = list(ts.breakpoints())[:-1]
    right = list(ts.breakpoints())[1:]
    width = np.subtract(right, left)
    total = np.zeros_like(width)
    colours = {g: COLOURS[i] for i, g in enumerate(groups)}

    A = gnn_m
    # A_norm = np.sum(A, axis=1)

    # plotting 1
    fig, ax = plt.subplots(1, figsize=(14, 4))
    for j, pop in enumerate(groups):
        ax.bar(left, A[:, j], bottom=total, width=width, align="edge",
               label=pop, color=colours[pop])
        total += A[:, j]
    ax.set_title(f"Chromosome painting ({foc})")
    ax.set_xticks([np.linspace(right[0], right[-1], 10)])
    ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    ax.set_xlim(0, np.max(right))
    ax.set_ylim(0, 1)
    ax.legend(bbox_to_anchor=(1.02, 0.76))
    ax.set_ylabel("GNN Fraction")
    ax.set_xlabel("Position")
    # mark outlier windows
    k_out = A[:, 0] > pfix
    left_arr = np.array(left)
    width_arr = np.array(width)
    k_left = left_arr[k_out]
    k_width = width_arr[k_out]
    if mark_outliers:
        for l, w in zip(k_left, k_width):
            p = mpl.patches.Rectangle(
                (l, 0), width=w, height=1, fill=False, linestyle="--", color="grey")
            ax.add_patch(p)
    fig.savefig(f"GNN_windows.{outfile}.{foc}.p1.pdf", bbox_inches='tight')




def plot_gnn_windows_p2(outfile, ts, gnn_m, groups, foc, pfix=0.95, mark_outliers=True):
    """Plot output from gnn windows.

    Parameters
    ----------
    outfile : TYPE
        DESCRIPTION.
    ts : TYPE
        DESCRIPTION.
    gnn_dict : TYPE
        DESCRIPTION.
    groups : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    left = list(ts.breakpoints())[:-1]
    right = list(ts.breakpoints())[1:]
    width = np.subtract(right, left)
    total = np.zeros_like(width)
    colours = {g: COLOURS[i] for i, g in enumerate(groups)}

    A = gnn_m
    # A_norm = np.sum(A, axis=1)

    # plotting 2
    gs = mpl.gridspec.GridSpec(2, 1, hspace=0.6)
    fig = plt.figure(figsize=(14, 4))
    ax_left = plt.subplot(gs[0])
    ax_right = plt.subplot(gs[1])
    for i, ax in enumerate([ax_left, ax_right]):
        group = groups[i]
        ax.bar(left, A[:, i], bottom=total, width=width, align="edge", label=group, color=colours[group])
        ax.set_title(f"Chromosome painting ({group})")
        ax.set_xticks([np.linspace(right[0], right[-1], 10)])
        ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
        ax.set_xlim(0, np.max(right))
        ax.set_ylim(0, 1)
        ax.legend(bbox_to_anchor=(1.02, 0.76))
        ax.set_ylabel("GNN Fraction")
        ax.set_xlabel("Position")
        # mark outlier windows
        k_out = A[:, i] > pfix
        left_arr = np.array(left)
        width_arr = np.array(width)
        k_left = left_arr[k_out]
        k_width = width_arr[k_out]
        if mark_outliers:
            for l, w in zip(k_left, k_width):
                p = mpl.patches.Rectangle(
                    (l, 0), width=w, height=1, fill=False, linestyle="--", color="grey")
                ax.add_patch(p)
    fig.savefig(f"GNN_windows.{outfile}.{foc}.p2.pdf", bbox_inches='tight')


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
    df = pd.read_csv(gnndf)
    # prepare data
    A = np.zeros((len(groups), len(df)))
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
    fig, ax = plt.subplots(1, figsize=(14, 4))
    x = np.arange(len(df))
    for j, region in enumerate(groups):
        ax.bar(x, A[j], bottom=np.sum(A[:j, :], axis=0), label=region, width=1,
               color=colours[region], align="edge")
    ax.set_xlim(0, len(df) - 1)
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
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


def parse_args(args_in):
    """Parse args."""
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("tree", type=str, help="tskit tree object")
    parser.add_argument("--ref", type=str, nargs="*", action="append",
                        help="reference sets to compare to focal node")
    parser.add_argument("--foc", type=str, nargs="*", action="append",
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
    gnn_win = args.gnn_windows
    gnn_time = args.time_windows
    # =========================================================================
    #  Loading and Checks
    # =========================================================================
    ts = tskit.load(tree)  # load tree
    # set refernce for comparison
    if ref_set:
        ref_set = list(map(int, args.ref[0]))
    else:
        ref_set = range(ts.num_populations)  # all populations
    # set target population
    if foc_set:
        foc_set = list(map(int, foc_set[0]))
        target_samples = []
        for i in foc_set:
            target_samples.extend(ts.samples(population=i))
    else:
        target_samples = ts.samples()
    # =========================================================================
    #  Main executions
    # =========================================================================
    ref_samples = [ts.samples(population=i) for i in ref_set]
    ref_groups = [json.loads(ts.population(i).metadata)["Group"] for i in ref_set]
    all_groups = [json.loads(ts.population(i).metadata)["Group"] for i in range(ts.num_populations)]

    if not gnn_win and not gnn_time:
        gnn_fx(outfile, ts, ref_samples, target_samples, ref_groups)
        plot_gnn_wg(f"GNN.{outfile}.csv", ref_groups, FOCAL_IND)
    else:
        assert len(foc_set) == 1, "windows option only works for 1 target set"
        target_group = all_groups[foc_set[0]]
        gnn_m = gnn_windows_fx(outfile, ts, ref_samples, target_samples, ref_groups,
                               target_group, gnn_win, gnn_time)
        #plot_gnn_windows_p1(outfile, ts, gnn_m, ref_groups, target_group)
        plot_gnn_windows_p2(outfile, ts, gnn_m, ref_groups, target_group)


if __name__ == "__main__":
    main()
