#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 17:01:00 2020
@author: Scott T. Small

This module is a direct implementation of:
    https://github.com/tskit-dev/tskit/issues/665

Example
-------

    $ python windowed_gnn.py --tree chrX.nosingle.tsinfer.trees --pop_ids K F
    $ python windowed_gnn.py --tree chrX.nosingle.tsinfer.trees --tar 0 --pop_ids K F --gnn_windows

Notes
-----


"""
import sys
from os import path
import argparse
import tskit
import json
import pandas as pd
import numpy as np


def gnn_fx(outfile, ts, ref_samples, target_samples, pop_ids):
    """Run mean GNN fx.

    Parameters
    ----------
    outfile : TYPE
        DESCRIPTION.
    ts : TYPE
        DESCRIPTION.
    ref_samples : TYPE
        DESCRIPTION.
    target_samples : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    # calc gnn
    gnn = ts.genealogical_nearest_neighbours(target_samples, ref_samples)
    
    # write out df
    sample_nodes = [ts.node(n) for n in ts.samples()]
    sample_ids = [n.id for n in sample_nodes]
    sample_names = [json.loads(ts.individual(n.individual).metadata)['Isolate'] for n in sample_nodes]
    sample_pops = [json.loads(ts.population(n.population).metadata)['Groups'] for n in sample_nodes]
    gnn_table = pd.DataFrame(data=gnn,
                             index=[pd.Index(sample_ids, name="Node")],
                             columns=pop_ids)
    gnn_table["Isolate"] = sample_names 
    gnn_table["Group"] = sample_pops
    gnn_table.to_csv(f"{outfile}.gnn.csv")


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


def windowed_gnn(ts, focal, reference_sets, windows=None, time_windows=None, span_normalise=True, time_normalise=True):
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


def gnn_windows_fx(outfile, ts, ref_samples, target_samples, pop_ids):
    """Calculate gnn in windows.

    Parameters
    ----------
    ts : Iterator
        tskit object, iterator of trees
    ref_samples : List
        list of reference sample nodes; [[0,2,3],[9,10,11]]
    target_sample : List
        list of target samples; [0,1,2,3]
    pop_ids : List
        list of population ids; ["K", "F"]

    Returns
    -------
    None.

    """
    windows = list(ts.breakpoints())  # all trees
    gnn_win = windowed_gnn(ts, target_samples, ref_samples, windows=windows)

    # save as df
    sample_nodes = [ts.node(n) for n in target_samples]
    sample_names = [json.loads(ts.individual(n.individual).metadata)['Isolate'] for n in sample_nodes]
    sample_names = list(dict.fromkeys(sample_names))
    col_names = [f"{n}_{i}" for n in sample_names for i in [0, 1]]
    iterables = [col_names, pop_ids]
    index = pd.MultiIndex.from_product(iterables, names=["Isolate", "Group"])
    gnn_table = pd.DataFrame(data=np.reshape(gnn_win,[len(gnn_win), np.product(gnn_win.shape[1:])]), columns=index)
    gnn_table.insert(loc=0, column="left_bp", value = list(ts.breakpoints())[:-1])
    gnn_table.insert(loc=1, column="right_bp", value = list(ts.breakpoints())[1:])
    
    gnn_table.to_csv(f"{outfile}.gnn_windows.csv")
    


def parse_args(args_in):
    """Parse args."""
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--tree", type=str, help="tskit tree object")
    parser.add_argument("--tar", type=str, default=None,
                        help="target nodes")
    parser.add_argument("--ref", type=str, default=None, 
                        help="reference nodes")
    parser.add_argument("--pop_ids", type=str, nargs="*", action="append",
                        help="pop ids for naming columns")
    parser.add_argument("--gnn_windows", action="store_true",
                        help="run gnn in windows mode")
    parser.add_argument("--outfile", type=str, default=None,
                        help="name for output file")
    return parser.parse_args(args_in)


def main():
    """Run main function."""
    args = parse_args(sys.argv[1:])
    # =========================================================================
    #  Gather args
    # =========================================================================
    tree = args.tree
    outfile = args.outfile
    if outfile is None:
        outfile = path.split(tree)[-1]
    ref_set = args.ref
    tar_set = args.tar
    gnn_win = args.gnn_windows
    pop_ids = args.pop_ids[0]
    # =========================================================================
    #  Loading and Checks
    # =========================================================================
    ts = tskit.load(tree)  # load tree
    print("tree loaded")
    # set refernce for comparison
    if ref_set:  # custom ref sets for comparison
        ref_nodes = []
        with open(ref_set) as f:
            for line in f:
                x = line.split(",")
                assert len(x) > 1, "recheck delimiter should be ,"
                ref_nodes.append(list(map(int, x)))
    else:  # all populations
        ref_nodes = [ts.samples(population=i) for i in range(ts.num_populations)]
    # set target population
    if tar_set is None:
        tar_nodes = ts.samples()
    elif tar_set.isnumeric():
        tar_nodes = ts.samples(population=int(tar_set))
    else:
        tar_nodes = []
        with open(tar_set) as f:
            for line in f:
                x = line.split(",")
                assert len(x) > 1, "recheck delimiter should be ,"
                tar_nodes.extend(list(map(int, x)))

    # =========================================================================
    #  Main executions
    # =========================================================================
    if gnn_win:
        gnn_windows_fx(outfile, ts, ref_nodes, tar_nodes, pop_ids)
    else:
        gnn_fx(outfile, ts, ref_nodes, tar_nodes, pop_ids)


if __name__ == "__main__":
    main()
