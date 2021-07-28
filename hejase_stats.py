# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 10:21:01 2021

@author: Scott T. Small

This module demonstrates documentation as specified by the `NumPy
Documentation HOWTO`_. Docstrings may extend over multiple lines. Sections
are created with a section header followed by an underline of equal length.

Example
-------
Examples can be given using either the ``Example`` or ``Examples``
sections. Sections support any reStructuredText formatting, including
literal blocks::

    $ python example_numpy.py


Section breaks are created with two blank lines. Section breaks are also
implicitly created anytime a new section starts. Section bodies *may* be
indented:

Notes
-----
    This is an example of an indented section. It's like any other section,
    but the body is indented to help it stand out from surrounding text.

If a section is indented, then a section break is created by
resuming unindented text.

Attributes
----------
module_level_variable1 : int
    Module level variables may be documented in either the ``Attributes``
    section of the module docstring, or in an inline docstring immediately
    following the variable.

    Either form is acceptable, but the two should not be mixed. Choose
    one convention to document module level variables and be consistent
    with it.

"""
import argparse
import sys
from os import path
import multiprocessing
import tskit
import pandas as pd
import numpy as np
import math
import functools
import networkx as nx
from itertools import product, combinations


def load_tree(tree):
    ts = tskit.load(tree)
    return ts


def tmrca_half_parallel_v2(tree_ix):
    mid = []
    tmrcah_rel = []
    time_rel = []
    sample_half = trees.num_samples / 2
    for ix in tree_ix:
        t = trees.at_index(ix)
        mid.append(((t.interval[1] - t.interval[0]) / 2) + t.interval[0])
        tmrcah = np.inf
        for n in t.nodes(order='timeasc'):
            if t.num_samples(n) >= p_half:
                count_pop = len(list(set(list(t.leaves(n))) & set(p_nodes)))
                if count_pop >= p_half:
                    tmrcah = t.time(n)
                    break
        for n in t.nodes(order='timeasc'):
            if t.num_samples(n) >= sample_half:
                time_r = t.time(n)
                break
        tmrcah_rel.append(tmrcah)
        time_rel.append(time_r)
        
    return mid, tmrcah_rel, time_rel


def tmrca_half_parallel_v1(tree_ix):
    mid = []
    tmrcah_rel = []
    time_rel = []
    for ix in tree_ix:
        t = trees.at_index(ix)
        mid.append(((t.interval[1] - t.interval[0]) / 2) + t.interval[0])
        tmrcah = np.inf
        for n in t.nodes(order='timeasc'):
            if t.num_samples(n) >= p_half:
                count_pop = len(list(set(list(t.leaves(n))) & set(p_nodes)))
                if count_pop >= p_half:
                    tmrcah = t.time(n)
                    break
        mrca = functools.reduce(t.mrca, p_nodes)
        tmrcah_rel.append(tmrcah)
        time_rel.append(t.time(mrca))
    return mid, tmrcah_rel, time_rel

def tmrca_half_parallel_v1_b(tree_ix):
    mid = []
    tmrcah_rel = []
    time_rel = []
    for ix in tree_ix:
        t = trees.at_index(ix)
        #t = t.subset(p_nodes)  # how do I take a subset of just a tree?
        mid.append(((t.interval[1] - t.interval[0]) / 2) + t.interval[0])
        tmrcah = np.inf
        for n in t.nodes(order='timeasc'):
            if t.num_samples(n) >= p_half:
                tmrcah = t.time(n)
                break  
        tmrcah_rel.append(tmrcah)
        time_rel.append(t.time(t.root))

    return mid, tmrcah_rel, time_rel


def tmrca_half(tree_str, pop_nodes, pop_ids, outfile="Out", nprocs=4, version=1):
    c_per_proc = 100  # chunk per processor
    # tmrcah from hejase and ref 44 therein    
    ts = load_tree(tree_str)
    print("tree loaded")
    global p_half
    global p_nodes
    global trees
    trees = ts
    n_trees = ts.num_trees
    tree_ix = list(range(0, n_trees))
    df_list = []
    for pop, nodes in zip(pop_ids, pop_nodes):
        mid = []
        tmrcah_rel = []
        time_rel = []
        p_half = len(nodes) / 2
        p_nodes = nodes

        if version == 1:
            # chunk and MP
            nk = nprocs * c_per_proc
            chunk_list = [tree_ix[i:i + nk] for i in range(0, n_trees, nk)]
            chunksize = math.ceil(nk/nprocs)
            with multiprocessing.Pool(nprocs) as pool:
                for i, args in enumerate(chunk_list):
                    mid_i, tmrcah_i, time_i = pool.map(tmrca_half_parallel_v1, args, chunksize=chunksize)
                    mid.extend(mid_i)
                    tmrcah_rel.extend(tmrcah_i)
                    time_rel.extend(time_i)
                    print(f"{100*(i/len(chunk_list))} percent complete")

        elif version == 2:
            # chunk and MP
            nk = nprocs * c_per_proc
            chunk_list = [tree_ix[i:i + nk] for i in range(0, n_trees, nk)]
            chunksize = math.ceil(nk/nprocs)
            with multiprocessing.Pool(nprocs) as pool:
                for i, args in enumerate(chunk_list):
                    mid_i, tmrcah_i, time_i = pool.map(tmrca_half_parallel_v2, args, chunksize=chunksize)
                    mid.extend(mid_i)
                    tmrcah_rel.extend(tmrcah_i)
                    time_rel.extend(time_i)
                    print(f"{100*(i/len(chunk_list))} percent complete")

        df_pop = pd.DataFrame({"population": pd.Series(pop*len(mid)),
                               "mid": pd.Series(mid), 
                               "tmrcah": pd.Series(tmrcah_rel), 
                               "time_rel": pd.Series(time_rel)})
        df_list.append(df_pop)
    df_pop_combine = pd.concat(df_list).reset_index(drop=True)
    df_pop_combine.to_csv(f"{outfile}.tmrca_half.csv", na_rep="NAN", index=False)
    
 
def cross_coal_10_parallel(tree_ix):
    mid = []
    cc10_rel = []
    time_rel = []
    sample_half = trees.num_samples / 2
    for ix in tree_ix:
        t = trees.at_index(ix)
        mid.append(((t.interval[1] - t.interval[0]) / 2) + t.interval[0])
        td = nx.DiGraph(t.as_dict_of_dicts())
        cc = list(nx.all_pairs_lowest_common_ancestor(td, list(product(p_nodes_cc[0], p_nodes_cc[1]))))
        cc10 =  np.mean(np.sort([t.time(i[1]) for i in cc])[:10])
        for n in t.nodes(order='timeasc'):
            if t.num_samples(n) > sample_half:
                cc10_rel.append(cc10)
                time_rel.append(t.time(n))
                break
    return mid, cc10_rel, time_rel


def cross_coal_10(tree_str, pop_nodes, pop_ids, outfile="Out", nprocs=4):
    c_per_proc = 10  # chunk per processor
    ts = load_tree(tree_str)
    print("tree loaded")
    n_trees = ts.num_trees
    tree_ix = list(range(0, n_trees))
    df_list = []
    global p_nodes_cc
    global trees
    trees = ts
    pop_node_pairs = combinations(pop_nodes, 2)
    pop_ids_pairs = combinations(pop_ids, 2)
    for pop, nodes in zip(pop_ids_pairs, pop_node_pairs):
        pop_pair = ["_".join(pop)]
        p_nodes_cc = nodes
        mid = []
        cc10_rel = []
        time_rel = []
        # chunk and MP
        nk = nprocs * c_per_proc
        chunk_list = [tree_ix[i:i + nk] for i in range(0, n_trees, nk)]
        chunksize = math.ceil(nk/nprocs)

        with multiprocessing.Pool(nprocs) as pool:
            for i, args in enumerate(chunk_list):
                mid_i, tmrcah_i, time_i = pool.map(cross_coal_10_parallel, args, chunksize=chunksize)
                mid.extend(mid_i)
                cc10_rel.extend(tmrcah_i)
                time_rel.extend(time_i)
                print(f"{100*(i/len(chunk_list))} percent complete")
        
        df_pop = pd.DataFrame({"population": pd.Series(pop_pair*len(mid)),
                       "mid": pd.Series(mid), 
                       "cross_coal10": pd.Series(cc10_rel), 
                       "time_rel": pd.Series(time_rel)})
        df_list.append(df_pop)
    
    df_pop_combine = pd.concat(df_list).reset_index(drop=True)
    df_pop_combine.to_csv(f"{outfile}.cross_coal10.csv", na_rep="NAN", index=False)


def parse_args(args_in):
    """Parse args."""
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--trees", type=str, required=True,
                        help="file containing ARGs")
    parser.add_argument("--outfile", type=str, default=None,
                        help="name for output file")
    parser.add_argument("--pop_ids", type=str, nargs="*", action="append",
                        help="pop ids for naming columns")
    parser.add_argument("--node_ids", type=str,
                        help="load pop nodes from this file")
    parser.add_argument("--np", type=int, default=4,
                        help="number of proicessors")
    parser.add_argument("--fx", type=str, default=None,
                        choices=("tmrca_half", "cross_coal_10"),
                        help="which fx to run ... since they both can take a long time")
    return(parser.parse_args(args_in))


def main():
    """Run main function."""
    args = parse_args(sys.argv[1:])
    # =========================================================================
    #  Gather args
    # =========================================================================
    args_file = args.trees
    outfile = args.outfile
    if outfile is None:
        outfile = path.split(args_file)[-1]
    pop_ids = args.pop_ids[0]
    node_file = args.node_ids
    nprocs = args.np
    fx = args.fx
    # load pop nodes
    pop_nodes = []
    with open(node_file) as f:
        for line in f:
            x = line.split(",")
            assert len(x) > 1, "recheck delimiter should be ,"
            pop_nodes.append(list(map(int, x)))
    assert len(pop_nodes) == len(pop_ids)
    # =========================================================================
    #  Main executions
    # =========================================================================
    if fx is None:
        tmrca_half(args_file, pop_nodes, pop_ids, outfile, nprocs)
        cross_coal_10(args_file, pop_nodes, pop_ids, outfile, nprocs)
    elif fx == "tmrca_half":
        tmrca_half(args_file, pop_nodes, pop_ids, outfile, nprocs)
    elif fx == "cross_coal_10":
        cross_coal_10(args_file, pop_nodes, pop_ids, outfile, nprocs)
    else:
        print("fx not recognized")


if __name__ == "__main__":
    main()
