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
from tqdm import tqdm
import tskit
import pandas as pd
import numpy as np
import functools
from itertools import combinations


def load_tree(tree):
    ts = tskit.load(tree)
    return ts


def calc_tmrcah(ts, p_nodes):
    mid = []
    tmrcah_rel = []
    time_rel = []
    time_rel2 = []
    p_half = len(p_nodes) / 2
    sample_half = ts.num_samples / 2
    iter1 = ts.trees(tracked_samples=p_nodes)
    for t in tqdm(iter1, total=ts.num_trees):
        mid.append(((t.interval[1] - t.interval[0]) / 2) + t.interval[0])
        tmrcah = None
        time_r = None
        for u in t.nodes(order='timeasc'):
            if t.num_tracked_samples(u) >= p_half and tmrcah is None:
                tmrcah = t.time(u)
            if t.num_samples(u) >= sample_half and time_r is None:
                time_r = t.time(u)
            if tmrcah is not None and time_r is not None:
                break
            
        tmrcah_rel.append(np.around(tmrcah))
        mrca = functools.reduce(t.mrca, p_nodes)
        time_rel.append(np.around(t.time(mrca)))
        time_rel2.append(np.around(time_r))
        
    return mid, tmrcah_rel, time_rel, time_rel2


def tmrca_half(tree_str, pop_nodes, pop_ids, outfile="Out"):
    # tmrcah from hejase and ref 44 therein    
    ts = load_tree(tree_str)
    print("tree loaded")
    df_list = []
    for pop, nodes in zip(pop_ids, pop_nodes):
        mid, tmrcah_rel, time_rel, time_rel2 = calc_tmrcah(ts, nodes)

        df_pop = pd.DataFrame({"population": pd.Series([pop]*len(mid)),
                               "mid": pd.Series(mid), 
                               "tmrcah": pd.Series(tmrcah_rel), 
                               "time_rel": pd.Series(time_rel),
                               "time_rel2": pd.Series(time_rel2)})
        df_list.append(df_pop)
    df_pop_combine = pd.concat(df_list).reset_index(drop=True)
    df_pop_combine.to_csv(f"{outfile}.tmrca_half.csv", na_rep="NAN", index=False)


def calc_cc10(ts, p_nodes_cc, cc_events=10):
    mid = []
    cc10_rel = []
    time_rel = []
    sample_half = ts.num_samples / 2    

    iter1 = ts.trees(tracked_samples=p_nodes_cc[0])
    iter2 = ts.trees(tracked_samples=p_nodes_cc[1])   
    for tree1, tree2 in tqdm(zip(iter1, iter2), total=ts.num_trees):
        mid.append(((tree1.interval[1] - tree1.interval[0]) / 2) + tree1.interval[0])
        cc10 = [] 
        rel = None
        n1 = 0
        n2 = 0
        i = 0
        for u in tree1.nodes(order='timeasc'):
            p1_n = tree1.num_tracked_samples(u)
            p2_n = tree2.num_tracked_samples(u)

            if p1_n > n1 and p2_n > n2 and i < cc_events:
                i += 1
                cc10.append(np.around(tree1.time(u)))
                n1 = p1_n
                n2 = p2_n

            if tree1.num_samples(u) > sample_half:
                if rel is None:
                    rel = np.around(tree1.time(u))
                if i == cc_events:
                    break

        time_rel.append(rel)
        cc10_rel.append(np.around(np.mean(cc10)))
        
    return mid, cc10_rel, time_rel


def cross_coal_10(tree_str, pop_nodes, pop_ids, outfile="Out"):
    ts = load_tree(tree_str)
    print("tree loaded")
    df_list = []
    pop_node_pairs = combinations(pop_nodes, 2)
    pop_ids_pairs = combinations(pop_ids, 2)
    for pop, nodes in zip(pop_ids_pairs, pop_node_pairs):
        pop_pair = ["_".join(pop)]
        mid, cc10_rel, time_rel = calc_cc10(ts, nodes)
        
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
    parser.add_argument("--fx", type=str, default=None,
                        choices=("tmrca_half", "cross_coal_10"),
                        help="which fx to run ... since they both can take a long time")
    return parser.parse_args(args_in)


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
        tmrca_half(args_file, pop_nodes, pop_ids, outfile)
        cross_coal_10(args_file, pop_nodes, pop_ids, outfile)
    elif fx == "tmrca_half":
        tmrca_half(args_file, pop_nodes, pop_ids, outfile)
    elif fx == "cross_coal_10":
        cross_coal_10(args_file, pop_nodes, pop_ids, outfile)
    else:
        print("fx not recognized")


if __name__ == "__main__":
    main()
