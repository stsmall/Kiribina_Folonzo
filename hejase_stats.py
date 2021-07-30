# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 10:21:01 2021

@author: Scott T. Small

This module calculates the tmrca half and cross coalescent 10 statistics
from the paper of Hejase et al. 2020 
(https://www.biorxiv.org/content/10.1101/2020.03.07.977694v2.full).

Example
-------
Examples using relate trees. Note that relate trees need to be converted to 
tree sequencing format using --mode ConvertToTreeSequence.

$ python hejase_stats.py --trees FOO.relate.trees --outfile FOO.1_2 
    --pop_ids 1 2 --node_ids pops.nodes.txt --fx tmrca_half

Notes
-----
The node file input was needed since Relate didnt transfer over information 
about the populations, meaning the populations were not stored in the tree seq.
The node file has a single line of comma delimited integers denoting the leaf
id associated with the desired population or group. The pop_ids need to be in 
the same order as the node file to ensure proper naming.

$ > head FOO.node.txt
1,2,3,4,5,6,7,8
11,12,13,14,15


"""
import argparse
from itertools import combinations
import functools
import numpy as np
from os import path
import pandas as pd
import sys
from tqdm import tqdm
import tskit

import pysnooper


def load_tree(tree):
    """Reads tree sequence from disk.    

    Parameters
    ----------
    tree : str
        file path to tree sequence

    Returns
    -------
    tskit tree sequencing object

    """

    return tskit.load(tree)

@pysnooper.snoop()
def calc_tmrcah(ts, p_nodes):
    """Calculate the tmraca half as defined in Hejase 2020.
    
    "...test on the time to the most recent common ancestor of half the haploid
    samples from a given species (TMRCAH). Requiring only half the samples
    allows us to consider partial sweeps and provides robustness to the
    inherent uncertainty in the inferred local trees."   

    Parameters
    ----------
    ts : Object
        object of type tskit tree seqeunce.
    p_nodes : List
        List of node ids as integers, [[0,1,2],[4,5,6]]

    Returns
    -------
    mid : List
        the center of the interval in base pairs position 
    tmrcah_rel : List
        the tmrca of half the population
    time_rel : List
        Full TMRCA of that population.
    time_rel2 : List
        Age of the youngest subtree that contains at least half of the samples.

    """
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


def tmrca_half(tree_str, pop_nodes, pop_ids, outfile):
    """Calculats the tmrca half fx from Hejase et al 2020.
    
        "...test on the time to the most recent common ancestor of half the haploid
    samples from a given species (TMRCAH). Requiring only half the samples
    allows us to consider partial sweeps and provides robustness to the
    inherent uncertainty in the inferred local trees."   
    
    Parameters
    ----------
    tree_str : str
        file path to tree sequence type.
    pop_nodes : List
        population leaves as integers loaded from file.
    pop_ids : List
        id of population nodes to be written in DataFrame.
    outfile : str
        base name of DataFrame file output.

    Returns
    -------
    None.

    """
    ts = load_tree(tree_str)
    print("tree loaded")
    df_list = []
    for pop, nodes in zip(pop_ids, pop_nodes):
        mid, tmrcah_rel, time_rel, time_rel2 = calc_tmrcah(ts, nodes)
        # set up DataFrame
        df_pop = pd.DataFrame({"population": pd.Series([pop]*len(mid)),
                               "mid": pd.Series(mid), 
                               "tmrcah": pd.Series(tmrcah_rel), 
                               "time_rel": pd.Series(time_rel),
                               "time_rel2": pd.Series(time_rel2)})
        df_list.append(df_pop)
    df_pop_combine = pd.concat(df_list).reset_index(drop=True)
    df_pop_combine.to_csv(f"{outfile}.tmrca_half.csv", na_rep="NAN", index=False)


@pysnooper.snoop()
def calc_cc10(ts, p_nodes_cc, cc_events=10):
    """Calculate the cross coalescent of two populations.
    
    Parameters
    ----------
    ts : Object
        object of type tskit tree seqeunce.
    p_nodes_cc : List
        List of node ids as integers, [[0,1,2],[4,5,6]]
    cc_events : int, optional
        the number of cross coalescent events to track. The default is 10.

    Returns
    -------
    mid : List
        the center of the interval in base pairs position 
    cc10_rel : List
        the cross coalescent of the population
    time_rel : List
        Full TMRCA of that population.

    """
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
        cc10_rel.append(cc10)
        
    return mid, cc10_rel, time_rel


def cross_coal_10(tree_str, pop_nodes, pop_ids, outfile):
    """Calculate the cross coalescent 10 stat from Hejase et al 2020.
    
    "...For a given local tree and pair of species, we considered the 10 most 
    recent cross coalescent events between the two species and normalized these 
    ages, as in test 2, by the age of the youngest subtree that contains at 
    least half of the total number of haploid samples."

    Parameters
    ----------
    tree_str : str
        file path to tree sequence type.
    pop_nodes : List
        population leaves as integers loaded from file.
    pop_ids : List
        id of population nodes to be written in DataFrame.
    outfile : str
        base name of DataFrame file output.

    Returns
    -------
    None.

    """
    ts = load_tree(tree_str)
    print("tree loaded")
    df_list = []
    pop_node_pairs = combinations(pop_nodes, 2)
    pop_ids_pairs = combinations(pop_ids, 2)
    for pop, nodes in zip(pop_ids_pairs, pop_node_pairs):
        mid, cc10_rel, time_rel = calc_cc10(ts, nodes)
        breakpoint()
        # prep df        
        cc10_dt = {f"cc_{i+1}":cc for i, cc in enumerate(zip(*cc10_rel))}
        cc10_cols = list(cc10_dt.keys())
        cc10_dt["time_rel"] = pd.Series(time_rel)
        pop_pair = ["_".join(pop)]
        cc10_dt["populations"] = pd.Series(pop_pair*len(mid))
        cc10_dt["mid"] = pd.Series(mid)
        # save df
        df_pop = pd.DataFrame(data=cc10_dt, columns=["population", "mid", "time_rel"]+cc10_cols) 
        df_list.append(df_pop)
    
    df_pop_combine = pd.concat(df_list).reset_index(drop=True)
    df_pop_combine.to_csv(f"{outfile}.cross_coal10.csv", na_rep="NAN", index=False)


def parse_args(args_in):
    """Parse args."""
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--trees", type=str, required=True,
                        help="file containing ARGs in tskit format")
    parser.add_argument("--outfile", type=str, default=None,
                        help="bas name for output file")
    parser.add_argument("--pop_ids", type=str, nargs="*", action="append",
                        help="pop ids for naming columns in output dataframe")
    parser.add_argument("--node_ids", type=str,
                        help="load pop nodes from this file, one per line and" 
                        "comma delimited")
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
    # load pop nodes from file
    pop_nodes = []
    with open(node_file) as f:
        for line in f:
            x = line.split(",")
            assert len(x) > 1, "recheck nodes file, delimiter should be ,"
            pop_nodes.append(list(map(int, x)))
    assert len(pop_nodes) == len(pop_ids), "some pop nodes dont have names"
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
