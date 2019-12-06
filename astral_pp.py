#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 17:28:18 2019
@author: Scott T. Small

astral_pp.py is a program to construct local species trees from collections
of gene trees distributed along the genome. This idea is similar to that
of Thawornwattana et al 2018, where they used BPP for this same thing. ASTRAL
is much faster and therefore it is easier to explore the clustering number for
fine scale analysis.

Example
-------

    $ python astral_pp.py --trees FOO.newick --coords BAR.coords.txt
    [--clust 50] [--astral_exe ./path2astral] [--scafs CHROM]

Notes
-----
    Since we are using ASTRAL this can utilize trees with >1 leaf per species,
    it also returns a species trees with all the same leaves. I then use the
    program twisst (Simon Martin) to represent the topologies along the
    chromosome. Alternately, if you set [groups] and [outgroup] the output
    will be a species tree with 1 leaf per species.

    ls -1v *.bestTree | cut -d"." -f4 > ../2L.cds.coords.txt

    for i in $( ls -1v *.bestTree);do
    cat $i >> 2L.cds.trees
    done


"""
import argparse
from os import path
import subprocess
import sys


def run_astral(tree_file: str,
               clust: int,
               astralexe: str,
               groups: str):
    """Take groups of trees by line and runs them through ASTRAL

    Parameters
    ----------
    tree_file: str, file
        file of newick trees, 1 per line
    clust: int
        number of loci to push into 1 file
    astralexe: str
        path to ASTRAL executable
    groups: str, file
        name of file with groups links

    Returns
    -------
    None

    """
    tree_list = []
    with open(tree_file) as trees:
        for line in trees:
            tree_list.append(line.strip())
    file1 = open("astral.tre", 'a')
    start = 0
    step = clust
    end = start + step
    tree_slice = tree_list[start:end]
    while tree_slice:
        with open("astral_tmp.tre", 'w') as file2:
            for tree in tree_slice:
                file2.write(f"{tree}\n")
        # run astral
        if groups:
            command = f"java -jar {astralexe} -i astral_tmp.tre -o astral.out -a {groups}"
        else:
            command = f"java -jar {astralexe} -i astral_tmp.tre -o astral.out"
        proc = subprocess.Popen(command, shell=True)
        proc.wait()
        # read astral file
        with open("astral.out", 'r') as astral_out:
            for line in astral_out:
                file1.write(line)
        start = end
        end += step
        tree_slice = tree_list[start:end]
    file1.close()
    return(None)


def make_windows(coord_file: str,
                 clust: int,
                 scaf: str):
    """Takes a list of coordinates corresponding to a file of trees, 1 per line
    and clusters them into groups for plotting

    Parameters
    ---------
    coordlist: str, file
        file with coordinates  corresponding to each line in tree_file
    clust: int
        number of loci to push into 1 file
    scaf: str
        name of scaffold or chromosome

    Returns
    -------
    None

    """
    start_list = []
    end_list = []
    step = clust
    with open(coord_file, 'r') as coords:
        for line in coords:
            s, e = line.strip().split("-")
            start_list.append(s)
            end_list.append(e)
    with open(f"{scaf}.windows.out", 'w') as file1:
        # clust coords
        s_ix = 0
        e_ix = s_ix + step
        while e_ix < len(end_list):
            start = start_list[s_ix]
            end = end_list[e_ix - 1]
            file1.write(f"{scaf}\t{start}\t{end}\n")
            s_ix = e_ix
            e_ix = s_ix + step
        else:
            start = start_list[s_ix]
            end = end_list[-1]
            file1.write(f"{scaf}\t{start}\t{end}\n")
    return(None)


def parse_args(args):
    parser = argparse.ArgumentParser(prog="astral_pp.py",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--trees", type=str, required=True,
                        help="trees file")
    parser.add_argument("--coords", type=str,
                        help="list of tree files, same order")
    parser.add_argument("--clust", type=int, default=100,
                        help="how many loci to cluster for astral")
    parser.add_argument("--astral_exe", type=str, required=True,
                        help="path to astral exec")
    parser.add_argument("--scafs", type=str, required=True,
                        help="scaf or chrom name for output")
    parser.add_argument("--groups", type=str,
                        help="file clustering tips into groups, use -a ASTRAL")
    parser.add_argument("--outgroup", type=str,
                        help="rooting using newick_utils")
    return(parser.parse_args(args))


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    TREES = args.trees
    CLUST = args.clust
    EXE = args.astral_exe
    GROUPS = args.groups
    OUT = args.outgroup
    SCAF = args.scafs
    COORD = args.coords
    run_astral(TREES, CLUST, EXE, GROUPS)
    make_windows(COORD, CLUST, SCAF)
    # reroot using newick utils
    NEWICK_UTILS_PATH = "/afs/crc.nd.edu/user/s/ssmall2/programs_that_work/newick-utils-1.6/src/nw_reroot"
    if OUT:
        if path.exists(NEWICK_UTILS_PATH):
            command = f"{NEWICK_UTILS_PATH} astral.tre {OUT} > {SCAF}.astral.{CLUST}.rooted.tre"
            proc = subprocess.Popen(command, shell=True)
            proc.wait()
        else:
            f"no newick-utils found in path, not rooting!!!!!"
            command = f"mv astral.tre {SCAF}.astral.{CLUST}.unrooted.tre"
            proc = subprocess.Popen(command, shell=True)
            proc.wait()
    else:
        command = f"mv astral.tre {SCAF}.astral.{CLUST}.unrooted.tre"
        proc = subprocess.Popen(command, shell=True)
        proc.wait()
    proc = subprocess.Popen(f"rm -f astral.out astral_tmp.tre", shell=True)
