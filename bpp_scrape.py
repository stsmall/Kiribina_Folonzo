#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 17:52:02 2019
@author: Scott T. Small

BppScrape.py -p nCDS.bpp. -s .out -c 100000 --scafs 2L 2R 3L 3R X
returns: Coordinate File, Twisst-style weights file

"""
import argparse
import glob
import re
import sys
from collections import defaultdict

import numpy as np


def scrape_bpp(prefix: str,
               suffix: str,
               chain: int,
               scafs: str):
    """Find and collect topo inforation from BPP file run with A01

    Parameters
    ----------
    prefix: str
        file name startswith()
    suffix: str
        file name endswith()
    chain: int
        mcmc chain length else will print proportions
    scafs: str
        name of scaffold to print to coords file

    Returns
    -------
    weights_dict: Dict[default(dict)]
        coords, topo, weights
    topo_list: List[str]
        list of topos found in the file

    """
    topo_list = []
    weights_ddict = defaultdict(dict)
    for sc in scafs:
        file_list = glob.glob(f"{prefix}{sc}*{suffix}")
        for bpp_out in file_list:
            coord = re.search(r"([0-9]+-[0-9]+)", bpp_out).group()
            start, stop = coord.split("-")
            keycoord = f"{sc}:{start}-{stop}"
            with open(bpp_out, "r") as bpp:
                for line in bpp:
                    if line.startswith("(A)"):
                        line = next(bpp)
                        while line.strip() != "":
                            x = line.split()
                            topo = "".join(x[3:])
                            if chain > 0:
                                weights_ddict[keycoord][topo] = float(x[0])
                            else:
                                weights_ddict[keycoord][topo] = float(x[1])
                            topo_list.append(topo)
                            line = next(bpp)
                        break
    return weights_ddict, topo_list


def write_weights(weights_ddict, topo_list):
    """Write weights in stype similar to twisst output

    Parameters
    ----------
    weights_ddict: default(dict)
        coords, topo, weights
    topo_list: List[str]
        list of topos found in the file

    Returns
    -------
    None

    """
    topo_set = list(set(topo_list))
    topo_header = ""
    with open("topos.out", "w") as topo_file:
        for i, topo in enumerate(topo_set):
            topo_file.write(f"#topo{i+1}\t{topo}\n")
            topo_header += f"topo{i+1}\t"
        topo_header.rstrip("\t")

    topo_count = len(topo_set)
    coord_count = len(weights_ddict)
    topo_freq = np.zeros([coord_count, topo_count])
    with open("weights.out", "w") as weights_file:
        weights_file.write(f"scaf\tstart\tstop\t{topo_header}\n")
        coord_list = list(weights_ddict.keys())
        coord_sort = sorted(coord_list,
                            key=lambda x: int(re.search(r"([0-9]+)\-", x).group(1)))
        for i, coord in enumerate(coord_sort):
            scaf, start, stop = re.findall(r"\w+", coord)
            topo_weights = [0] * topo_count
            for topo in weights_ddict[coord]:
                ix = topo_set.index(topo)
                topo_weights[ix] = weights_ddict[coord][topo]
            topo_freq[i] = topo_weights
            pweights = '\t'.join(map(str, topo_weights))
            weights_file.write(f"{scaf}\t{start}\t{stop}\t{pweights}\n")
    np.savetxt("topo_freq.out", np.sum(topo_freq, axis=0)/coord_count)


def parse_args(args_in):
    parser = argparse.ArgumentParser(prog="sys.argv[0].py",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-p', "--prefix", default="",
                        help="prefix of input file")
    parser.add_argument('-s', "--suffix", default=".out",
                        help="suffix of input file")
    parser.add_argument("--scafs", required=True, nargs='+', action="append",
                        help="scaffold or chromosome")
    parser.add_argument('-c', "--chainLen", default=0, type=int,
                        help="chain length, default or value of 0 will return proportions")
    return(parser.parse_args(args_in))


if __name__ == "__main__":
    # =========================================================================
    #  Gather args
    # =========================================================================
    args = parse_args(sys.argv[1:])
    PREFIX = args.prefix
    SUFFIX = args.suffix
    CHAIN_LEN = args.chainLen
    SCAF = args.scafs[0]
    # =========================================================================
    #  Main executions
    # =========================================================================
    WEIGHTS_DDICT, TOPO_LIST = scrape_bpp(PREFIX, SUFFIX, CHAIN_LEN, SCAF)
    write_weights(WEIGHTS_DDICT, TOPO_LIST)
