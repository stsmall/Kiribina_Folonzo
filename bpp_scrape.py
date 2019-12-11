#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 17:52:02 2019
@author: Scott T. Small

BppScrape.py -p nCDS.bpp. -s .out -c 100000 --scafs 2L 2R 3L 3R X
returns: Coordinate File, Twisst-style weights file

"""
import sys
import argparse
import glob
from collections import defaultdict
import re


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
            coord = bpp_out.replace(prefix, "").replace(suffix, "")
            scaf, start, stop = re.findall(r'\w+', coord)
            keycoord = f"{sc}:{start}-{stop}"
            with open(bpp_out, "r") as bpp:
                for line in bpp:
                    if line.startswith("(A)"):
                        for line in bpp:
                            if line.strip() == "":
                                break
                            else:
                                x = line.split()
                                topo = "".join(x[3:])
                                if chain > 0:
                                    weightsDict[keycoord][topo] = float(x[1]) * chain
                                else:
                                    weightsDict[keycoord][topo] = float(x[1])
                                topoList.append(topo)
    return (weights_dict, topo_list)


def write_weights(weights_dict, topo_list):
    """Write weights in stype similar to twisst output

    Parameters
    ----------
    weights_dict: default(dict)
        coords, topo, weights
    topo_list: List[str]
        list of topos found in the file

    Returns
    -------
    None

    """
    topo_file = open("topos.out", "w")
    topo_set = list(set(topo_list))
    topo_header = ""
    topo_count = len(topo_set)
    for i, topo in enumerate(topo_set):
        topo_file.write(f"#topo{i+1}\t{topo}\n")
        topo_header += f"topo{i+1}\t"
    topo_header.rstrip("\t")
    topo_file.close()

    weights_file = open("weights.out", "w")
    weights_file.write(f"scaf\tstart\tstop\t{topo_header}\n")
    coord_list = list(weights_dict.keys())
    coord_sort = sorted(coord_list, key=lambda x: int(re.search(r"([0-9]+)\-", x).group(1)))
    for coord in coord_sort:
        scaf, start, stop = re.findall(r"\w+", coord)
        topo_weights = [0] * topo_count
        for topo in weights_dict[coord]:
            ix = topo_set.index(topo)
            topo_weights[ix] = weights_dict[coord][topo]
        weights_file.write(f"{scaf}\t{start}\t{stop}\t{'\t'.join(map(str, topo_weights))}\n")
    weights_file.close()


def parse_args(args_in):
    parser = argparse.ArgumentParser(prog="sys.argv[0].py",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-p', "--prefix", required=True,
                        help="prefix of input file")
    parser.add_argument('-s', "--suffix", required=True,
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
    WEIGHTS_DICT, TOPO_LIST = scrape_bpp(PREFIX, SUFFIX, CHAIN_LEN, SCAF)
    write_weights(WEIGHTS_DICT, TOPO_LIST)