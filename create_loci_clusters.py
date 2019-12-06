#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 14:09:42 2019
@author: Scott T. Small

This code is designed to take a set of fasta or phylip files each with 1 locus
and cluster them into 1 file containing many loci. This specific iteration of
the code was designed to generate input files for BPP.

Example
-------

  $ python create_loci_clusters.py --aln_path tests/fasta_files --aln_format fa


Notes
-----

The file name is expected to carry information.

$chrom.$cds_type.$coords.aln.33.fa = 2L.cds.102218-102369.aln.33.fa


"""
import argparse
import glob
import os
import random
import re
import sys
from typing import Dict, List

from tqdm import tqdm


class Sequence:
    """The Sequence object has a string *header* and
    various representations."""

    def __init__(self, header, seq):
        self.header = re.findall(r"^>(\S+)", header)[0]
        self.seq = seq

    def __len__(self):
        return len(self.seq)

    @property
    def phylip(self):
        return self.header + " " + self.seq.replace(".", "-") + "\n"

    @property
    def fasta(self):
        return ">" + self.header + "\n" + self.seq + "\n"


def fasta_parse(fasta_file):
    """Reads the file at *path* and yields Sequence objects in a lazy fashion

    Parameters
    ----------
    fasta_file: str
        single fasta file with 1 locus

    Returns
    -------
    Sequence: obj
        object of class Sequence containing header and sequence data

    """
    header = ""
    seq = ""
    with open(fasta_file, "r") as fa:
        for line in fa:
            line = line.strip("\n")
            if line.startswith(">"):
                if header:
                    yield Sequence(header, seq)
                header = line
                seq = ""
                continue
            seq += line
    yield Sequence(header, seq)


def convert_fa_to_phy(file_list):
    """Convert fasta to a phylip format file

    Parameters
    ----------
    file_list: List[str]
        list of files and paths containing FOO.fa

    Returns
    -------
    phy_list: List[str]
        list of files now if phylip format

    """
    phy_list = []
    file_path = os.path.dirname(file_list[0])
    pbar = tqdm(total=len(file_list))
    for fasta in file_list:
        if not os.path.exists(fasta):
            raise RuntimeError(f"No file at: {fasta}")
        pbar.update(1)
        base = os.path.basename(fasta).replace("fa", "phy")
        phylip = os.path.join(file_path, base)
        seqs = fasta_parse(fasta)
        tmp_path = phylip + "." + str(random.randint(100000000, 1000000000))
        count = 0
        with open(tmp_path, "w") as tmp_file:
            for seq in seqs:
                tmp_file.write(seq.phylip)
                count += 1
                len_seq = len(seq)
        with open(tmp_path, "r") as old, open(phylip, "w") as new:
            new.write(f" {count} {len_seq}\n")
            new.writelines(old)
        os.remove(tmp_path)
        phy_list.append(phylip)
    pbar.close()
    return phy_list


def write_to_clust(aln_path: str,
                   clust_files: List[str],
                   map_dict: Dict[str, str],
                   bpp: bool,
                   strict: bool = False):
    """Writes alignments to a single files

    Parameters
    ----------
    clust_files: List[str]
        list of files that are part of the clustering in the new single file
    map_dict
    strict: bool
        write strict phylip format, 10 characters justified

    Returns
    -------
    None

    """

    count = len(clust_files)
    first_file = clust_files[0]
    last_file = clust_files[-1]
    chrom = first_file.split(".")[0]
    cds_type = first_file.split(".")[1]
    assert (cds_type.isalpha() is True), "check file name, expects : $chrom.$cds_type.$coords.aln.33.fa"
    start = re.findall(r"[\.\-]?([0-9]+)[\.\-]", first_file)[0]
    end = re.findall(r"[\.\-]?([0-9]+)[\.\-]", last_file)[1]
    with open(f"{chrom}.{cds_type}.{start}-{end}.{count}.txt", "w") as f:
        for aln in clust_files:
            with open(f"{os.path.join(aln_path, aln)}", "r") as locus:
                if aln.endswith(".fa"):
                    header = next(locus).strip("\n")
                    seq = ""
                    for line in locus:
                        line = line.strip("\n")
                        if line.startswith(">"):
                            if map_dict:
                                header = map_dict[header]
                            f.write(f"{header}\n")
                            f.write(f"{seq}\n")
                            header = line
                            seq = ""
                        else:
                            seq += line
                    if aln is not clust_files[-1]:
                        f.write("\n")
                else:
                    header = next(locus).strip("\n")
                    f.write(f"{header}\n")
                    for line in locus:
                        seq_phy = line.split()
                        spname = seq_phy[0]
                        if map_dict:
                            spname = map_dict[spname]
                        dna = seq_phy[1]
                        if bpp:
                            f.write(f"^{spname[:10]:<10}{dna}\n")
                        else:
                            if strict:
                                f.write(f"{spname[:10]:<10}{dna}\n")
                            else:
                                f.write(f"{spname} {dna}\n")
                    if aln is not clust_files[-1]:
                        f.write("\n")
        f.write("\n")
    return (chrom, cds_type)


def cluster_alnments(file_list: List[str],
                     cluster: int,
                     map_file: str,
                     bpp: bool):
    """Cluster alignment files

    Parameters
    ----------
    aln_list: List[str]
        list of file paths
    cluster: int
        how many to cluster in a single file

    Returns
    -------
    None

    """
    if map_file:
        map_dict = {}
        with open(map_file, "r") as new_ids:
            for line in new_ids:
                ids, name = line.split()
                map_dict[name] = ids
    else:
        map_dict = None
    pattern = r".([0-9]+)-"
    aln_list = [os.path.basename(x) for x in file_list]
    aln_path = os.path.dirname(file_list[0])
    lsorted = sorted(aln_list, key=lambda x: int(re.search(pattern, x).group(1)))
    s_ix = 0
    e_ix = cluster
    while (e_ix < len(lsorted)):
        clust_files = lsorted[s_ix:e_ix]
        chrom, cds_type = write_to_clust(aln_path, clust_files, map_dict, bpp)
        s_ix = e_ix
        e_ix = s_ix + cluster
    # last window
    if (s_ix < len(lsorted)):
        e_ix = len(lsorted)
        clust_files = lsorted[s_ix:e_ix]
        chrom, cds_type = write_to_clust(aln_path, clust_files, map_dict, bpp)
    return (chrom, cds_type)


def make_bpp(chrom: str,
             cds_type: str,
             ctl_file: str):
    """Creates input file for BPP

    Parameters
    ----------
    chrom: str
        name of chromosome or scaffold
    cds_type: str
        name of cds or ncds
    ctl_file: str
        input formatted control file

    Returns
    -------
    writes out a file

    """
    clust_files = glob.glob(f"{chrom}.{cds_type}*.txt")
    for loci in clust_files:
        coords = loci.split(".")[2]
        assert (coords.split("-")[0].isdigit() is True), "check file name, expects : $chrom.$cds_type.$coords.aln.33.fa"

        loci_count = loci.split(".")[-2]
        out_file = f"{chrom}.{cds_type}.{coords}"
        with open(f"A01.bpp.{out_file}.ctl", "w") as f:
            with open(ctl_file, "r") as ctl:
                for line in ctl:
                    if line.startswith("seqfile"):
                        f.write(f"seqfile = {loci}\n")
                    elif line.startswith("outfile"):
                        f.write(f"outfile = bpp.{out_file}.out\n")
                    elif line.startswith("mcmcfile"):
                        f.write(f"mcmcfile = bpp.{out_file}.mcmc.txt\n")
                    elif line.startswith("nloci"):
                        f.write(f"nloci = {loci_count}\n")
                    else:
                        f.write(line)


def parse_args(args_in):
    """Arg parser"""
    parser = argparse.ArgumentParser(prog=f"{sys.argv[0]}",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--aln_path", type=str, required=True,
                        help="path to files")
    parser.add_argument("--aln_format", type=str, required=True,
                        choices=("fa", "phy"),
                        help="format of alignment file, fasta or phylip")
    parser.add_argument("--cluster", type=int, default=100,
                        help="how many loci to cluster in sinlge file")
    parser.add_argument("--convert", action="store_true",
                        help="convert fa to phy")
    parser.add_argument("--bpp", action="store_true",
                        help="make clusters for BPP program, requires example control file")
    parser.add_argument("--control_file", type=str,
                        help="control file for BPP")
    parser.add_argument("--map_file", type=str,
                        help="map file for changing IDs")
    return (parser.parse_args(args_in))


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    # =============================================================================
    #  Gather args
    # =============================================================================
    ALN_PATH = args.aln_path
    ALN_FORMAT = args.aln_format
    CLUSTER = args.cluster
    MAKE_BPP = args.bpp
    CNTRL_BPP = args.control_file
    MAP_IDS = args.map_file
    CONVERT = args.convert
    # =============================================================================
    #  Main executions
    # =============================================================================
    FILE_LIST = glob.glob(f"{ALN_PATH}/*.{ALN_FORMAT}")
    if len(FILE_LIST) == 0:
        raise ValueError("file list is empty, check path and extensions")
    if MAKE_BPP is True:
        if ALN_FORMAT == "fa":
            print(f"converting to phylip format for BPP")
            PHYLIST = convert_fa_to_phy(FILE_LIST)
        if CNTRL_BPP is not None:
            CHROM, CDS_TYPE = cluster_alnments(PHYLIST, CLUSTER, MAP_IDS, MAKE_BPP)
            make_bpp(CHROM, CDS_TYPE, CNTRL_BPP)
        else:
            raise ValueError("BPP argument requires control and map files")
    else:
        if ALN_FORMAT == "fa":
            if CONVERT is True:
                print(f"converting to phylip format for BPP")
                PHYLIST = convert_fa_to_phy(FILE_LIST)
                cluster_alnments(PHYLIST, CLUSTER, MAP_IDS, MAKE_BPP)
        else:
            cluster_alnments(FILE_LIST, CLUSTER, MAP_IDS, MAKE_BPP)
