#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 14:09:42 2019
@author: Scott T. Small


Example
-------
Examples can be given using either the ``Example`` or ``Examples``
sections. Sections support any reStructuredText formatting, including
literal blocks::

    $ python make


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
import glob
import sys
from typing import Dict, List
import os
import random
import re
import argparse
from tqdm import tqdm


class Sequence():
    """The Sequence object has a string *header* and
    various representations."""

    def __init__(self, header, seq):
        self.header = re.findall(r"^>(\S+)", header)[0]
        self.seq = seq

    def __len__(self):
        return len(self.seq)

    @property
    def phylip(self):
        return self.header + " " + self.seq.replace('.', '-') + "\n"

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
    header = ''
    seq = ''
    with open(fasta_file, 'r') as fa:
        for line in fa:
            line = line.strip('\n')
            if line.startswith('>'):
                if header:
                    yield Sequence(header, seq)
                header = line
                seq = ''
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
        tmp_path = phylip + '.' + str(random.randint(100000000, 1000000000))
        count = 0
        with open(tmp_path, 'w') as tmp_file:
            for seq in seqs:
                tmp_file.write(seq.phylip)
                count += 1
                len_seq = len(seq)
        with open(tmp_path, 'r') as old, open(phylip, 'w') as new:
            new.write(f" {count} {len_seq}\n")
            new.writelines(old)
        os.remove(tmp_path)
        phy_list.append(phylip)
    pbar.close()
    return(phy_list)


def write_to_clust(clust_files: List[str],
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
    start = re.findall(r"[\.\-]?([0-9]+)[\.\-]", first_file)[0]
    end = re.findall(r"[\.\-]?([0-9]+)[\.\-]", last_file)[1]
    with open(f"{chrom}.{cds_type}.{start}-{end}.{count}.txt", 'w') as f:
        for aln in clust_files:
            with open(aln, 'r') as locus:
                if aln.endswith(".fa"):
                    header = locus.next()
                    seq = ''
                    for line in locus:
                        line = line.strip('\n')
                        if line.startswith('>'):
                            if map_dict:
                                header = map_dict[header]
                            f.write(f"{header}\n")
                            f.write(f"{seq}\n")
                            header = line
                            seq = ''
                        else:
                            seq += line
                else:
                    header = locus.next()
                    f.write(f"\n{header}\n")
                    for line in locus:
                        seq = line.split()
                        spname = seq[0]
                        if map_dict:
                            spname = map_dict[spname]
                        dna = seq[1]
                        if bpp:
                            f.write(f"^{spname[:10]:<10}{dna}\n")
                        else:
                            if strict:
                                f.write(f"{spname[:10]:<10}{dna}\n")
                            else:
                                f.write(f"{spname} {dna}\n")
        f.write("\n")
    return(chrom, cds_type)


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
    pbar = tqdm(total=len(file_list))

    map_dict = {}
    with open(map_file, 'r') as new_ids:
        for line in new_ids:
            ids, name = line.split()
            map_dict[name] = ids
    pattern = r".([0-9]+)-"
    aln_list = [os.path.basename(x) for x in file_list]
    lsorted = sorted(aln_list, key=lambda x: int(re.search(pattern, x).group(1)))
    s_ix = 0
    e_ix = cluster
    while e_ix < len(lsorted):
        pbar.update(cluster)
        clust_files = lsorted[s_ix:e_ix]
        chrom, cds_type = write_to_clust(clust_files, map_dict, bpp)
        s_ix = e_ix
        e_ix = s_ix + cluster
    # last window
    if s_ix < len(lsorted):
        e_ix = len(lsorted)
        clust_files = lsorted[s_ix:e_ix]
        chrom, cds_type = write_to_clust(clust_files, map_dict, bpp)
    pbar.close()
    return(chrom, cds_type)


def make_bpp(chrom: str,
             cds_type: str,
             ctl_file: str):
    """
    """
    clust_files = glob.glob(f"{chrom}.{cds_type}*.txt")
    for loci in clust_files:
        coords = loci.split(".")[2]
        loci_count = loci.split(".")[-2]
        out_file = f"{chrom}.{cds_type}.{coords}"
        with open(f"A01.bpp.{out_file}.ctl", 'w') as f:
            with open(ctl_file, 'r') as ctl:
                for line in ctl:
                    if line.startswith("seqfile"):
                        f.write(f"seqfile = {loci}\n")
                    elif line.startswith("outfile"):
                        f.write(f"outfile = bpp.{out_file}.out\n")
                    elif line.startswith("mcmcfile"):
                        f.write("mcmcfile = bpp.{out_file}.msmc.txt\n")
                    elif line.startswith("nloci"):
                        f.write(f"nloci = {loci_count}\n")
                    else:
                        f.write(line)
    return(None)


def parse_args(args_in):
    parser = argparse.ArgumentParser(prog=f"{sys.argv[0]}", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--aln_path", type=str, required=True,
                        help="path to files")
    parser.add_argument("--aln_format", type=str, required=True,
                        choices=("fa", "phy"),
                        help="format of alignment file, fasta or phylip")
    parser.add_argument("--cluster", type=int, default=100,
                        help="how many loci to cluster in sinlge file")
    parser.add_argument("--bpp", action="store_true",
                         help="make clusters for BPP program, requires example control file")
    parser.add_argument("--control_file", type=str,
                        help="control file for BPP")
    parser.add_argument("--map_file", type=str,
                        help="map file for changing IDs")
    return(parser.parse_args(args_in))


def main(args_in):
    args_out = parse_args(args_in)
    return(args_out)


if __name__ == "__main__":
    args = main(sys.argv[1:])
# =============================================================================
#  Gather args
# =============================================================================
    ALN_PATH = args.aln_path
    ALN_FORMAT = args.aln_format
    CLUSTER = args.cluster
    MAKE_BPP = args.bpp
    CNTRL_BPP = args.control_file
    MAP_IDS = args.map_file
# =============================================================================
#  Main executions
# =============================================================================
    FILE_LIST = glob.glob(f"{ALN_PATH}/*.{ALN_FORMAT}")
    if len(FILE_LIST) == 0:
        raise ValueError("file list is empty, check path and extensions")
    if MAKE_BPP:
        if ALN_FORMAT == "fa":
            print(f"converting to phylip format for BPP")
            PHYLIST = convert_fa_to_phy(FILE_LIST)
        if (CNTRL_BPP is not None):
            CHROM, CDS_TYPE = cluster_alnments(PHYLIST, CLUSTER, MAP_IDS, MAKE_BPP)
            make_bpp(CHROM, CDS_TYPE, CNTRL_BPP)
        else:
            raise ValueError("BPP argument requires control and map files")
    else:
        cluster_alnments(FILE_LIST, CLUSTER, MAP_IDS, MAKE_BPP)
