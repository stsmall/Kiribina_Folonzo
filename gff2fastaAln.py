#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 12:21:45 2019
@author: Scott T. Small

This module creates input files for the methods of Thawornattana 2018.
Using a gff file it parses out individual alignments of exons and non-coding
sequences.

Example
-------

    $ python gff2fastaAln.py --gff FOO.gff --aln FOO.aln.fa --chromlen int
                            [--minCDS] [--prct] [--distance int] [--mxlength int]
                            [--mnlength] [--bpp] [--clust int]

"""

from Bio import SeqIO
from dataclasses import dataclass
from typing import Dict, List
import sys
import argparse


@dataclass
class CDS:
    __slots__ = ["start", "end"]
    start: int
    end: int


def get_cds(gff_file: str,
            min_cds: int,
            key: str = "CDS"):
    """Converts gff file to dict with coordinates of each CDS

    Parameters
    ----------
    gffFile: File
        file in gff or gff3 format
    minCDS: int
        min length of coding loci

    Returns
    -------
    cds: dict
        dict with values of dataclass
    chrom: str
        value of chromosome

    """
    cds = {}  # type: Dict
    i = 0  # type: int
    with open(gff_file, 'r') as gff:
        for line in gff:
            if not line.startswith("#"):
                x_list = line.split()
                chrom_gff = x_list[0]  # type: str
                feature = x_list[2]  # type: str
                if key in feature:
                    start = int(x_list[3])  # type: int
                    end = int(x_list[4])  # type: int
                    if start > end:  # reverse strand
                        start = int(x_list[4])
                        end = int(x_list[3])
                    if end - start >= min_cds:
                        cds[f"cds_{i}"] = CDS(start-1, end)
                        i += 1
    return(cds, chrom_gff)


def get_ncds(cdsdict: Dict[str, object],
             mxlen: int,
             mnlen: int,
             distance: int,
             chromlen: int):
    """Returns the non-coding loci around the coding sequence

    Parameters
    ----------
    cdsdict: dict
        dict with values of dataclass: CDS.start, CDS.stop
    mxlen: int
        maximum length of non-coding loci
    mnlen: int
        min length of non-coding loci
    distance: int
        distance between non-coding loci
    chromlen: int
        length of chromosome

    Returns
    -------
    non_cds: dict
        dict with values of dataclass

    """
    non_cds = {}
    loci = 0
    for i in range(-1, len(cdsdict)):
        if i == -1:  # start to first CDS
            start = 0
            end = cdsdict[f"cds_0"].start
        else:
            start = cdsdict[f"cds_{i}"].end
            try:
                end = cdsdict[f"cds_{i+1}"].start
            except KeyError:
                end = chromlen
        nstart = start
        nend = nstart + mxlen
        if nend < end:
            while nend < end:
                non_cds[f"ncds_{loci}"] = CDS(nstart, nend)
                nstart = nend + distance
                nend = nstart + mxlen
                loci += 1
            if (end - nstart) >= mnlen:
                non_cds[f"ncds_{loci}"] = CDS(nstart, end)
                loci += 1
        else:
            if (end - start) >= mnlen:
                non_cds[f"ncds_{loci}"] = CDS(start, end)
                loci += 1

#    # end of chrom
#    nstart = end
#    nend = end + mxlen
#    if nend < chromlen:
#        while nend < chromlen:
#            non_cds[f"ncds_{loci}"] = CDS(nstart, nend)
#            nstart += distance
#            nend = nstart + mxlen
#            loci += 1
#        if (chromlen - nstart) >= mnlen:
#            non_cds[f"ncds_{loci}"] = CDS(nstart, end)
#    else:
#        if (chromlen - end) >= mnlen:
#            non_cds[f"ncds_{loci}"] = CDS(end, chromlen)
    return(non_cds)


def format_fasta(fname: str,
                 gff_dict: Dict[str, object],
                 fasta_file: str,
                 clust: int,
                 chrom: str,
                 prct: float,
                 bpp: bool,
                 just: int = 10):
    """Creates BPP input files or if clust=1 and BPP=False, then a fasta file
    for each CDS and non-CDS.

    Parameters
    ----------
    fname: str
        cds or ncds
    gff_dict: dict
        dict with values of dataclass: fname.start, fname.stop
    fastaFile: file
        alignment in fasta format
    clust: int
        number of loci to push into 1 file
    chrom: str
        chromosome name
    prct: float
        cutoff for percent of missing data in alignment
    bpp: bool
        format for BPP True or False
    just: int
        justified for BPP

    Returns
    -------
    None

    """

    fasta_sequences = list(SeqIO.parse(fasta_file, 'fasta'))
    skip_gaps = 0
    loci = 0
    s = gff_dict[f"{fname}_0"].start
    try:
        e = gff_dict[f"{fname}_{clust-1}"].end
    except KeyError:
        e = gff_dict[f"{fname}_{len(gff_dict)-1}"].end
    out_file = open(f"{fname}.bpp.{chrom}.{s}-{e}.txt", 'w')
    for i in range(len(gff_dict)):
        k = f"{fname}_{str(i)}"
        loci_list = []
        header_list = []
        if loci >= clust:
            out_file.close()
            try:
                s = gff_dict[k].start
                e = gff_dict[f"{fname}_{i + clust-1}"].end
            except KeyError:
                e = gff_dict[f"{fname}_{len(gff_dict)-1}"].end
            out_file = open(f"{fname}.bpp.{chrom}.{s}-{e}.txt", 'w')
            loci = 0
        else:
            for fasta in fasta_sequences:
                header, sequence = fasta.id, str(fasta.seq)
                loci_list.append(sequence[gff_dict[k].start:gff_dict[k].end])
                header_list.append(header)
            samples = len(header_list)
            seqlen = len(loci_list[0])
            # Ns check point
            try:
                if any((seqX.count("N")/seqlen) > prct for seqX in loci_list):
                    # print("skipping, too many Ns")
                    skip_gaps += 1
                else:
                    if bpp is True:
                        out_file.write(f"{samples} {seqlen}\n\n")
                    else:
                        out_file.write(f"\n")
                    for head, seq in zip(header_list, loci_list):
                        if bpp is True:
                            out_file.write(f"^{head}{' '*(just-len(head))}{seq}\n")
                        else:
                            out_file.write(f">{head}\n{seq}\n")
                    loci += 1
            except ZeroDivisionError:
                import ipdb; ipdb.set_trace()
    out_file.close()
    print(f"{skip_gaps} regions skipped due to excess N's")
    return(None)


def write_to_bed(fname: str,
                 gff_dict: Dict[str, object],
                 chrom: int):
    """Write the contents of the dicts to file in bed format

    Parameters
    ----------
    fname: str
        cds or ncds
    gff_dict: dict
        dict with values of dataclass: fname.start, fname.stop
    chrom: str
        chromosome name

    Returns
    -------
    None
    """
    with open(f"{fname}.bed", 'w') as out_bed:
        for i in len(gff_dict):
            start = gff_dict(f"{fname}_{i}").start
            end = gff_dict(f"{fname}_{i}").end
            out_bed.write(f"{chrom}\t{start}\t{end}\n")
    return(None)


def parse_args(args):
    """Argument parser
    """
    parser = argparse.ArgumentParser(prog="gff2fastaAln.py",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--gff", type=str, required=True,
                        help="gff file")
    parser.add_argument("--fasta", type=str, required=True,
                        help="alignment in fasta format")
    parser.add_argument("--minCDS", type=int, default=100,
                        help="min lenght for including coding loci")
    parser.add_argument("--prct", type=float, default=0.5,
                        help="percent sequence length to allow that is Ns")
    parser.add_argument("--distance", type=int, default=2000,
                        help="min distance between non-coding loci")
    parser.add_argument("--mxlength", type=int, default=1000,
                        help="max length for non-coding loci")
    parser.add_argument("--mnlength", type=int, default=100,
                        help="min length for non-coding loci")
    parser.add_argument("--chromlen", type=int, required=True,
                        help="length of chromosome or contig")
    parser.add_argument("--bpp", action="store_true",
                        help="bpp input files")
    parser.add_argument("--clust", type=int, default=100,
                        help="how many loci to cluster into 1 file")
    return(parser.parse_args(args))


def main(args):
    """Calls argument parser
    """
    args = parse_args(args)
    return(args)


if __name__ == "__main__":
    args = main(sys.argv[1:])
    GFF_FILE = args.gff
    FASTA_FILE = args.fasta
    MIN_LEN_CDS = args.minCDS
    PRCT_MISS = args.prct
    DIST_BETW = args.distance
    MAX_LEN = args.mxlength
    MIN_LEN = args.mnlength
    CHROM_LEN = args.chromlen
    BPP = args.bpp
    CLUST = args.clust
    cds_dict, n_chrom = get_cds(GFF_FILE, MIN_LEN_CDS)
    ncds_dict = get_ncds(cds_dict, MAX_LEN, MIN_LEN, DIST_BETW, CHROM_LEN)
    format_fasta("cds", cds_dict, FASTA_FILE, CLUST, n_chrom, PRCT_MISS, BPP)
    format_fasta("ncds", ncds_dict, FASTA_FILE, CLUST, n_chrom, PRCT_MISS, BPP)
    write_to_bed("ncds", ncds_dict, n_chrom)
    write_to_bed("cds", cds_dict, n_chrom)
