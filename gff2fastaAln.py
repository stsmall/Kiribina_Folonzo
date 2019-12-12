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
from tqdm import tqdm
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
    print(f"calculating non-cds locations\n")
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
    return(non_cds)


def write_outfile(s_ix: int,
                  e_ix: int,
                  fname: str,
                  chrom: str,
                  bpp: bool,
                  header_list,
                  loci_list,
                  just: int = 10):
    """Writes outfile for format fasta

    Parameters
    ----------
    s_ix: int
        start coordinates
    e_ix: int
        end coordinates
    fname: str
        cds or ncds
    chrom: str
        chromosome name
    bpp: bool
        format for BPP True or False
    header_list: List[List[str]]
        lists of headers from fasta for each locus
    loci_list: List[List[str]]
        list of sequences from fasta for each locus
    just: int
        justified for BPP

    Returns
    -------
    None

    """
    with open(f"{fname}.{chrom}.{s_ix}-{e_ix}.{len(loci_list)}.txt", 'w') as out_file:
        for item in zip(header_list, loci_list):
            headers, seqs = item
            if bpp is True:
                out_file.write(f"{len(headers)} {len(seqs[0])}\n\n")
            for sample, dna in zip(headers, seqs):
                if bpp is True:
                    out_file.write(f"^{sample}{' '*(just-len(sample))}{dna}\n\n")
                else:
                    out_file.write(f">{sample}\n{dna}\n")
    return(None)


def get_fastaseq(fasta_sequences,
                 locus_k):
    """Retrieves sequence from fasta file by coordinates

    Parameters
    ----------
    fasta_sequences: obj of SeqIO
        contains fasta sequences and headers
    locus_k: obj of dataclass
        contains start and end coordinates for loci

    Returns
    -------
    header_l: list[str]
        list of headers from fasta
    loci_l: list[str]
        list of sequences from fasta

    """
    header_l = []
    loci_l = []
    for fasta in fasta_sequences:
        header, sequence = fasta.id, str(fasta.seq)
        loci_l.append(sequence[locus_k.start:locus_k.end])
        header_l.append(header)
    return(header_l, loci_l)


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
    total_loci = len(gff_dict)
    pbar = tqdm(total=total_loci)
    loci_list = []
    header_list = []
    k_dict = {}
    print(f"\n{fname}: formatting files from alignments\n")
    while loci < total_loci:
        while len(loci_list) < clust:
            pbar.update(1)
            k = f"{fname}_{str(loci)}"
            try:
                header_l, loci_l = get_fastaseq(fasta_sequences, gff_dict[k])
            except KeyError:
                k = f"{fname}_{str(loci-1)}"
                break
            seqlen = len(loci_l[0])
            # Ns check point
            miss_list = [(seqX.count("N")/seqlen) for seqX in loci_l]
            if any((seqX.count("N")/seqlen) > prct for seqX in loci_l):
                skip_gaps += 1
            else:
                loci_list.append(loci_l)
                header_list.append(header_l)
                if len(loci_list) == 1:
                    s_ix = gff_dict[k].start
            k_dict[k] = max(miss_list)
            loci += 1
        else:
            e_ix = gff_dict[k].end
            write_outfile(s_ix, e_ix, fname, chrom, bpp, header_list, loci_list)
            loci_list = []
            header_list = []
    if len(loci_list) > 0:
        e_ix = gff_dict[k].end
        write_outfile(s_ix, e_ix, fname, chrom, bpp, header_list, loci_list)
        loci_list = []
        header_list = []
    pbar.close()
    print(f"\n{fname}: {skip_gaps} regions skipped due to excess N's\n")
    return(k_dict)


def write_to_bed(fname: str,
                 gff_dict: Dict[str, object],
                 chrom: int,
                 kdict: Dict[str, float]):
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
        print(f"writing to bed files\n")
        for i in range(len(gff_dict)):
            k = f"{fname}_{str(i)}"
            start = gff_dict[k].start
            end = gff_dict[k].end
            out_bed.write(f"{chrom}\t{start}\t{end}\t{kdict[k]:.2f}\n")
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
    # CDS
    cds_dict, n_chrom = get_cds(GFF_FILE, MIN_LEN_CDS)
    kdict = format_fasta("cds", cds_dict, FASTA_FILE, CLUST, n_chrom, PRCT_MISS, BPP)
    write_to_bed("cds", cds_dict, n_chrom, kdict)
    # Non-CDS
    ncds_dict = get_ncds(cds_dict, MAX_LEN, MIN_LEN, DIST_BETW, CHROM_LEN)
    kdict = format_fasta("ncds", ncds_dict, FASTA_FILE, CLUST, n_chrom, PRCT_MISS, BPP)
    write_to_bed("ncds", ncds_dict, n_chrom, kdict)
