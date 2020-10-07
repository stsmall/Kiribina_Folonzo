#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 12:04:20 2020
@author: Scott T. Small

Use the est-sfs program to add probabilistic ancestral states to a VCF file.
est-sfs will also produce a unfolded site frequence spectrum.

Example
-------

    $ python est-sfs_format.py -i ingroup -o outgroup1 outgroup2


Notes
-----
    1) requires the allele counts from output of vcftools --counts
    2) count files must be zipped

"""
import sys
import gzip
import argparse
from collections import defaultdict


def count_allele(ancdict, counts_line):
    """Count alleles.

    Parameters
    ----------
    ancdict : TYPE
        DESCRIPTION.
    counts_line : TYPE
        DESCRIPTION.

    Returns
    -------
    ancdict : dict
        dict with building states
    hap : int
        number of haplotypes

    """
    bp_order = ["A", "C", "G", "T"]
    anc_list = [0, 0, 0, 0]  # A C G T

    ref, ref_count = counts_line[4].split(":")
    if len(ref) == 1 and int(ref_count) > 0:
        bp_ix = bp_order.index(ref)
        anc_list[bp_ix] = 1
    else:
        try:
            alt, alt_count = counts_line[5].split(":")
            if len(alt) == 1 and int(alt_count) > 0:
                bp_ix = bp_order.index(alt)
                anc_list[bp_ix] = 1
        except IndexError:
            pass

    return anc_list


def estsfs_format(fileIngroup, fileOutgroup):
    """Read in allele counts for est-sfs input.

    Parameters
    ----------
    fileIngroup : str
        ingroup counts
    fileOutgroup : list
        outgroup counts

    Returns
    -------
    anc_dict : dict
        dictionary of anc sites
    hap_list : list
        returns the max of the haplist

    """
    anc_dict = defaultdict(list)
    ingroup = fileIngroup
    outgroups = fileOutgroup
    anc_dict = {}
    # get ingroup counts
    with gzip.open(ingroup, 'r') as counts:
        line = next(counts)  # skip header
        for line in counts:
            counts_line = line.split()
            chrom = counts_line[0]
            pos = counts_line[1]
            site = f'{chrom}_{pos}'
            counts = count_allele(counts_line)
            anc_dict[site].append(counts)
    # get outgroup counts
    for file in outgroups:
        with gzip.open(file, 'r') as counts:
            line = next(counts)  # skip header
            for line in counts:
                counts_line = line.split()
                chrom = counts_line[0]
                pos = counts_line[1]
                site = f'{chrom}_{pos}'
                if site in anc_dict:
                    counts = count_allele(counts_line)
                    anc_dict[site].append(counts)

    return anc_dict


def estsfs_infiles(anc_dict):
    """Run est-sfs.

    Parameters
    ----------
    anc_dict : dict
        DESCRIPTION.

    Returns
    -------
    None.

    """
    # create input file
    first = next(iter(anc_dict.keys()))
    chrom, pos = first.split("_")
    with open(f"{chrom}.est.infile", 'w') as est:
        for key in anc_dict:
            counts = [",".join(x) for x in anc_dict[key]]
            est.write(f'{" ".join(counts)}\n')
    # create config file
    n_outgroups = len(counts) - 1
    config = open(f"{chrom}.config.file", 'w')
    config.write(f'n_outgroup={n_outgroups}\nmodel 1\nnrandom 1')

    return None


def parse_args(args_in):
    """Parse args."""
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', "--ingroup", type=str, required=True,
                        help="ingroup/focalgroup counts")
    parser.add_argument('-o', "--outgroup", type=str, nargs='+', required=True,
                        help="outgroup counts")
    return(parser.parse_args(args_in))


def main():
    """Run main function."""
    # =========================================================================
    #  Gather args
    # =========================================================================
    args = parse_args(sys.argv[1:])

    fileIngroup = args.ingroup
    fileOutgroup = args.outgroup

    # =========================================================================
    #  Main executions
    # =========================================================================
    anc_dict = estsfs_format(fileIngroup, fileOutgroup)
    estsfs_infiles(anc_dict)


if __name__ == "__main__":
    main()
