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


def count_allele(counts_line, ingroup):
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
    ref_count = int(ref_count)
    if len(ref) == 1 and ref_count > 0:
        bp_ix = bp_order.index(ref)
        if ingroup:
            anc_list[bp_ix] += ref_count
        else:
            anc_list[bp_ix] = 1
    try:
        alt, alt_count = counts_line[5].split(":")
        alt_count = int(alt_count)
        if len(alt) == 1 and alt_count > 0:
            try:
                bp_ix = bp_order.index(alt)
                if ingroup:
                    anc_list[bp_ix] += alt_count
                else:
                    if sum(anc_list) == 0:
                        anc_list[bp_ix] = 1
            except ValueError:
                pass
    except IndexError:
        pass

    return anc_list


def estsfs_format(file_ingroup, file_outgroup):
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
    ingroup = file_ingroup
    outgroups = file_outgroup
    # get ingroup counts
    with gzip.open(ingroup, 'r') as counts:
        line = next(counts)  # skip header
        for line in counts:
            line = line.decode()
            line = line.split()
            chrom = line[0]
            pos = line[1]
            site = f'{chrom}_{pos}'
            anc_counts = count_allele(line, ingroup=True)
            anc_dict[site].append(anc_counts)
    # get outgroup counts
    for file in outgroups:
        with gzip.open(file, 'r') as counts:
            line = next(counts)  # skip header
            for line in counts:
                line = line.decode()
                line = line.split()
                chrom = line[0]
                pos = line[1]
                site = f'{chrom}_{pos}'
                if site in anc_dict:
                    anc_counts = count_allele(line, ingroup=False)
                    anc_dict[site].append(anc_counts)

    return anc_dict


def estsfs_infiles(anc_dict, n_outgroup):
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
    chrom = first.split("_")[0]
    out = open(f"{chrom}.pos.txt", 'w')
    with open(f"{chrom}.est.infile", 'w') as est:
        for key in anc_dict:
            chrom, pos = key.split("_")
            counts = [",".join(map(str, x)) for x in anc_dict[key]]
            while len(counts) < (n_outgroup + 1):
                counts.append('0,0,0,0')
            est.write(f'{" ".join(counts)}\n')
            out.write(f"{chrom}\t{pos}\n")
    out.close()
    # create config file
    n_outgroups = len(counts) - 1
    config = open(f"{chrom}.config.file", 'w')
    config.write(f'n_outgroup={n_outgroups}\nmodel 1\nnrandom 1')
    config.close()


def parse_args(args_in):
    """Parse args."""
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', "--ingroup", type=str, required=True,
                        help="ingroup/focalgroup counts")
    parser.add_argument('-o', "--outgroup", type=str, nargs='+', required=True,
                        help="outgroup counts")
    return parser.parse_args(args_in)


def main():
    """Run main function."""
    # =========================================================================
    #  Gather args
    # =========================================================================
    args = parse_args(sys.argv[1:])

    file_ingroup = args.ingroup
    file_outgroup = args.outgroup

    # =========================================================================
    #  Main executions
    # =========================================================================
    anc_dict = estsfs_format(file_ingroup, file_outgroup)
    estsfs_infiles(anc_dict, len(file_outgroup))


if __name__ == "__main__":
    main()
