#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 12:04:20 2020
@author: Scott T. Small

Use the est-sfs program to add probabilistic ancestral states to a VCF file.
est-sfs will also produce a unfolded site frequence spectrum.

Example
-------

    $ python polarize_vcf.py -v FOO.vcf.gz -i ingroup -o outgroup1 outgroup2 -h 20 2 2


Notes
-----
    1) requires the allele counts from output of vcftools --counts
    2) vcf files must be gzipped

"""
import sys
import gzip
import argparse


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
    chrom = counts_line[0]
    pos = counts_line[1]
    hap = counts_line[3]
    anc_list = [0, 0, 0, 0]  # A C G T
    # first allele
    ref, ref_count = counts_line[4].split(":")
    bp_ix = bp_order.index(ref)
    anc_list[bp_ix] += int(ref_count)
    # second allele
    alt, alt_count = counts_line[5].split(":")
    bp_ix = bp_order.index(alt)
    anc_list[bp_ix] += int(alt_count)

    site = f'{chrom}_{pos}'
    try:
        ancdict[site].append(anc_list)
    except KeyError:
        pass
    return ancdict, hap


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
    anc_dict = {}
    ingroup = fileIngroup
    outgroups = fileOutgroup

    # get ingroup counts
    with gzip.open(ingroup, 'r') as counts:
        line = next(counts)  # skip header
        for line in counts:
            counts_line = line.split()
            chrom = counts_line[0]
            pos = counts_line[1]
            site = f'{chrom}_{pos}'
            anc_dict[site] = []
            anc_dict, hap = count_allele(anc_dict, counts_line)
    # get outgroup counts
    for file in outgroups:
        with gzip.open(file, 'r') as counts:
            line = next(counts)  # skip header
            for line in counts:
                counts_line = line.split()
                anc_dict, hap = count_allele(anc_dict, counts_line)

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
    contig = list(anc_dict.keys()[0]).split("_")[0]
    with open(f"{contig}.est.infile", 'w') as est:
        for key in anc_dict:
            counts = [",".join(x) for x in anc_dict[key]]
            est.write(f'{" ".join(counts)}\n')
    # create config file
    n_outgroups = len(counts)
    config = open("config.file", 'w')
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
