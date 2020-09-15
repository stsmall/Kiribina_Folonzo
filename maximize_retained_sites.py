#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 12:03:48 2020
@author: Scott T. Small

Many programs do not allow for missing data. Keeping highly missing
individual results in regions of false negatives (for hets). We can remove a sub-
set of individuals to minimize these regions. This script tries to find a way to
maximize both the number of sites and retained individual.

Example
-------

    $ python maximize_retained_indvs.py -m FOO.miss.indv.gz -t 0.10 \
        [-u FOO.uncalled.gz] -f FOO.filtered.gz

Notes
-----
    requires four files:
            KirFol.AfunF3.{chroms}.filter.mask.txt.gz from apply_filter_to_vcf.py
            KirFol.AfunF3.{chroms}.uncalled.mask.txt.gz from get_uncalled.py
            KirFol.AfunF3.{chroms}.miss.site.txt.gz  #  prior on missingness

    1) remove highly missing individuals by using --missing-indv and cutoff
    2) run vcftools site missingness : KirFol.AfunF3.{chroms}.miss.site.txt
    3) run apply_filter_to_vcf.py and get_uncalled.py
    4) try to maximize sites and individuals within target, default 5%.


"""
import sys
import numpy as np
import gzip
import argparse


def make_missarray(vcfFile, uncalled, filtered):
    """Load missing data to dictionary.

    Parameters
    ----------
    uncalled : TYPE
        DESCRIPTION.
    filtered : TYPE
        DESCRIPTION.

    Returns
    -------
    miss_dict : dict
        dictionary of missing sites. '5':['indv1', 'indv2', 'indv3']
    """
    pos_list = []
    with gzip.open(vcfFile, 'rb') as vcf:
        for line in vcf:
            line = line.decode()
            if line.startwith("##"):
                pass
            elif line.startswith("#CHROM"):
                samples = line.split()[9:]
            else:
                vcf_list = line.split()
                chrom = vcf_list[0]
                pos = vcf_list[1]
                pos_list.append(pos)

    rows = len(samples)
    columns = len(pos_list)
    pos_arr = np.array(pos_list)
    miss_arr = np.ones([rows, columns])

    with gzip.open(filtered, 'rb') as filt:
        for ix, line in enumerate(filt):
            line = line.decode()
            indv, chrom, *sites = line.split()
            site_arr = np.array(sites)
            miss_arr[ix, ] = np.in1d(pos_arr, site_arr)

    return miss_arr


def check_missingness(miss_arr, missSite):
    """Check that missing values are consistent.

    Parameters
    ----------
    miss_arr : TYPE
        DESCRIPTION.
    missSite : TYPE
        DESCRIPTION.

    Returns
    -------
    site_missing : TYPE
        DESCRIPTION.
    indv_missing : TYPE
        DESCRIPTION.

    """
    site_missing = np.sum(miss_arr, axis=0) / miss_arr.shape[0]
    indv_missing = np.sum(miss_arr, axis=1) / miss_arr.shape[1]

    with gzip.open(missSite, 'rb') as miss:
        for line in miss:
            line = line.decode()

    # TODO: check that these make sense
    return site_missing, indv_missing


def iterate_target_missingness(miss_arr, target):
    # TODO: set up minimize or something
    keep = open("keep.txt", 'w')

    return None


def parse_args(args_in):
    """Parse args."""
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v", "--vcf", type=str, required=True,
                        help="vcf to iterate for missingness")
    parser.add_argument("-t", "--target", type=float, required=True, default=0.05,
                        help="target of site missingness")
    parser.add_argument("-u", "--uncalled", type=str,
                        help="file records of uncalled sites")
    parser.add_argument("-f", "--filtered", type=str, required=True,
                        help="file records of filtered sites")
    parser.add_argument("-m", "--miss", type=str,
                        help="output of vcftools --missing-site")
    return(parser.parse_args(args_in))


def main():
    """Run main function."""
    args = parse_args(sys.argv[1:])
    # =========================================================================
    #  Gather args
    # =========================================================================
    vcfFile = args.vcfFile
    miss = args.miss
    target = args.target
    uncalled = args.uncalled
    filt = args.filtered
    # =========================================================================
    #  Main executions
    # =========================================================================
    miss_arr = make_missarray(vcfFile, uncalled, filt)
    site, indv = check_missingness(miss_arr, miss)
    iterate_target_missingness(miss_arr, target)


if __name__ == "__main__":
    main()
