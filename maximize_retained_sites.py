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
from tqdm import tqdm
import allel
import copy
import argparse


def max_sites(sample_mask, gt, perc):
    """Minimize fx.

    Returns
    -------
    None.

    """
    gtt = gt[:, sample_mask]
    n_sites = gtt.shape[0]
    # count sites w/out indv
    missing = gtt.is_missing()
    sites = np.sum(missing, axis=1) / n_sites
    # count sites lost at percentage
    p_sites = sites < perc
    np_sites = np.sum(p_sites)
    return sum(sample_mask) * (1/np_sites)


def missing(gt, sample_mask):
    """Count original missing.

    Parameters
    ----------
    gt : TYPE
        DESCRIPTION.
    sample_mask : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    gtt = gt[:, sample_mask]
    n_sites = gtt.shape[0]
    # count sites w/out indv
    missing = gtt.is_missing()
    sites = np.sum(missing, axis=0) / n_sites
    return(sites)


def load_vcf(vcfFile, perc):
    """Load vcf.

    Parameters
    ----------
    vcfFile : TYPE
        DESCRIPTION.

    Returns
    -------
    gt : TYPE
        DESCRIPTION.

    """
    print(f"loading vcf:{vcfFile}")
    callset = allel.read_vcf(vcfFile)
    sample_id = callset['samples']
    gt = allel.GenotypeArray(callset['calldata/GT'])

    sample_mask = np.array([True] * len(sample_id))
    miss_list = missing(gt, sample_mask)
    miss_list_save = copy.copy(miss_list)

    print("running...")
    pbar = tqdm(total=len(sample_mask))
    mx = []
    samplecfg = []
    while sum(sample_mask) > 0:
        mx.append(max_sites(sample_mask, gt, perc))
        samplecfg.append(copy.copy(sample_mask))
        ix = np.where(miss_list == max(miss_list))
        sample_mask[ix] = False
        miss_list[ix] = 0
        pbar.update(1)
    pbar.close()

    mxgrad = np.gradient(mx)
    mxgrad_ix = np.where(mxgrad == max(mxgrad))[0][0]
    final_cfg = sample_id[samplecfg[mxgrad_ix]]
    return final_cfg


def parse_args(args_in):
    """Parse args."""
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v", "--vcfFile", type=str, required=True,
                        help="vcf to iterate for missingness")
    parser.add_argument("-t", "--target", type=float, required=True, default=0.05,
                        help="target of site missingness")
    return(parser.parse_args(args_in))


def main():
    """Run main function."""
    args = parse_args(sys.argv[1:])
    # =========================================================================
    #  Gather args
    # =========================================================================
    vcfFile = args.vcfFile
    target = args.target
    # =========================================================================
    #  Main executions
    # =========================================================================
    sample_cfg = load_vcf(vcfFile, target)
    with open(f"keep_samples.{target}.txt", 'wt') as out:
        for s in sample_cfg:
            out.write(f"{s}\n")
    print(f"{sample_cfg}")
    print(f"{len(sample_cfg)}")


if __name__ == "__main__":
    main()
