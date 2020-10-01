#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 11:44:33 2020
@author: Scott T. Small

This script return positions thinned for linkage using scikit-allel. An alternative
to using something like PLINK

Example
-------

    $ python example_numpy.py


"""
import allel
import zarr
import numcodecs
import sys
import numpy as np
import argparse
from dataclasses import dataclass


@dataclass
class AllelData:
    """Class to hold allel data."""

    __slots__ = ["name", "gt", "pos"]
    name: str
    gt: list
    pos: list


def vcf2zarr(chroms, zarr_path, vcf_path):
    """Convert vcf to zarr.

    Parameters
    ----------
    chroms : TYPE
        DESCRIPTION.
    zarr_path : TYPE
        DESCRIPTION.
    vcf_path : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    for c in chroms:
        allel.vcf_to_zarr(vcf_path, zarr_path, group=c,
                          fields='*', alt_number=2, log=sys.stdout,
                          compressor=numcodecs.Blosc(cname='zstd', clevel=1, shuffle=False))
    return None


def load_zarr(chroms, zarr_path):
    """Load zarr to GenotypeArray.

    Parameters
    ----------
    chroms : TYPE
        DESCRIPTION.
    zarr_path : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    callset = zarr.open_group(zarr_path, mode='r')
    chrom_dict = {}
    for c in chroms:
        pos = allel.SortedIndex(callset[f'{c}/variants/POS'])
        gt = allel.GenotypeArray(callset[f'{c}/calldata/GT'])
        chrom_dict[c] = AllelData(c, gt, pos)

    return chrom_dict


def ld_prune(gn, pos, size=500, step=200, threshold=.1, n_iter=5):
    """Remove sites in LD.

    Parameters
    ----------
    gn : TYPE
        DESCRIPTION.
    pos : TYPE
        DESCRIPTION.
    size : TYPE, optional
        DESCRIPTION. The default is 500.
    step : TYPE, optional
        DESCRIPTION. The default is 200.
    threshold : TYPE, optional
        DESCRIPTION. The default is .1.
    n_iter : TYPE, optional
        DESCRIPTION. The default is 5.

    Returns
    -------
    TYPE
        DESCRIPTION.
    gn : TYPE
        DESCRIPTION.

    """
    for i in range(n_iter):
        loc_unlinked = allel.locate_unlinked(gn, size=size, step=step, threshold=threshold)
        n = np.count_nonzero(loc_unlinked)
        n_remove = gn.shape[0] - n
        print(f"iteration {i+1} retaining {n} removing {n_remove} variants")
        gn = gn.compress(loc_unlinked, axis=0)
        pos = pos[loc_unlinked]

    return allel.SortedIndex(pos), gn


def filter_fx(gt, pos):
    """Filter GenotypeArray on singletons and biallelic.

    Parameters
    ----------
    gt : TYPE
        DESCRIPTION.
    pos : TYPE
        DESCRIPTION.

    Returns
    -------
    gf : TYPE
        DESCRIPTION.

    """
    # filter on singleton and >2 allele
    ac = gt.count_alleles()
    print(np.count_nonzero((ac.max_allele() == 1) & ac.is_singleton(1)))
    flt = (ac.max_allele() == 1) & (ac[:, :2].min(axis=1) > 1)
    posflt = pos[flt]
    gf = gt.compress(flt, axis=0)

    return gf, posflt


def get_ldthin_pos(chrom_dict, n_snps, threshold, n_iter, size=500, step=200):
    """Return a list of sites that are unlinked.

    Parameters
    ----------
    chrom_dict : TYPE
        DESCRIPTION.
    n_snps : TYPE
        DESCRIPTION.
    threshold : TYPE
        DESCRIPTION.
    n_iter : TYPE
        DESCRIPTION.

    Returns
    -------
    thinld_dict : TYPE
        DESCRIPTION.

    """
    thinld_dict = {}
    for chrom in chrom_dict.keys():
        pos = chrom_dict[chrom].pos
        gt = chrom_dict[chrom].gt
        gt, pos = filter_fx(gt, pos)  # possibly not needed
        gn = gt.to_n_alt()
        # random downsample to reduce LD
        n = n_snps  # number of SNPs to choose randomly
        vidx = np.random.choice(gn.shape[0], n, replace=False)
        vidx.sort()
        gnr = gn.take(vidx, axis=0)
        posflt = pos[vidx]
        # thinld
        thinpos, gnu = ld_prune(gnr, posflt, size=size, step=step, threshold=threshold, n_iter=n_iter)
        thinld_dict[chrom] = thinpos

    return thinld_dict


def write_ldthin_pos(thinpos_dict):
    """Write positions to file.

    Parameters
    ----------
    thinpos_dict : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    for nchr in thinpos_dict.keys():
        with open(f"{nchr}.thin.pos.txt", 'w') as tpos:
            for pos in thinpos_dict[nchr]:
                tpos.write(f"{nchr}\t{pos}\n")

    return None


def parse_args(args_in):
    """Parse args."""
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--chroms", nargs='+', action="append")
    parser.add_argument("--zarr", type=str, help="where to write zarr")
    parser.add_argument("--vcf", type=str, help="where to get vcf, will use glob")
    parser.add_argument("--n_snp", type=int, default=500000, help="number of snps"
                        " to retain.")
    parser.add_argument("--n_iter", type=int, default=5, help="number of iterations"
                        " of thinning")
    parser.add_argument("--r2", type=float, default=0.1, help="remove snps above"
                        " this correlation.")
    return(parser.parse_args(args_in))


def main():
    """Run main function."""
    args = parse_args(sys.argv[1:])
    # =========================================================================
    #  Gather args
    # =========================================================================
    chroms = args.chroms  # ["2", "3", "X"]
    vcf_path = args.vcf
    zarr_path = args.zarr
    n_snps = args.n_snp
    threshold = args.r2
    n_iter = args.n_iter
    # =========================================================================
    #  Main executions
    # =========================================================================
    vcf2zarr(chroms, zarr_path, vcf_path)
    chrom_dict = load_zarr(chroms, zarr_path)
    thinpos_dict = get_ldthin_pos(chrom_dict, n_snps, threshold, n_iter)
    write_ldthin_pos(thinpos_dict)


if __name__ == "__main__":
    main()
