#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 10:59:28 2019
@author: Scott T. Small

This program creates a masking file that is in matrix form from a vcf with
missing data. It will also create a bed file of sites that are missing from
all samples.

Example
-------

    $ python missmask.py foo.vcf.gz


Notes
-----
    This should work on phased and unphased genotypes.

"""
import argparse
from collections import defaultdict
import re
import gzip
parser = argparse.ArgumentParser()
parser.add_argument('-f',"--INvcf", type=str,
                    help='path to vcf IN file')
parser.add_argument("--filter", action="store_true")
args = parser.parse_args()


def miss_mask(vcfFile, IX=9):
    """Creates a matrix file with missing data

    Parameters
    ----------
    vcfFile : file
        vcf type file
    IX : int
        column with first sample id
    Returns
    -------
    mask_dict : dict
        dict of masked sites
    chrom : str
        chromosome of vcf

    """
    mask_dict = defaultdict(list)
    chrom = ""
    with open("Genome.mask.bed", 'w') as fout:
        with gzip.open(vcfFile, 'rb') as vcf:
            for line in vcf:
                line = line.decode()
                if "#CHROM" in line:
                    indv_list = line.split()[IX:]
                else:
                    if not line.startswith("#"):
                        var_list = line.split()
                        if chrom:
                            if chrom != var_list[0]:
                                raise Exception("Expects only 1 chromosome, {}".format(line))
                        chrom = var_list[0]
                        pos = int(var_list[1])
                        filt = var_list[6]
                        miss = [i for i, s in enumerate(var_list[IX:]) if re.search(r'\./|/\.|\.\||\|\.', s)]
                        if len(miss) == len(indv_list):
                            fout.write("{}\t{}\t{}\n".format(chrom, pos-1, pos))
                        elif filter:
                            if filt not in ["PASS", "."]:
                                if var_list[4] != ".":
                                    fout.write("{}\t{}\t{}\t{}\n".format(chrom, pos-1, pos, filt))
                        else:
                            for gt_ix, gt in enumerate(var_list[IX:]):
                                sample = indv_list[gt_ix]
                                if "." in gt.split(":")[0]:
                                    mask_dict[sample].append(str(pos))
    return(mask_dict, chrom)


if __name__ == '__main__':
    mask_dict, chrom = miss_mask(args.INvcf)
    with gzip.open("Individual.mask.txt.gz", 'wt') as f:
        for ind in mask_dict.keys():
            f.write("{}\t{}\t{}\n".format(ind, chrom, "\t".join(mask_dict[ind])))
