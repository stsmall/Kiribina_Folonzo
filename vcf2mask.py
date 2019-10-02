#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 13:19:47 2019
@author: Scott T. Small

Example
-------
VCF file with FILTER columns
X 10 ID REF ALT QUAL FILTER

    $ python vcf2mask.py -m foo.mask -o foo.out.mask --keyword PASS CNN2D


Notes
-----
    (1)This script will not read beyond the FILTER column, so make certain that
    enteries with a PASS have a genotype.
    (2) It expects 1 chromosome per file it will exit with an error but will
    still print what is already masked


"""

import argparse
import glob
import os.path
import gzip

parser = argparse.ArgumentParser()
parser.add_argument('-m', "--mask", type=str, help="masking file"
                    "can have filter field or any number but must be in the same"
                    "heading order as a vcf")
parser.add_argument('-o', "--output", type=str, help="name of "
                    "output file, it exists will append")
parser.add_argument("--keyword", type=str, nargs='*', default=[".", "PASS"],
                    required=False, help="list of keywords to condition for the"
                    "filter columns")
args = parser.parse_args()


def readmask(mask_file, flags):
    """Takes a file with a vcf like header and records positions that have
    a value in the FILTER column not in flags, default ["Pass", "."]

    Parameters
    ----------
    mask_file : file
        file with chr, pos etc ...
    flags : list
        list of strings to interpret as passing, defaults ["PASS", "."]

    Returns
    -------
    mask_list : list
        list of masked positions
    sample_name : str
        name of sample taken from header
    chromosome : str
        name of chromosome taken from header


    """
    sample_name = ""
    chromosome = ""
    mask_list = []
    with gzip.open(mask_file, 'rb') as mf:
        for line in mf:
            line = line.decode()
            if line.startswith("#CHROM"):
                header_line = line.split()
                sample_name = header_line[-1]
            elif not line.startswith("#"):
                var_line = line.split()
                if chromosome:
                    if chromosome != var_line[0]:
                        raise Exception("Expects only 1 chromosome")
                        break
                chromosome = var_line[0]
                position = var_line[1]
                flt_col = var_line[6]
                if flt_col not in flags:
                    mask_list.append(position)

    return(sample_name, chromosome, mask_list)


def writemask(sample, chrom, mask_list, output_name):
    """Write out the masking list to a file, append if it exists

    Parameters
    ----------
    mask_list : list
        list of masked positions
    sample_name : str
        name of sample taken from header
    chromosome : str
        name of chromosome taken from header
    output_name : str
        file name to write, if exists will append

    Returns
    -------
    None

    """

    if os.path.isfile(output_name):
        f = open(output_name, 'a')
        f.write("{}\t{}\t{}\n".format(sample, chrom, "\t".join(mask_list)))
        f.close()
    else:
        f = open(output_name, 'w')
        f.write("{}\t{}\t{}\n".format(sample, chrom, "\t".join(mask_list)))
        f.close()

    return(None)


if __name__ == "__main__":
    mask = args.mask
    flags = args.keyword
    mask_file_list = glob.glob('*{}'.format(mask))
    for mask_file in mask_file_list:
        sample, chrom, mask_list = readmask(mask_file, flags)
        writemask(sample, chrom, mask_list, args.output)
