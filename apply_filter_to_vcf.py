#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 12:03:09 2020
@author: Scott T. Small

This script takes a vcf file run on an individual sample, here filtered with the
CNN module of GATK, and transfers the filtered sites to a population VCF called
with genotype vcfs

Example
-------

    $ python apply_filet_to_vcf.py -v VCF.unfiltered -d individual/filteredVcfs \
        -o output_file

Notes
-----
    The names of the samples are expected to the same. The input file should be
    gzipped and the output file will be gzipped. All filtered files for
    individuals should be in the same directory and gzipped.

"""
import sys
import gzip
import glob
from collections import defaultdict
import argparse


def build_filter(filter_dir):
    """Read vcfs and record filtered sites in dict.

    Parameters
    ----------
    filter_dir : TYPE
        DESCRIPTION.

    Returns
    -------
    filt_dict : TYPE
        DESCRIPTION.

    """
    filt_dict = defaultdict(list)
    vcf_files = glob.glob(f"{filter_dir}/*vcf.gz")
    for file in vcf_files:
        with gzip.open(file, 'rb') as f:
            for line in f:
                line = line.decode()
                if line.startswith("##"):
                    pass
                elif line.startswith("#CHROM"):
                    sample = line.split()[-1]
                else:
                    vcf_line = line.split()
                    chrom = vcf_line[0]
                    pos = vcf_line[1]
                    filt = vcf_line[6]
                    if filt == ".":
                        pass
                    else:
                        site = f'{chrom}_{pos}'
                        filt_dict[site].append(sample)
    return filt_dict


def write_filter(filt_dict, output_name):
    """Write individual filters as file.

    Parameters
    ----------
    filt_dict : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    indv_dict = defaultdict(list)
    for site in filt_dict.keys():
        chrom, pos = site.split("_")
        for indv in filt_dict[site]:
            indv_dict[f'{indv}:{chrom}'].append(pos)
    contig_chrom = output_name.split(".")[2]
    with gzip.open(f"KirFol.AfunF3.{contig_chrom}.filter.mask.txt.gz", 'wt') as filt:
        for indv in indv_dict.keys():
            site_list = "\t".join(indv_dict[indv])
            filt.write(f'{indv}\t{site_list}\n')

    return None


def add_filter(vcfFile, output_name, filt_dict):
    """Add filter to vcf and write new.

    Parameters
    ----------
    vcfFile : TYPE
        DESCRIPTION.
    output_name : TYPE
        DESCRIPTION.
    filt_dict : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    with gzip.open(f'{output_name}', 'wt') as out:
        with gzip.open(vcfFile, 'rb') as vcf:
            for line in vcf:
                line = line.decode()
                if line.startswith("##"):
                    out.write(line)
                elif line.startswith("#CHROM"):
                    sample_list = line.split()
                    out.write(line)
                else:
                    vcf_line = line.split()
                    chrom = vcf_line[0]
                    pos = vcf_line[1]
                    site = f"{chrom}_{pos}"
                    masked_samples = filt_dict[site]
                    for indv in masked_samples:
                        s_ix = sample_list.index(indv)
                        vcf_line[s_ix] = "./."
                    vcf_tab = "\t".join(vcf_line)
                    out.write(f"{vcf_tab}\n")
    return None


def parse_args(args_in):
    """Parse args."""
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', "--vcf", type=str, required=True,
                        help="vcf file to apply filter")
    parser.add_argument('-d', "--filter_directory", type=str, required=True,
                        help="directory that contains individual vcfs with filter "
                        "information")
    parser.add_argument('-o', "--output_name", type=str, help="name to give "
                        "output file. gzip will be added")
    return(parser.parse_args(args_in))


def main():
    """Run main function."""
    args = parse_args(sys.argv[1:])
    # =========================================================================
    #  Gather args
    # =========================================================================
    vcfFile = args.vcf
    indv_dir = args.filter_directory
    output_name = args.output_name
    if output_name is None:
        output_name = vcfFile.rstrip(".vcf.gz")
    # =========================================================================
    #  Main executions
    # =========================================================================
    filt_dict = build_filter(indv_dir)
    write_filter(filt_dict, output_name)
    add_filter(vcfFile, output_name, filt_dict)


if __name__ == "__main__":
    main()
