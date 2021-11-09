#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 12:03:09 2020
@author: Scott T. Small

Filter triallelic sites by checking the allele frequencies. If the 3rd allele
is either a singleton or doubleton, replace as missing and count as a sequencing
error. Return a VCF with those individual genotypes masked and a count/perctage
of errors

Example
-------

    $ python calc_err_from_triallelic.py -v FOO.vcf

"""
import sys
import gzip
import argparse


def add_filter(vcfFile):
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
    with gzip.open(f'{vcfFile}.triallelic', 'wt') as out:
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
                    ref = vcf_line[3]
                    alt = vcf_line[4]
                    fields = vcf_line[7]
                    if len(ref) > 2:
                        print(vcf_line)
                        continue
                    if len(alt) > 1:
                        AC = fields.split(";")[-1].split("=")[1].split(",")
                        AC = map(int, AC)
                        if AC[-1] < 3:
                            if AC[-1] == 1:
                                idx = [i+9 for i, j in enumerate(vcf_line[9:]) if "2" in j.split(":")[0]]
                                vcf_line[idx] = "./."
                            elif AC[-1] == 2:
                                idx = [i+9 for i, j in enumerate(vcf_line[9:]) if j.split(":")[0].count("2") == 2]
                                if idx is not None:
                                    vcf_line[idx] = "./."
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
                        "output file")
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
