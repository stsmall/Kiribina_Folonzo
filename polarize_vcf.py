#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 13:44:12 2018

@author: scott
"""

import argparse
import gzip
import sys

def read_estsfs(ingroup, estOut):
    """Read est-sfs outfiles.

    Parameters
    ----------
    anc_dict : TYPE
        DESCRIPTION.

    Returns
    -------
    polat_dict : TYPE
        DESCRIPTION.

    """
    est_dict = {}
    with open(estOut, 'r') as est:
        for line in est:
            line = line.split()
            chrom = line[0]
            pos = line[1]
            prob_maj = line[2]
            est_dict[f"{chrom}_{pos}"] = [prob_maj]

    with open(ingroup, 'r') as counts:
        for line in counts:
            line = line.split()
            # find maj count
            # try to add to est_dict
            # if not exists then add as [0.0, maj]

    return est_dict


def polarize_vcf(vcfFile, est_dict):
    """Add AA and MajProb to INFO column.

    Parameters
    ----------
    vcfFile ; str
        vcf file to add AA
    polar_dict : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    outfile = vcfFile.rstrip("vcf.gz")
    pvcf = gzip.open(f"{outfile}.anc.vcf.gz", 'wb')
    for k in est_dict.keys():
        chrom, pos = k.split("_")
        break
    ancbed = gzip.open(f"{chrom}.anc.bed.gz", 'wb')

    with gzip.open(vcfFile, 'rb') as vcf:
        for line in vcf:
            if line.startswith("#"):
                pvcf.write(line)
            else:
                line = line.split()
                chrom = line[0]
                pos = line[1]
                site = f'{chrom}_{pos}'
                aa = line[8].split(";")
                maj_prob, maj_allele = est_dict[site]
                aa.insert(0, f'AA={maj_allele};AAProb={maj_prob}')
                line[8] = ";".join(aa)
                ancbed.write(f'{chrom}\t{pos-1}\t{pos}\t{maj_allele}\t{maj_prob}\n')
                vcf_tab = "\t".join(line)
                pvcf.write(f'{vcf_tab}\n')
    pvcf.close()
    ancbed.close()

    return None


def parse_args(args_in):
    """Parse args."""
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', "--vcfFile", type=str, required=True,
                        help="vcf file")
    parser.add_argument('-i', "--ingroup", type=str, help="counts file for "
                        "ingroup")
    parser.add_argument('-e', "--estFile", type=str, required=True,
                        help="est-sfs output with 1st and 2nd columns having CHR"
                        "POS")
    return(parser.parse_args(args_in))


def main():
    """Run main function."""
    args = parse_args(sys.argv[1:])
    # =========================================================================
    #  Gather args
    # =========================================================================
    estFile = args.estFile
    vcfFile = args.vcfFile
    counts = args.ingroup
    # =========================================================================
    #  Main executions
    # =========================================================================
    polar_dict = read_estsfs(estFile, counts)
    polarize_vcf(vcfFile, polar_dict)


if __name__ == "__main__":
    main()
