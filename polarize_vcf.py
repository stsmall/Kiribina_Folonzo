#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 12:04:20 2020
@author: Scott T. Small

Use the est-sfs program to add probabilistic ancestral states to a VCF file.
est-sfs will also produce a unfolded site frequence spectrum.

Example
-------

    $ python polarize_vcf.py -v VCF -i ingroup -e estsfsFile


Notes
-----
    1) requires the allele counts from output of vcftools --counts
    2) count files must be zipped

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
            site = f'{chrom}_{pos}'
            est_dict[site] = [prob_maj]

    with gzip.open(ingroup, 'r') as counts:
        line = next(counts)  # skip header
        for line in counts:
            line = line.decode()
            line = line.split()
            chrom = line[0]
            pos = line[1]
            site = f'{chrom}_{pos}'
            ref, ref_count = line[4].split(":")
            alt, alt_count = line[5].split(":")
            if int(ref_count) > int(alt_count):
                maj = ref
            else:
                maj = alt
            try:
                est_dict[site].extend(maj)
            except KeyError:
                est_dict[site] = [0.0, maj]

    return est_dict


def addAA2vcf(vcfFile, est_dict):
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
    pvcf = gzip.open(f"{outfile}.anc.vcf.gz", 'wt')

    first = next(iter(est_dict.keys()))
    chrom, pos = first.split("_")
    ancbed = gzip.open(f"{chrom}.anc.bed.gz", 'wt')

    with gzip.open(vcfFile, 'r') as vcf:
        for line in vcf:
            line = line.decode()
            if line.startswith("#"):
                if line.startswith("##FORMAT"):
                    pvcf.write('##INFO=<ID=AA,Number=1,Type=String,Description="Anc Allele">')
                    pvcf.write('##INFO=<ID=AAProb,Number=A,Type=Float,Description="Prob Maj is Anc">')
                pvcf.write(line)
            else:
                line = line.split()
                chrom = line[0]
                pos = line[1]
                pos = int(pos)
                ref = line[3]
                alt = line[4]
                site = f'{chrom}_{pos}'
                aa = line[7].split(";")
                maj_prob, maj_allele = est_dict[site]
                if maj_allele == ref:
                    min_allele = alt
                elif maj_allele == alt:
                    min_allele = ref
                if maj_prob > 0.60:
                    aa.insert(0, f'AA={maj_allele};AAProb={maj_prob}')
                else:
                    aa.insert(0, f'AA={min_allele};AAProb={1-maj_prob}')
                line[7] = ";".join(aa)
                ancbed.write(f'{chrom}\t{pos-1}\t{pos}\t{maj_allele}\t{min_allele}\t{maj_prob}\n')
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
    ingroup_counts = args.ingroup
    # =========================================================================
    #  Main executions
    # =========================================================================
    est_dict = read_estsfs(ingroup_counts, estFile)
    addAA2vcf(vcfFile, est_dict)
    #polarize_vcf()


if __name__ == "__main__":
    main()
