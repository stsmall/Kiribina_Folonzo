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
    site_count = 0
    est_dict = {}
    contig = list(anc_dict.keys()[0]).split("_")[0]
    ancFiles = glob.glob(f"{contig}.estsfs*.anc")
    for ancFile in ancFiles:
        with open(ancFile, 'r') as est:
            for line in est:
                if line.startswith('0'):
                    pass
                else:
                    est_line = line.split()
                    est_dict[site_count] = float(est_line[2])
                    site_count += 1

    polar_dict = {}
    for i, site in enumerate(anc_dict.keys()):
        polar_dict[site] = est_dict[str(i)]

    return polar_dict


def polarize_vcf(vcfFile, polar_dict):
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
    pvcf = gzip.open(f"{outfile}.AA.vcf.gz", 'wb')

    ancbed = gzip.open(f"{outfile}.anc.bed.gz", 'w')

    with gzip.open(vcfFile, 'rb') as vcf:
        for line in vcf:
            if line.startswith("#"):
                pvcf.write(line)
            else:
                vcf_line = line.split()
                chrom = vcf_line[0]
                pos = vcf_line[1]
                ref = vcf_line[3]
                alt = vcf_line[4]
                site = f'{chrom}_{pos}'
                aa = vcf_line[8].split(";")
                try:
                    maj_prob = polar_dict[site]
                except KeyError:
                    aa.insert(0, 'AA=NA;MajProb=NA')
                    vcf_line[8] = ";".join(aa)
                    vcf_tab = "\t".join(vcf_line)
                    pvcf.write(f'{vcf_tab}\n')
                    ancbed.write(f'{chrom}\t{pos-1}\t{pos}\tNA\tNA\n')
                    continue
                # AC=X in INFO field
                ac, af, an, *_ = vcf_line[8].split(";")
                c_alt = int(ac.split("=")[1])
                c_ref = int(an.split("=")[1]) - int(ac.split("=")[1])
                if c_ref >= c_alt:
                    maj = ref
                else:
                    maj = alt
                aa.insert(0, f'AA={maj};MajProb={maj_prob}')
                vcf_line[8] = ";".join(aa)
                ancbed.write(f'{chrom}\t{pos-1}\t{pos}\t{maj}\t{maj_prob}\n')
                vcf_tab = "\t".join(vcf_line)
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
    parser.add_argument('-i', "--ingroup", type=str, help="counts file for ingroup")
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
