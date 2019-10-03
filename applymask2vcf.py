#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 17:26:30 2019
@author: Scott T. Small

This script takes a VCF and a masking file, created by filt_vcf2mask.py and
returns a new VCF where the filtered sites are coded as missing.

Example
-------
The below example takes a gzipped VCF containing mulitple individuals and a
masking file (BAR.mask.txt). --geno tells it to only print genotypes and not
long INFO columns

    $ python applymask2vcf.py -f FOO.pop.vcf.gz -m BAR.mask.txt -o outFile --geno


Notes
-----
    If a sample is in the VCF but not in masking, then it will not be masked
    if the sample is in the mask but not VCF then it will be skipped.


"""
import argparse
from collections import defaultdict
import gzip

parser = argparse.ArgumentParser()
parser.add_argument('-f', "--vcfFile", type=str, help="vcf file gz")
parser.add_argument('-m', "--maskFile", type=str, help="mask file gz")
parser.add_argument('-o', "--outFile", type=str, help="output name")
parser.add_argument("--geno", action="store_true", help="if added then only "
                    "returns genotypes")
args = parser.parse_args()


def mask2dict(mask_file):
    """Translates mask file into dictionary

    Parameters
    ----------
    mask_file : file
        name of the file with masking data

    Returns
    -------
    mask_dict : dict
        dict of masked positions {"1001":["KirFol1", "KirFol2"]}

    """
    chromosome = ""
    mask_dict = defaultdict(list)
    with gzip.open(mask_file, 'rb') as mf:
        for line in mf:
            line = line.decode()
            if line:
                mask_line = line.split()
                sample_name = mask_line[0]
                if chromosome:
                    if chromosome != mask_line[1]:
                        raise Exception("Expects only 1 chromosome, {}".format(line))
                chromosome = mask_line[1]
                mask_pos = mask_line[2:]
                for pos in mask_pos:
                    mask_dict[pos].append(sample_name)

    return(mask_dict)

def applymask(vcf_file, mask_dict, out_file, geno):
    """Masks a vcf with a mask dict

    Parameters
    ----------
    vcf_file : file
        gzipped file of VCF
    mask_dict : dict
        dict of mask file
    out_file : file
        write new vcf
    geno : bool
        if true then only genotype data is printed

    Returns
    -------
    None

    """
    with gzip.open("{}.gz".format(out_file), 'wb') as f:
        with gzip.open(vcf_file, 'rb') as mf:
           for line in mf:
               line = line.decode()
               if line.startswith("##"):
                   f.write("{}".format(line).encode())
               elif line.startswith("#CHROM"):
                   f.write("{}".format(line).encode())
                   sample_ix = line.split()
               else:
                   var_list = line.split()
                   position = var_list[1]
                   try:
                       mask_samples = mask_dict[position]
                       for ms in mask_samples:
                           if ms in sample_ix:
                               ms_ix = sample_ix.index(ms)
                               var_list[ms_ix] = "./."
                       if geno:
                           var_list[8] = "."
                       f.write("{}\n".format("\t".join(var_list)).encode())
                   except KeyError:
                       f.write("{}\n".format("\t".join(var_list)).encode())
    return(None)


if __name__ == "__main__":
    vcf_file = args.vcfFile
    mask_file = args.maskFile
    out_file = args.outFile
    mask_dict = mask2dict(mask_file)
    applymask(vcf_file, mask_dict, out_file, args.geno)
