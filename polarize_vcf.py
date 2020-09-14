#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 12:04:20 2020
@author: Scott T. Small

Use the est-sfs program to add probabilistic ancestral states to a VCF file.
est-sfs will also produce a unfolded site frequence spectrum.

Example
-------

    $ python polarize_vcf.py -v FOO.vcf.gz -i ingroup -o outgroup1 outgroup2 -h 20 2 2


Notes
-----
    1) requires the allele counts from output of vcftools --counts
    2) vcf files must be gzipped

"""
import sys
import os
from path import join
import gzip
import glob
from subprocess import run
import argparse


def count_allele(ancdict, counts_line):
    """Count alleles.

    Parameters
    ----------
    ancdict : TYPE
        DESCRIPTION.
    counts_line : TYPE
        DESCRIPTION.

    Returns
    -------
    ancdict : dict
        dict with building states
    hap : int
        number of haplotypes

    """
    bp_order = ["A", "C", "G", "T"]
    chrom = counts_line[0]
    pos = counts_line[1]
    hap = counts_line[3]
    anc_list = [0, 0, 0, 0]  # A C G T
    # first allele
    ref, ref_count = counts_line[4].split(":")
    bp_ix = bp_order.index(ref)
    anc_list[bp_ix] += int(ref_count)
    # second allele
    alt, alt_count = counts_line[5].split(":")
    bp_ix = bp_order.index(alt)
    anc_list[bp_ix] += int(alt_count)

    site = f'{chrom}_{pos}'
    try:
        ancdict[site].append(anc_list)
    except KeyError:
        pass
    return ancdict, hap


def estsfs_format(fileIngroup, fileOutgroup):
    """Read in allele counts for est-sfs input.

    Parameters
    ----------
    fileIngroup : str
        ingroup counts
    fileOutgroup : list
        outgroup counts

    Returns
    -------
    anc_dict : dict
        dictionary of anc sites
    hap_list : list
        returns the max of the haplist

    """
    anc_dict = {}
    hap_list = []
    ingroup = fileIngroup
    outgroups = fileOutgroup

    with open(ingroup, 'r') as counts:
        line = next(counts)  # skip header
        for line in counts:
            counts_line = line.split()
            chrom = counts_line[0]
            pos = counts_line[1]
            site = f'{chrom}_{pos}'
            anc_dict[site] = []
            anc_dict, hap = count_allele(anc_dict, counts_line)
            hap_list.append(hap)

    for file in outgroups:
        with open(file, 'r') as counts:
            line = next(counts)  # skip header
            for line in counts:
                counts_line = line.split()
                anc_dict, hap = count_allele(anc_dict, counts_line)

    return anc_dict, max(hap_list)


def run_estsfs(anc_dict, estsfs_dir, haps):
    """Run est-sfs.

    Parameters
    ----------
    anc_dict : dict
        DESCRIPTION.

    Returns
    -------
    None.

    """
    contig = list(anc_dict.keys()[0]).split("_")[0]
    with open(f"{contig}.est.infile", 'w') as est:
        for key in anc_dict:
            counts = [",".join(x) for x in anc_dict[key]]
            est.write(f'{" ".join(counts)}\n')

    n_outgroups = len(counts)
    config = open("config.file", 'w')
    config.write(f'n_outgroup={n_outgroups}\nmodel 1\nnrandom 1')

    cmd = f'split -l 100000 {contig}.est.infile {contig}.estsfs'
    proc = run(cmd, shell=True, encoding='ascii', check=True)

    estsfs_exe = os.path.join(estsfs_dir, "est-sfs")
    infiles = glob.glob(f"{contig}.estsfs*")
    for infile in infiles:
        cmd = f"{estsfs_exe} config.file {infile} seed.file {infile}.sfs {infile}.anc"
        proc = run(cmd, shell=True, encoding='ascii', check=True)
    return None


def read_estsfs_out(anc_dict):
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
                        help="vcfFile to add AA column")
    parser.add_argument('-i', "--ingroup", type=str, required=True,
                        help="ingroup/focalgroup counts")
    parser.add_argument('-h', "--haplotypes", type=int, nargs='+',
                        action="append", required=True,
                        help="number of haplotypes in focal group")
    parser.add_argument('-o', "--outgroup", type=str, nargs='+', required=True,
                        help="outgroup counts")
    parser.add_argument("--estsfs_exe", type=str,
                        help="path to executable for est-sfs")

    return(parser.parse_args(args_in))


def main():
    """Run main function."""
    # =========================================================================
    #  Gather args
    # =========================================================================
    args = parse_args(sys.argv[1:])

    fileIngroup = args.ingroup
    fileOutgroup = args.outgroup
    haplist = args.haplotypes
    estsfs_dir = args.estsfs_exe
    assert len(haplist) == (len(fileOutgroup) + 1)
    vcfFile = args.vcfFile
    # =========================================================================
    #  Main executions
    # =========================================================================
    # TODO: add a check to filter on number of haplotypes if est-sfs doesnt take missing
    anc_dict, haps = estsfs_format(fileIngroup, fileOutgroup)
    run_estsfs(anc_dict, estsfs_dir, haps)
    polar_dict = read_estsfs_out(anc_dict)
    polarize_vcf(vcfFile, polar_dict)


if __name__ == "__main__":
    main()
