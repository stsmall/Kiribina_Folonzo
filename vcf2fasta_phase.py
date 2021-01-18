# -*- coding: utf-8 -*-
"""
Created on Sun Jan 25 11:25:05 2015
@author: stsmall

Take an input VCF for a single sample and makes a diploid consensus
(ambiguities) from a reference fasta file. Note it also checks for instances
of ./. in the VCF file.vcf2fasta will do this with a phased VCF, appropriate
for something like DnaSP. This implementation was entirely designed for PSMC
(Li and Durbin 2011)

usage: python vcf2fasta.py VCF FASTA-REF
"""

from Bio import SeqIO
from collections import defaultdict
import numpy as np
import argparse
import tqdm
import re
import sys


def seq_mask(mask_file, name):
    with open(mask_file) as mf:
        for line in mf:
            if line.startswith(name):
                ma = np.array(line.strip().split()[2:], dtype=np.int) - 1
                break
    return ma


def vcf2fasta(fasta_file, vcfdict, bed_coords, mask_file):
    """Reads in fast, changes base using dictionary entry.

    note that fasta will be 0 based so position wll be -1
    SeqIO.index(fasta_file, 'fasta') which builds a dictionary without putting
    the sequences in memory
    """
    for name in vcfdict.keys():
        fastadict = SeqIO.index(fasta_file, 'fasta')
        with open(name + ".fasta", 'w') as out_file:
            for chrom in fastadict.keys():
                header = fastadict[chrom].id
                sequence = str(fastadict[chrom].seq)   # strings are immutable
                # add mask
                if mask_file:
                    breakpoint()
                    mask_ls = seq_mask(mask_file, name)
                    m_seq = np.char.array(sequence)
                    m_seq[mask_ls] = 'N'
                    seq = list(m_seq)
                    seq2 = list(m_seq)
                else:
                    seq = list(sequence)
                    seq2 = list(sequence)
                # add SNPs
                for items in vcfdict[name]:
                    pos, allele = items
                    assert len(allele) > 1
                    if seq[pos-1] != "N":
                        # first haplotype
                        if allele[0] == "N":
                            allele1 = "N"
                        elif seq[pos-1].islower():  # soft masked ref
                            allele1 = allele[0].lower()
                        else:
                            allele1 = allele[0]
                        seq[pos-1] = allele1
                        # second haplotype
                        if allele[1] == "N":
                            allele2 = "N"
                        elif seq2[pos-1].islower():
                            allele2 = allele[1].lower()
                        else:
                            allele2 = allele[1]
                        seq2[pos-1] = allele2
                # when done with the header, write
                if bed_coords:
                    with open(bed_coords) as bed:
                        for line in bed:
                            if line.startswith("chrom"):
                                pass
                            else:
                                chrom, start, end = line.split()
                                start = int(start)
                                end = int(end)
                                out_file.write(">{}_0:{}\n{}\n".format(name, f"{chrom}:{start}_{end}", ''.join(seq[start:end])))
                                out_file.write(">{}_1:{}\n{}\n".format(name, f"{chrom}:{start}_{end}", ''.join(seq2[start:end])))
                else:
                    out_file.write(">{}_0:{}\n{}\n".format(name, header, ''.join(seq)))
                    out_file.write(">{}_1:{}\n{}\n".format(name, header, ''.join(seq2)))


def vcfsample(vcf_file):
    """Read a vcf file and stores info in a dictionary.

    Here using a tuple for the key
    """
    snp_count = 0
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if not line.startswith("#"):
                snp_count += 1

    progressbar = tqdm.tqdm(total=snp_count, desc="Read VCF", unit='snp')
    vcfdict = defaultdict(list)
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if line.startswith("#CHROM"):
                samples = line.strip().split()
            if not line.startswith("#"):
                progressbar.update(1)
                x = line.strip().split()
                POS = int(x[1])
                alleles = x[3] + x[4]
                for sample in range(9, len(samples)):
                    gt = re.split("/|\|", x[sample].split(":")[0])
                    ALLELE = ''
                    for al in gt:
                        try:
                            ALLELE += alleles[int(al)]  # if missing ./. will be .
                        except ValueError:
                            ALLELE += "N"  # make . an 'N'
                    vcfdict[samples[sample]].append([POS, ALLELE])

    progressbar.close()
    return vcfdict


def parse_args(args_in):
    """Parse args."""
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', "--vcfFile", type=str, required=True,
                        help='path to vcf file')
    parser.add_argument('-f', "--fasta", required=True, type=str,
                        help='path to fasta file')
    parser.add_argument('-B', '--bedfile', type=str,
                        help='bed file of coordinates')
    parser.add_argument('-M', '--maskfile', type=str,
                        help='mask file of coordinates')
    return(parser.parse_args(args_in))


def main():
    """Run main function."""
    args = parse_args(sys.argv[1:])
    # =========================================================================
    #  Gather args
    # =========================================================================
    in_vcf = args.vcfFile
    in_fasta = args.fasta
    mask_file = args.maskfile
    bed_coords = args.bedfile
    # =========================================================================
    #  Main executions
    # =========================================================================
    vcfdict = vcfsample(in_vcf)
    vcf2fasta(in_fasta, vcfdict, bed_coords, mask_file)


if __name__ == "__main__":
    main()
