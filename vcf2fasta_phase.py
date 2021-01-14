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
import tqdm
import argparse
import re
import sys


def vcf2fasta(fastaFile, vcfdict, bed_coords, mask_dt):
    """reads in fasta, changes base using dictionary entry. note that fasta
       will be 0 based so position wll be -1
       notes on Biopython: SeqIO.to_dict() which builds all sequences into a
       dictionary and save it in memory
       SeqIO.index() which builds a dictionary without putting the sequences
       in memory
    """
    for name in vcfdict.keys():
        fasta_sequences = SeqIO.parse(fastaFile, 'fasta')
        # fastadict = SeqIO.to_dict(fasta, 'fasta')
        with open(name + ".fasta", 'w') as out_file:
            # for sample, seq in fastadict.items():
            for fasta in fasta_sequences:
                # read in header and sequence
                header, sequence = fasta.id, str(fasta.seq)
                # retrieve position from dictionary of VCF matching header
                seq = list(sequence)  # strings are immutable
                seq2 = list(sequence)
                # add SNPs
                for items in tqdm(vcfdict[name][header]):
                    pos, allele = items
                    if len(allele) > 1:  # phased
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
                # add mask
                if mask_dt:
                    for k in mask_dt[name].keys():
                        seq[k-1] = "N"
                        seq2[k-1] = "N"
                # when done with the header, write
                if bed_coords:
                    with open(bed_coords) as bed:
                        for line in bed:
                            if line.startswith("chrom"):
                                continue
                            else:
                                chrom, start, end = line.split()
                                start = int(start)
                                end = int(end)
                                out_file.write(">{}_0:{}\n{}\n".format(name, f"{chrom}:{start}_{end}", ''.join(seq[start:end])))
                                out_file.write(">{}_1:{}\n{}\n".format(name, f"{chrom}:{start}_{end}", ''.join(seq2[start:end])))
                else:
                    out_file.write(">{}_0:{}\n{}\n".format(name, header, ''.join(seq)))
                    out_file.write(">{}_1:{}\n{}\n".format(name, header, ''.join(seq2)))


def vcfsample(vcf, fastaFile, bed_coords):
    """reads a vcf file and stores info in a dictionary. Here using a tuple
       for the key
    """
    vcfdict = defaultdict(lambda: defaultdict(list))
    with open(vcf, 'r') as vcf:
        for line in vcf:
            if line.startswith("#CHROM"):
                samples = line.strip().split()
            if not line.startswith("#"):
                x = line.strip().split()
                CHR = x[0]
                POS = int(x[1])
                alleles = x[3] + x[4]
                for sample in range(9, len(samples)):
                    gt = re.split("/|\|", x[sample].split(":")[0])
                    ALLELE = ''
                    for al in gt:
                        try:
                            ALLELE += alleles[int(al)]
                        except ValueError:
                            ALLELE += "N"
                    vcfdict[samples[sample]][CHR].append([POS, ALLELE])

    return vcfdict


def mask_dict(mask_file):
    mask_dt = {}
    with open(mask_file, 'r') as mask:
        for line in mask:
            m_lin = line.split()
            id = m_lin[0]
            pos = m_lin[2:]
            mask_dt[id] = {k: '' for k in pos}

    return mask_dt


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
    vcfdict = vcfsample(in_vcf, in_fasta, bed_coords)
    mask_dt = mask_dict(mask_file)
    vcf2fasta(in_fasta, vcfdict, bed_coords, mask_dt)


if __name__ == "__main__":
    main()