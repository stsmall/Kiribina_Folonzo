# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 12:47:16 2018

@author: scott
"""
from Bio import SeqIO
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-f', "--fasta", required=True, type=str,
                    help='path to fasta file')
parser.add_argument('-aa', '--ancestral', required=True, type=str,
                    help='ancestral fasta')
args = parser.parse_args()


def vcf2fasta(fastaFile, ancbed):
    fastadict = {}
    with open(ancbed, 'r') as bed:
        for line in bed:
            x = line.split()
            chrom = x[0]
            pos = int(x[1])  #  0 based
            maj_allele = x[3]
            min_allele = x[4]
            maj_prob = float(x[5])
            if maj_prob > .65:
                fastadict[pos] = maj_allele
            else:
                fastadict[pos] = min_allele

    fasta_sequences = SeqIO.parse(fastaFile, 'fasta')
    with open(f"{chrom}.AA.fasta", 'w') as out_file:
        for fasta in fasta_sequences:
            # read in header and sequence
            header, sequence = fasta.id, str(fasta.seq)
            assert header == chrom, "only 1 chrom entry allowed"
            seq = list(sequence)  # strings are immutable
            for pos in fastadict.keys():
                allele = fastadict[pos]
                # replace base w/ allele at pos-1
                if seq[pos].islower():
                    allele = allele.lower()
                seq[pos] = allele
            # when done with the header, write
            out_file.write(">{}\n{}\n".format(header, ''.join(seq)))
    return(None)


if __name__ == "__main__":
    vcf2fasta(args.fasta, args.ancestral)
