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
import argparse
import random
import re
import sys


def IUPAC(ALLELE):
    """ALLELE needs to be transformed into IUPAC
    """
    # IUPAC transformation
    if "A" in ALLELE and "G" in ALLELE:
        ALLELE = "R"
    elif "C" in ALLELE and "T" in ALLELE:
        ALLELE = "Y"
    elif "G" in ALLELE and "C" in ALLELE:
        ALLELE = "S"
    elif "A" in ALLELE and "T" in ALLELE:
        ALLELE = "W"
    elif "G" in ALLELE and "T" in ALLELE:
        ALLELE = "K"
    elif "A" in ALLELE and "C" in ALLELE:
        ALLELE = "M"
    else:
        ALLELE = "N"
    return(ALLELE)


def vcf2fasta(fastaFile, vcfdict, phased, NaRef, bed_coords):
    """reads in fasta, changes base using dictionary entry. note that fasta
       will be 0 based so position wll be -1
       notes on Biopython: SeqIO.to_dict() which builds all sequences into a
       dictionary and save it in memory
       SeqIO.index() which builds a dictionary without putting the sequences
       in memory
    """

    for name in vcfdict.keys():
        fasta_sequences = SeqIO.parse(fastaFile, 'fasta')
        # fastadict = SeqIO.to_dict(fasta, 'fasta'
        with open(name + ".fasta", 'w') as out_file:
            # for sample, seq in fastadict.items():
            for fasta in fasta_sequences:
                # read in header and sequence
                header, sequence = fasta.id, str(fasta.seq)
                # retrieve position from dictionary of VCF matching header
                seq = list(sequence)  # strings are immutable
                if phased:
                    seq2 = list(seq)
                for items in vcfdict[name][header]:
                    pos, allele = items
                    if phased and len(allele) > 1:  # phased
                        # first haplotype
                        if allele[0] == "N":
                            if NaRef:
                                if seq[pos-1] != "N":
                                    allele1 = seq[pos-1]
                            else:
                                allele1 = "N"
                        elif seq[pos-1].islower():  # soft masked ref
                            allele1 = allele[0].lower()
                        else:
                            allele1 = allele[0]
                        seq[pos-1] = allele1
                        # second haplotype
                        if allele[1] == "N":
                            if NaRef:
                                if seq2[pos-1] != "N":
                                    allele2 = seq2[pos-1]
                            else:
                                allele2 = "N"
                        elif seq2[pos-1].islower():
                            allele2 = allele[1].lower()
                        else:
                            allele2 = allele[1]
                        seq2[pos-1] = allele2
                    else:
                        if allele == "N":
                            if NaRef:
                                if seq[pos-1] != "N":
                                    allele = seq[pos-1]
                            else:
                                allele = "N"
                        elif seq[pos-1].islower():
                            allele = allele.lower()
                        seq[pos-1] = allele
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
                                if phased:
                                    out_file.write(">{}_0:{}\n{}\n".format(name, f"{chrom}:{start}_{end}", ''.join(seq[start:end])))
                                    out_file.write(">{}_1:{}\n{}\n".format(name, f"{chrom}:{start}_{end}", ''.join(seq2[start:end])))
                                else:
                                    out_file.write(">{}:{}\n{}\n".format(name, f"{chrom}:{start}_{end}", ''.join(seq[start:end])))
                else:
                    if phased:
                        out_file.write(">{}_0:{}\n{}\n".format(name, header, ''.join(seq)))
                        out_file.write(">{}_1:{}\n{}\n".format(name, header, ''.join(seq2)))
                    else:
                        out_file.write(">{}:{}\n{}\n".format(name, header, ''.join(seq)))



def vcfsample(vcf, fastaFile, NaRef, iupac, major, randallele, phased, bed_coords):
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
                for sample in range(9, len(samples)):
                    gt = re.split("/|\|", x[sample].split(":")[0])
                    if len(gt) < 2:
                        gt.append(gt[0])
                    if "," in x[4]:
                        alt = x[4].split(",")
                        alleles = x[3] + alt[0] + alt[1]
                        ALLELE = ''
                        for al in gt:
                            try:
                                ALLELE += alleles[int(al)]
                            except ValueError:
                                if NaRef:
                                    ALLELE += x[3]
                                else:
                                    ALLELE += "N"
                    else:
                        alleles = x[3] + x[4]
                        ALLELE = ''
                        for al in gt:
                            try:
                                ALLELE += alleles[int(al)]
                            except ValueError:
                                if NaRef:
                                    ALLELE += x[3]
                                else:
                                    ALLELE += "N"
                    if ALLELE[0] != ALLELE[1]:
                        if iupac:
                            ALLELE = IUPAC(ALLELE)  # returns single base
                        elif major:  # returns single base
                            # find AD field
                            af = map(int, gt[1].split(","))
                            # max of AD
                            af_ix = af.index(max(af))
                            # stupid tri-allelic
                            if len(af) > 2:
                                if af_ix == 0:
                                    ALLELE = x[3]
                                elif af_ix == 1:
                                    ALLELE = x[4].split(",")[0]
                                else:
                                    ALLELE = x[4].split(",")[1]
                            else:
                                ALLELE = x[af_ix + 3]
                        else:
                            if randallele:  # returns single base
                                ALLELE = random.choice(ALLELE)
                            elif phased:
                                pass  # keep both bases
                    else:
                        if phased:
                            pass  # keep both bases
                        else:
                            ALLELE = ALLELE[0]  # returns single base
                    vcfdict[samples[sample]][CHR].append([POS, ALLELE])
    vcf2fasta(fastaFile, vcfdict, phased, NaRef, bed_coords)
    return(None)


def vcfAAsample(vcf, fasta):
    """ a specific iteration for making an ancestral fasta from a vcf output
        by Maffilter vcf.
    """
    vcf_sample = defaultdict(lambda: defaultdict(list))
    with open(vcf, 'r') as vcf:
        for line in vcf:
            if "#" in line:
                sample = line.strip()[1:]
            if "#" not in line:
                x = line.split()
                CHR = x[0]
                POS = int(x[1])
                allele = x[2]
                vcf_sample[sample][CHR].append([POS, allele])
    vcf2fasta(fasta, vcf_sample)


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
    parser.add_argument('-aa', '--ancestral', action='store_true',
                        help='assume ancestral vcf')
    parser.add_argument('-N', '--NaRef', action='store_true',
                        help='add Ref for missing, default is N')
    parser.add_argument('-U', '--iupac', action='store_true',
                        help='add IUPAC base for het')
    parser.add_argument('-M', '--major', action='store_true',
                        help='call major allel if het based on coverage')
    parser.add_argument('-R', '--rand', action='store_true',
                        help='choose random allele at hets')
    parser.add_argument('-P', '--phased', action='store_true',
                        help='output 2 instead of 1 seq per ind sample')
    return(parser.parse_args(args_in))


def main():
    """Run main function."""
    args = parse_args(sys.argv[1:])
    # =========================================================================
    #  Gather args
    # =========================================================================
    ancestral = args.ancestral
    in_vcf = args.vcfFile
    in_fasta = args.fasta
    na_ref = args.NaRef
    iupac = args.iupac
    major = args.major
    rand = args.rand
    phased = args.phased
    bed_coords = args.bedfile
    # =========================================================================
    #  Main executions
    # =========================================================================
    if ancestral:
        vcfAAsample(in_vcf, in_fasta)
    else:
        if args.phased:
            assert phased > (iupac + major + rand)
        vcfsample(in_vcf, in_fasta, na_ref, iupac, major, rand, phased, bed_coords)


if __name__ == "__main__":
    main()