#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 09:14:16 2020
@author: Scott T. Small

Reorder VCF to the ancestral allele. The ancestral allele is expected to be
specified in the INFO line as AA=str. Optional to another field with the probability.

*Only reorients the gt not the full FORMAT. So 1/1 will be 0/0 but the remainder
of the fields will not be changed

Example
-------

    $ python vcf2derived_allele.py VCF

"""
import sys
import argparse
import gzip
import tqdm

CHROMDICT = {}
print("GLOBAL declared for chrom lengths, change line 49 or comment to infer from VCF")
CHROMDICT = {"3": 93833897, "2": 99479438, "X": 17661987}


def reset_genotypes(genotypes, allele_index):
    """Reset genotype to ancestral as reference.

    Parameters
    ----------
    genotypes : List
        DESCRIPTION.
    allele_index : Dict
        DESCRIPTION.
    phased : bool
        DESCRIPTION.

    Returns
    -------
    genotypes : List
        DESCRIPTION.

    """
    for i, gt in enumerate(genotypes):
        if ":" in gt:
            gt = gt.split(":")[0]
        assert len(gt) == 3, f"misformed gt field sample {i}, {gt}"
        try:
            gt[0] = allele_index[gt[0]]
            gt[2] = allele_index[gt[2]]
        except KeyError:
            # ./.
            pass
        genotypes[i] = gt
    return genotypes


def repolarize(vcf_file, ancprob=0.90):
    """Reorient VCF to ancestral allele.

    Parameters
    ----------
    vcf_file : str
        DESCRIPTION.
    phased : bool, optional
        DESCRIPTION. The default is False.
    ancprob : float, optional
        DESCRIPTION. The default is 0.90.

    Returns
    -------
    low_count : int
        DESCRIPTION.
    site_count : int
        DESCRIPTION.

    """
    if vcf_file.endswith(".gz"):
        fopen = gzip.open
    else:
        fopen = open
    outfile = open(f"{vcf_file}.derived", 'w')

    # set progress bar and chrom lengths
    with fopen(vcf_file, 'rt') as pgb:
        for line in pgb:
            if line.startswith("##contig"):
                if CHROMDICT:
                    pass
                else:
                    chrom, length = line.split(",")
                    chrom = chrom.split("=")[2]
                    length = int(length.split("=")[1][:-1])
                    CHROMDICT[chrom] = length
            elif not line.startswith("#"):
                v_cols = line.split()
                chrom = v_cols[0]
                pos = int(v_cols[1])
                break
    progressbar = tqdm.tqdm(total=CHROMDICT[chrom], desc="Read VCF", unit='sites')

    # read vcf
    low_count = 0
    site_count = 0
    with fopen(vcf_file, 'rt') as vcf:
        for line in vcf:
            if line.startswith("#"):
                outfile.write(line)
            elif not line.startswith("#"):
                v_cols = line.split()
                chrom = v_cols[0]
                variant_pos = int(v_cols[1])
                progressbar.update(variant_pos - pos)
                variant_pos = pos
                site_count += 1
                # alleles
                ref = v_cols[3]
                assert len(ref) == 1, "indels can not be polarized"
                alt = v_cols[4].split(",")
                alleles = [ref] + alt
                # ancestral state
                INFO = v_cols[7].split(";")
                try:
                    for field in INFO:  # AA
                        if "AA" in field:
                            ancestral = field.split("=")[1]
                        if "AAProb" in field:
                            probability = float(field.split("=")[1])
                            if probability < ancprob:
                                low_count += 1
                except IndexError:
                    print("Improper AA field")
                    sys.exit()
                # reset genotype ID
                if ref != ancestral:
                    # rearrange so derived is ref
                    ordered_alleles = [ancestral] + list(set(alleles) - {ancestral})
                    allele_index = {old_index: ordered_alleles.index(allele)
                                    for old_index, allele in enumerate(alleles)}
                    # reset all genotype fields
                    genotypes = v_cols[8:]
                    genotypes = reset_genotypes(genotypes, allele_index)
                    # replace vcf fields
                    v_cols[3] = ancestral
                    v_cols[4] = ",".join(list(set(alleles)-{ancestral}))
                    v_cols[8:] = genotypes
                    v_cols_line = "\t".join(v_cols)
                    # write back to vcf line
                    outfile.write(f"{v_cols_line}\n")
                else:
                    outfile.write(line)
            else:
                print(line)

    return low_count, site_count


def parse_args(args_in):
    """Parse args."""
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("vcf", type=str, help="VCF type file to reorder")
    return(parser.parse_args(args_in))


def main():
    """Run main function."""
    args = parse_args(sys.argv[1:])
    # =========================================================================
    #  Gather args
    # =========================================================================
    vcf_file = args.vcf
    # =========================================================================
    #  Main executions
    # =========================================================================
    repolarize(vcf_file)


if __name__ == "__main__":
    main()
