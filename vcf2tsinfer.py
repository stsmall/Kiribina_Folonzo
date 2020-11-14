#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 17:53:34 2020
@author: Scott T. Small

This code follows the example at:
    https://tsinfer.readthedocs.io/en/latest/tutorial.html
for loading data from a VCF into tsinfer

Example
-------

    $ python example_numpy.py

Notes
-----
    This is an example of an indented section. It's like any other section,
    but the body is indented to help it stand out from surrounding text.

If a section is indented, then a section break is created by
resuming unindented text.

"""

import cyvcf2
import tsinfer
import tqdm
import pandas as pd
import sys
import argparse


def add_metadata(vcf, samples, meta, label_by):
    """Add tsinfer meta data.

    Parameters
    ----------
    vcf : TYPE
        DESCRIPTION.
    samples : TYPE
        DESCRIPTION.
    meta : TYPE
        DESCRIPTION.
    label_by : TYPE, optional
        DESCRIPTION. The default is "Karyotype".

    Returns
    -------
    None.

    """
    pop_lookup = {}
    for pop in meta[label_by].unique():
        pop_lookup[pop] = samples.add_population(metadata={label_by: pop})
    for indiv in vcf.samples:
        meta_dict = meta.loc[indiv].to_dict()
        pop = pop_lookup[meta.loc[indiv][label_by]]
        samples.add_individual(ploidy=2, metadata=meta_dict, population=pop)


def add_diploid_sites(vcf, samples):
    """Read the sites in the vcf and add them to the samples object.

    Reordering the alleles to put the ancestral allele first, if it is available.

    Parameters
    ----------
    vcf : TYPE
        DESCRIPTION.
    samples : TYPE
        DESCRIPTION.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    None.

    """
    progressbar = tqdm.tqdm(total=samples.sequence_length, desc="Read VCF", unit='bp')
    pos = 0
    for variant in vcf:
        progressbar.update(variant.POS - pos)
        if pos == variant.POS:
            raise ValueError("Duplicate positions for variant at position", pos)
        else:
            pos = variant.POS
        if any([not phased for _, _, phased in variant.genotypes]):
            raise ValueError("Unphased genotypes for variant at position", pos)

        alleles = [variant.REF] + variant.ALT
        ancestral = variant.INFO.get('AA', variant.REF)
        ordered_alleles = [ancestral] + list(set(alleles) - {ancestral})
        allele_index = {old_index: ordered_alleles.index(allele) for old_index, allele in enumerate(alleles)}
        genotypes = [allele_index[old_index] for row in variant.genotypes for old_index in row[0:2]]

        samples.add_site(pos, genotypes=genotypes, alleles=alleles)
    progressbar.close()


def chrom_len(vcf):
    """Get chromosome length."""
    assert len(vcf.seqlens) == 1
    return vcf.seqlens[0]


def parse_args(args_in):
    """Parse args."""
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--vcf", required=True, type=str)
    parser.add_argument("--outfile", required=True, type=str)
    parser.add_argument("--meta", required=True, type=str,
                        help="metadata for names and populations. Columns must "
                        "include individualID and be tab delimited")
    parser.add_argument('-t', "--threads", type=int, default=1)
    parser.add_argument("--pops_header", type=str, default="Karyotype")

    return(parser.parse_args(args_in))


def main():
    """Run main function."""
    args = parse_args(sys.argv[1:])
    # =========================================================================
    #  Gather args
    # =========================================================================
    vcf_path = args.vcf
    outfile = args.outfile
    threads = args.threads
    label_by = args.pops_header
    meta = pd.read_csv(args.meta, sep="\t", index_col="sampleID", dtype=object)
    # =========================================================================
    #  Main executions
    # =========================================================================

    vcf = cyvcf2.VCF(vcf_path)
    with tsinfer.SampleData(path=f"{outfile}.samples",
                            sequence_length=chrom_len(vcf),
                            num_flush_threads=threads,
                            max_file_size=2**37) as samples:

        add_metadata(vcf, samples, meta, label_by)
        add_diploid_sites(vcf, samples)

    print(f"Sample file created for {samples.num_samples} samples ({samples.num_individuals}) with {samples.num_sites} variable sites.", flush=True)


if __name__ == "__main__":
    main()
