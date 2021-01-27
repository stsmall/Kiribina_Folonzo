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
print("GLOBAL declared for chrom lengths, change line 26 or comment to infer from VCF")
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
            gt = allele_index[gt[0]] + gt[1] + allele_index[gt[2]]
        except KeyError:
            # ./.
            pass
        genotypes[i] = gt
    return genotypes


def repolarize(vcf_file, mask, ancprob=0.90):
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

    # set progress bar and chrom lengths
    with fopen(vcf_file, 'rt') as pgb:
        for line in pgb:
            if line.startswith("##contig"):
                if CHROMDICT:
                    add_contig = True
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

    if mask:
        maskfile = open(f"chr{chrom}.polarize.mask.bed", 'w')
        outfile = open(f'{vcf_file.rstrip(".vcf.gz")}.derived.noerr.vcf', 'w')
    else:
        outfile = open(f'{vcf_file.rstrip(".vcf.gz")}.derived.vcf', 'w')
    # read vcf
    low_count = 0
    site_count = 0
    reorder_count = 0
    with fopen(vcf_file, 'rt') as vcf:
        for line in vcf:
            if line.startswith("##contig"):
                if add_contig:
                    for c, l in CHROMDICT.items():
                        outfile.write(f"##contig=<ID={c},length={l}>\n")
                    add_contig = False
            elif not line.startswith("#"):
                v_cols = line.split()
                chrom = v_cols[0]
                var_pos = int(v_cols[1])
                progressbar.update(var_pos - pos)
                pos = var_pos
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
                        if "AA=" in field:
                            ancestral = field.split("=")[1]
                        if "AAProb=" in field:
                            probability = float(field.split("=")[1])
                    if probability < ancprob:
                        low_count += 1
                        if mask is True:
                            maskfile.write(f"{chrom}\t{var_pos - 1}\t{var_pos}\n")
                            continue
                except IndexError:
                    print("Improper AA field")
                    sys.exit()
                # reset genotype ID
                if ref != ancestral:
                    reorder_count += 1
                    # rearrange so derived is ref
                    ordered_alleles = [ancestral] + list(set(alleles) - {ancestral})
                    allele_index = {str(old_index): str(ordered_alleles.index(allele))
                                    for old_index, allele in enumerate(alleles)}
                    # reset all genotype fields
                    genotypes = v_cols[9:]
                    genotypes = reset_genotypes(genotypes, allele_index)
                    # replace vcf fields
                    v_cols[3] = ancestral
                    v_cols[4] = ",".join(list(set(alleles)-{ancestral}))
                    v_cols[9:] = genotypes
                    v_cols_line = "\t".join(v_cols)
                    # write back to vcf line
                    outfile.write(f"{v_cols_line}\n")
                else:
                    outfile.write(line)
            else:
                outfile.write(line)
    outfile.close()
    if mask:
        maskfile.close()
    return low_count, site_count


def parse_args(args_in):
    """Parse args."""
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("vcf", type=str, help="VCF type file to reorder")
    parser.add_argument("--ancprob", type=float, help="cutoff for HQ")
    parser.add_argument("--mask", action="store_true")
    return(parser.parse_args(args_in))


def main():
    """Run main function."""
    args = parse_args(sys.argv[1:])
    # =========================================================================
    #  Gather args
    # =========================================================================
    vcf_file = args.vcf
    ancprob = args.ancprob
    mask = args.mask
    # =========================================================================
    #  Main executions
    # =========================================================================
    low, site = repolarize(vcf_file, mask, ancprob)
    print(f"low prob ancestral state = {low}, {low/site}")


if __name__ == "__main__":
    main()
