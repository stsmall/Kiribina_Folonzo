# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 13:42:55 2021
@author: Scott T. Small

make input files for 3PCLR

Notes
-----

get derived counts for pop1:
    vcftools --vcf species.vcf --keep pop1.list --counts --derived --out pop1

get derived counts for pop2:
    vcftools --vcf species.vcf --keep pop2.list --counts --derived --out pop2

get derived counts for outgroup:
    vcftools --vcf species.vcf --keep outgroup.list --counts --derived --out outgroup

get freq of snps in outgroup:
    vcftools --vcf species.vcf --keep outgroup.list --freq2 --out outgroup

get positions that are polymorphic in outrgoup:
    awk '$5 < 1' outgroup.frq | cut -f2 > outgroup.polymorphic.txt

Keep these position for pop1 and pop2:
    fgrep -wf outgroup.polymorphic.txt pop1.count > pop1.polymorphic.count
    fgrep -wf outgroup.polymorphic.txt pop2.count > pop2.polymorphic.count
    fgrep -wf outgroup.polymorphic.txt outgroup.count > outgroup.polymorphic.count

python format_3PCLR.py --pop1 pop1.polymorphic.count --pop2 pop2.polymorphic.count
    --out outgroup.polymorphic.count --cM chr.cM.map.txt

"""
import argparse
import sys
import numpy as np
from collections import Counter
import tqdm


def parse_counts(file):
    dd = {}
    size = []
    with open(file) as f:
        for line in f:
            x = line.split()
            dd[x[1]] = x[3:]
            size.append(int(x[3]))
    return max(size), dd


def parse_cm(file):
    pos = []
    cm = []
    with open(file) as f:
        for line in f:
            c, p, m = line.split()
            pos.append(p)
            cm.append(m)

    pos_arr = np.array(pos, dtype=np.int64)
    cm_arr = np.array(cm, dtype=np.float64)
    return pos_arr, cm_arr


def input_3PCLR(pop1, pop2, out, cm):
    # load cm for interpolate
    pp, gp = parse_cm(cm)
    # load counts
    np1_n, pop1_dt = parse_counts(pop1)
    np2_n, pop2_dt = parse_counts(pop2)
    nout_n, out_dt = parse_counts(out)
    # set file and header
    names = pop1.split(".")
    chrom = names[0]
    pop1_n = names[1]
    pop2_n = pop2.split(".")[1]
    out_n = out.split(".")[1]
    f = open(f"{chrom}.3pclr.input", 'w')
    f.write(f"chr\tphyspos\tgenpos\tm{pop1_n}\tn{pop1_n}\tm{pop2_n}\tn{pop2_n}\tm{out_n}\tn{out_n}\n")
    # write outfile
    drift_ls = []
    for pos in tqdm(out_dt.keys()):
        cm_inp = np.interp(pos, pp, gp)
        anc = out_dt[pos][1].split(":")[0]
        mout = out_dt[pos][2].split(":")[1]
        nout = out_dt[pos][0]
        mp1 = "0"
        np1 = str(np1_n)
        mp2 = "0"
        np2 = str(np2_n)

        if pos in pop1_dt:
            mp1 = pop1_dt[pos][2].split(":")[1]
            np1 = pop1_dt[pos][0]
            anc = pop1_dt[pos][1].split(":")[0]
        if pos in pop2_dt:
            mp2 = pop2_dt[pos][2].split(":")[1]
            np2 = pop2_dt[pos][0]
            anc = pop2_dt[pos][1].split(":")[0]

        # get outgroup polarize
        if anc != out_dt[pos][1].split(":")[0]:
            mout = out_dt[pos][1].split(":")[1]
        # get concat line
        counts = "\t".join([mp1, np1, mp2, np2, mout, nout])
        drift_ls.append(counts)
        f.write(f"{chrom}\t{pos}\t{cm_inp}\t{counts}\n")
    f.close()
    return drift_ls


def input_drift(drift_ls):
    f = open("drift.input", 'w')
    f.write("DerivedAllelesPopA\tSampledAllelesPopA\tDerivedAllelesPopB\tSampledAllelesPopB\tDerivedAllelesPopC\tSampledAllelesPopC\tnum\n")
    allele_counts = Counter(drift_ls)
    for counts in allele_counts.keys():
        num = allele_counts[counts]
        f.write(f"{counts}\t{num}\n")
    f.close()


def parse_args(args_in):
    """Parse args."""
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser = argparse.ArgumentParser()
    parser.add_argument("--pop1", required=True,
                        help="derived counts pop 1")
    parser.add_argument("--pop2", required=True,
                        help="derived counts pop 2")
    parser.add_argument('--out', required=True,
                        help="derived counts outgroup")
    parser.add_argument('--cM', required=True,
                        help="genomic position in cM")
    return(parser.parse_args(args_in))


def main():
    """Run main function."""
    args = parse_args(sys.argv[1:])
    # =========================================================================
    #  Gather args
    # =========================================================================
    pop1 = args.pop1
    pop2 = args.pop2
    outpop = args.out
    cm = args.cM
    # =========================================================================
    #  Main executions
    # =========================================================================
    drift = input_3PCLR(pop1, pop2, outpop, cm)
    input_drift(drift)


if __name__ == "__main__":
    main()
