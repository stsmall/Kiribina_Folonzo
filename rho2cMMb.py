#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 13:52:26 2019
@author: Scott T. Small

This script converts estimates of Rho from common programs into cMMb. Using a
standard conversion of Ne or the methods of Booker et al. 2017 for wild mice
https://doi.org/10.1534/genetics.117.300063

Correspondence from Tom Booker
-------
LD will potentially be influenced by a range of non-neutral processes
(e.g. selective sweeps, background selection, past demographic processes),

I did it this way was because
    a) 4Ner estimates are very noisy
    b) the distances between SNPs varies.
If one were to take the arithmetic mean of the rho estimates obtained from LD,
a physically lose pair of SNPs that, for one reason or another, have yielded a
spuriously high rho estimate could pull you mean rho for the chromosome up.
Accounting for the distances between SNPs by dealing with the cumulative rho
(what I called the frequency weighted mean in the Genetics paper).

Assumes:
    independent estimate of the total map length (in cM) for your chromosome
    else the code will express the cumulative rho as a proportion of the total
    rho.

# =============================================================================
# totalMapLength = 100 # e.g. mouse chromosome 1 has a map length of ~100cM.
# SNPpos = [1, 100, 1000]
# rho = [0.001, 0.005, 0.03]
# rhocum = []  # store the cumulative rho for each SNP position
#
# for i in range(len(SNPpos))
#      rhoTemp = (rho[i] * (SNPpos[i] - SNPpos[i-1]))
#      rhocum.append(rhocum[-1] + rhoTemp)
# cM = []
# for  j in rhocum:
#     cMperSNP =  (j / rhocum[-1])
#     cM.append(cMperSNP)
#
# =============================================================================

For each pair of adjacent SNPs, between which you have an estimate of 4Ner
1. calculate the recombination distance
    rhoTemp = (rho[i] * (SNPpos[i] - SNPpos[i-1]))
    # This is the contribution that the recombination distance
    # between a pair of SNPs makes to the total map length
2. add this to the total rho
    rhocum[-1] + rhoTemp
    # Take the cumulative sum by adding rhoTemp to the total
3. Keep a record of how much each successive SNP is adding to the total
    rhocum.append()
4. divide each position by the total rhocum and multiply that value by the
    total map length of the chromosome.
    cMperSNP = (j / rhocum[-1]) * totalMapLength
5. The list 'cM' will then contain estimates of the genetic positions of each
    SNP (excluding the first one) on your chromosome
    * If you didn't have an estimate of cM, ignore the multiplication.
    then the map as a proportion of the total rho for each chromosome.

Notes
-----
Map sizes for Afunestus from Wondji et al. 2005
https://doi.org/10.1534/genetics.105.044800
    X mapsize = 44.7
    2RL mapsize = 158
    3RL mapsize = 180
    3R mapsize = 90.9
    3L mapsize = 89.1
    2L mapsize = 63.2
    2R mapsize = 94.8


"""
import argparse
import sys


def read_helmet(helm_file):
    """Parses LDHelmet output and parses SNP position and Rho estimates

    Parameters
    ----------
    helm_file: str
        file containing data from LDHelmet

    Returns
    -------
    snp_list: List[int]
        list of snp positions
    rho_list: List[float]
        list of rho estimates between snps

    """
    snp_list = []
    rho_list = []
    with open(helm_file) as rhomap:
        for nline in rhomap:
            if nline.split()[0].isdigit():
                for line in rhomap:
                    try:
                        x = line.split()
                        snp = x[0]
                        rho = x[2]
                        snp_list.append(int(snp))
                        rho_list.append(float(rho))
                    except ValueError:
                        continue
    return snp_list, rho_list


def read_jump(jump_file):
    """Parses LDJump output and parses SNP position and Rho estimates

    Parameters
    ----------
    jump_file: str
        file containing data from LDJump

    Returns
    -------
    snp_list: List[int]
        list of snp positions
    rho_list: List[float]
        list of rho estimates between snps

    """
    snp_list = []
    rho_list = []
    with open(jump_file) as rhomap:
        for line in rhomap:
            try:
                x = line.split(",")
                snp = x[0]
                rho = x[1]
                snp_list.append(int(snp))
                rho_list.append(float(rho))
            except ValueError:
                continue
    return (snp_list, rho_list)


def read_ismc(ismc_file):
    """Parses a bedgraph output and parses SNP position and Rho estimates

    Parameters
    ----------
    ismc_file: str
        file containing data from LDJump

    Returns
    -------
    snp_list: List[int]
        list of snp positions
    rho_list: List[float]
        list of rho estimates between snps

    """
    snp_list = []
    rho_list = []
    with open(ismc_file) as rhomap:
        for line in rhomap:
            try:
                x = line.split()
                start = int(x[1])
                rho = float(x[3])
                snp_list.append(start)
                rho_list.append(rho)
            except ValueError:
                continue
    return snp_list, rho_list


def write_recomb(pos_list, cMMb_list, cM_list, prho_list):
    """Write output to file

    Parameters
    ----------
    pos_list: List[int]
        positions for recomb breakpoint
    cMMb_list: List[float]
        in cM per Mb
    cM_list: List[float]
        in cM for easy conversion

    Returns
    -------
    file

    """
    with open("cMMb.LD.out", 'w') as cmb:
        for i, j in zip(pos_list, cMMb_list):
            cmb.write(f"{i} {j}\n")
    with open("cM.LD.out", 'w') as cm:
        for i, j in zip(pos_list, cM_list):
            cm.write(f"{i} {j}\n")
    with open("proh.LD.out", 'w') as pr:
        for i, j in zip(pos_list, prho_list):
            pr.write(f"{i} {j}\n")


def recomb_map(snp_list, rho_list, Ne, map_size):
    """Standard recombination map conversion

    Parameters
    ----------
    snp_list: List[int]
        list of snp positions
    rho_list: List[float]
        list of rho estimates between snps
    Ne: int
        effective population size
    map_size: float
        map size of chromosome

    Returns
    -------
    pos_list: List[int]
        positions for recomb breakpoint
    cMMb_list: List[float]
        in cM per Mb
    cM_list: List[float]
        in cM for easy conversion

    """
    pos_list = []
    cMMb_list = []
    cM_list = []
    cM = 0
    # cumRho = 0
    for i, pos in enumerate(snp_list):
        cMMb_avg = ((rho_list[i]*100)/(4*Ne)) * 1E6
        # cMMb_avg = 50*log(1/(1-2*(rholist[i]/4*Ne))) * 1E6
        pos_list.append(pos)
        cMMb_rho = rho_list[i] * cMMb_avg  # average rate from Chan 2012
        cMMb_list.append(cMMb_rho)  # rate between SNPs
        if (i == 0):
            cM += (cMMb_rho * (pos)) / map_size
            # cumRho += (pos * rholist[i])/size
        else:
            cM += (cMMb_rho * (pos - snp_list[i-1])) / map_size
        cM_list.append(cM)
    return pos_list, cMMb_list, cM_list


def recomb_map_booker(snp_list, rho_snp_list, map_size):
    """Following the conversion method of Booker et al 2017 in genetics

    Parameters
    ----------
    snp_list: List[int]
        list of snp positions
    rho_list: List[float]
        list of rho estimates between snps
    map_size: float
        map size of chromosome

    Returns
    -------
    pos_list: List[int]
        positions for recomb breakpoint
    cMMb_list: List[float]
        in cM per Mb
    cM_list: List[float]
        in cM for easy conversion

    """
    pos_list = []
    rho_list = []
    total_rho = 0
    prho_list = []
    cM_list = []
    cMMb_list = []
    for i, rho in enumerate(rho_snp_list):
        if i == 0:
            rho_temp = (rho * snp_list[0])
            rho_list.append(rho_temp)
        else:
            rho_temp = rho * (snp_list[i] - snp_list[i-1])
            rho_list.append(total_rho + rho_temp)
        pos_list.append(snp_list[i])
        total_rho += rho_temp
    # rho to cM
    for i, rho in enumerate(rho_list):
        prho_snp = (rho / total_rho)
        prho_list.append(prho_snp)
        cM_list.append(prho_snp * map_size)
        if i == 0:
            cMMb_list.append(0)
        else:
            cMMb_list.append(((prho_snp - prho_list[i-1]) * map_size) /
                             ((snp_list[i] - snp_list[i-1])/1E6))
    return pos_list, cMMb_list, cM_list, prho_list


def parse_args(args_in):
    parser = argparse.ArgumentParser(prog="sys.argv[0].py",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--ld_file", required=True, help="recombination file")
    parser.add_argument("--map_size", type=float, default=1,
                        help="chrom map size else will use prop of total rho")
    parser.add_argument("--format", type=str, required=True,
                        choices=("LDhelmet", "LDJump", "iSMC", "ReLERNN"),
                        help="format of LD file")
    parser.add_argument("--EffectivePopSize", type=int)
    parser.add_argument("--booker", action="store_true",
                        help="use method in Booker et al. 2017")
    return (parser.parse_args(args_in))


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    # =========================================================================
    #  Gather args
    # =========================================================================
    MAP_SIZE = args.map_size
    FORMAT = args.format
    LD_FILE = args.ld_file
    BOOKER = args.booker
    NE = args.EffectivePopSize
    # =========================================================================
    #  Main executions
    # =========================================================================
    if FORMAT == "LDhelmet":
        SNP_LIST, RHO_LIST = read_helmet(LD_FILE)
    elif FORMAT == "LDJump":
        SNP_LIST, RHO_LIST = read_jump(LD_FILE)
    elif FORMAT == "iSMC":
        SNP_LIST, RHO_LIST = read_ismc(LD_FILE)
    elif FORMAT == "ReLERNN":
        pass
        # SNP_LIST, RHO_LIST = read_relernn(LD_FILE)
    if BOOKER is True:
        POS, CMMB, CM, PRHO = recomb_map_booker(SNP_LIST, RHO_LIST, MAP_SIZE)
    else:
        POS, CMMB, CM, PRHO = recomb_map(SNP_LIST, RHO_LIST, NE, MAP_SIZE)
    write_recomb(POS, CMMB, CM, PRHO)
