#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 12:03:48 2020
@author: Scott T. Small

This script tries to find a way to maximize both the number of sites and retained
individual

Example
-------

    $ python maximaize_retained_indvs.py

Notes
-----
    requires four files:
            KirFol.AfunF3.{chroms}.filter.mask.txt from apply_filter_to_vcf.py
            KirFol.AfunF3.{chroms}.uncalled.mask.txt from get_uncalled.py
            KirFol.AfunF3.{chroms}.miss.indv.txt

    1) remove highly missing individuals
    2) check site missingness : KirFol.AfunF3.{chroms}.miss.site.txt
    3) remove individuals that maximize keeping sites, only 10% missing/site


"""

# remove indv using miss.indv.txt
# calculate miss.site
# stepwise remove individuals from High-low to minimize the average site_missingness
# phasing will work w/ 20% missing data, instead of removing sites, remove indvs


import sys
import argparse


def parse_args(args_in):
    """Parse args."""
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument()
    return(parser.parse_args(args_in))


def main():
    """Run main function."""
    args = parse_args(sys.argv[1:])
    # =========================================================================
    #  Gather args
    # =========================================================================

    # =========================================================================
    #  Main executions
    # =========================================================================


if __name__ == "__main__":
    main()
