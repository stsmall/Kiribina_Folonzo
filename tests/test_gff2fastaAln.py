#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 20 15:51:27 2019
@author: Scott T. Small

?can github run tests and check coverage upon push?

unit testing for functions

# main unit testing
conda install -c anaconda pytest

# allows use in spyder
conda install -c spyder-ide spyder-unittest

# tests code against versions of python and modules
conda install -c conda-forge tox

# tests assertions beyond listed
conda install -c conda-forge hypothesis

# how much of application is tested by unit tests
conda install -c anaconda coverage

# mock values useful for multi-level testing where func is dependent on another
# MagicMock
conda install mock

# decorator with pytest to avoid loading heavy examples
@pytest.fixture(scope='module')
pytest.raises(RuntimeError)
pytest.warns(RuntimeWarning)
    warnings.warn("")

**Note**

"""

from gff2fastaAln import get_cds
from gff2fastaAln import get_ncds
from gff2fastaAln import format_fasta
from dataclasses import dataclass
import filecmp


@dataclass
class CDS:
    __slots__ = ["start", "end"]
    start: int
    end: int


def test_getCDS():
    """Test of getCDS
    """
    gffFile = "test_gff2fastaAln.gff"
    cds, chrom = get_cds(gffFile, 100, key="CDS")
    tcds = {
         "cds_0": CDS(80183, 81392),
         "cds_1": CDS(88350, 91331),
         "cds_2": CDS(112086, 112923),
         "cds_3": CDS(156134, 157030),
         "cds_4": CDS(174681, 180181),
         "cds_5": CDS(297373, 302451)
         }
    for k in cds.keys():
        assert(tcds[k].start == cds[k].start)
        assert(tcds[k].end == cds[k].end)
    assert(chrom == "3R")


def test_getNonCDS():
    """Test of getNonCDS
    """
    tcds = {
     "cds_0": CDS(50, 500),
     "cds_1": CDS(1000, 1200),
     "cds_2": CDS(10800, 11800),
     "cds_3": CDS(15000, 16000)
     }

    tncds = {
        "ncds_0": CDS(500, 1000),
        "ncds_1": CDS(1200, 2200),
        "ncds_2": CDS(4200, 5200),
        "ncds_3": CDS(7200, 8200),
        "ncds_4": CDS(10200, 10800),
        "ncds_5": CDS(11800, 12800),
        "ncds_6": CDS(14800, 15000),
        "ncds_7": CDS(16000, 17000),
        "ncds_8": CDS(19000, 20000)
        }
    ncds = get_ncds(tcds, 1000, 100, 2000, 20000)
    for k in ncds.keys():
        assert(tncds[k].start == ncds[k].start)
        assert(tncds[k].end == ncds[k].end)


def test_formatFasta():
    """Compare files
    """
    tcds = {
     "cds_0": CDS(5, 10),
     "cds_1": CDS(15, 25),
     "cds_2": CDS(40, 42),
     "cds_3": CDS(43, 50)
     }

    format_fasta("cds", tcds, "test_fastaFile.fa", 5, "3R", 0.5, False, just=10)
    assert(filecmp.cmp("test_CDS.bpp.3R.5-50.txt", "CDS.bpp.3R.5-50.txt") == True)


if __name__ == "__main__":
    test_getCDS()
    test_getNonCDS()
    test_formatFasta()