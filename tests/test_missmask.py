#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 15:51:27 2019
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
import filecmp
from missmask import miss_mask as miss_mask


def test_mask():
    """Test of masking function for missing data
    """
    vcfFile = "test_missmask.vcf.gz"
    mask_dict, chrom = miss_mask(vcfFile)
    assert(chrom == "Contig1")
    x = {"Kir1":["2", "6", "7", "8", "14"],"Fol1":["1", "2", "3", "6", "7", "8", "14", "15"]}
    assert(x == mask_dict)
    # bed file compare2files
    assert(filecmp.cmp("Genome.mask.test.bed", "Genome.mask.bed") == True)

if __name__ == "__main__":
    test_mask()
