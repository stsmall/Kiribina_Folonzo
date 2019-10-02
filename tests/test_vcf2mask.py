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
from vcf2mask import readmask as readmask

def test_mask():
    """Test of masking function
    """
    mask_file = "1sample.vcf.gz"
    sample = "KirFol1"
    chrom = "Contig0"
    test_list = ["500", "505", "11931", "12897", "13218", "17642"]
    sample_name, chrom_name, mask_list = readmask(mask_file, [".", "PASS"])
    assert(sample == sample_name)
    assert(chrom == chrom_name)
    assert(test_list == mask_list)

if __name__ == "__main__":
    test_mask()
