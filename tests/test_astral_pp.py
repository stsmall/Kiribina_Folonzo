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

from ..astral_pp import run_astral
from ..astral_pp import make_windows
import filecmp


def test_run_astral():
    """Test of getCDS
    """
    tree_file = "test_tree_file.txt"
    run_astral(tree_file, 5, "~/programs_that_work/ASTRAL/Astral/astral.5.6.1.jar", False)
    assert(filecmp.cmp("test_astral.tre", "astral.tre") == True)


def test_make_windows():
    """Compare files
    """
    coord_list = "test_coords.txt"
    make_windows(coord_list, 5, "3R")
    assert(filecmp.cmp("test_3R.windows.out", "3R.windows.out") == True)


if __name__ == "__main__":
    test_run_astral()
    test_make_windows()
