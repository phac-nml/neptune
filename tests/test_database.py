#!/usr/bin/env python

"""
# =============================================================================

Copyright Government of Canada 2015

Written by: Eric Marinier, Public Health Agency of Canada,
    National Microbiology Laboratory

Funded by the National Micriobiology Laboratory and the Genome Canada / Alberta
    Innovates Bio Solutions project "Listeria Detection and Surveillance
    using Next Generation Genomics"

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this file except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

# =============================================================================
"""

import os
import sys
import StringIO

from TestingUtility import *
prepareSystemPath()

from neptune.Database import *
from neptune.Utility import *

import unittest

"""
# =============================================================================

QUERY DATABASE

# =============================================================================
"""
class TestQueryDataBase(unittest.TestCase):

    """ 
    # =============================================================================

    test_simple

    PURPOSE:
        Tests a simple query.

    INPUT:
        0:

        database constructed from:
        ACTGAACCTTGGAAACCCTTTGGGAAAACCCCTTTTGGGG\
        AAAAACCCCCTTTTTGGGGGAAAAAACCCCCCTTTTTTGGGGGG

        query:
        AAACCCTTTGGGAAAACCCCTTTTGGGGAAAAA

    EXPECTED:
        0: "long.query\t33\tlong\t33\t100.00\t33\n" in result

    # =============================================================================
    """
    def test_simple(self):

        databaseLocation = "tests/data/filter/long.database/LONG"
        queryLocation = "tests/data/filter/long.query"
        outputLocation = "tests/output/filter/temp.out"
        filterPercent = 0.50
        seedSize = 11

        queryOutputLocation = queryDatabase(databaseLocation, queryLocation, outputLocation,
            filterPercent, seedSize)

        with open (queryOutputLocation, "r") as myfile:
            result = myfile.read()

        expected = "long.query\t33\tlong\t33\t100.00\t33\n"

        self.assertTrue(expected in result)

        os.remove(queryOutputLocation)


    """ 
    # =============================================================================

    test_missing

    PURPOSE:
        Tests a query that does not exist.

    INPUT:
        0:

        database constructed from:
        ACTGAACCTTGGAAACCCTTTGGGAAAACCCCTTTTGGGGAAAAACCCCCTTTTTGGGGGAAAAAACCCCCCTTTTTTGGGGGG

        missing.query:
        ATATATATATATATATATATATATATAT

    EXPECTED:
        0: ""

    # =============================================================================
    """
    def test_missing(self):

        databaseLocation = "tests/data/filter/long.database/LONG"
        queryLocation = "tests/data/filter/missing.query"
        outputLocation = "tests/output/filter/temp.out"
        filterPercent = 0.50
        seedSize = 11

        queryOutputLocation = queryDatabase(databaseLocation, queryLocation, outputLocation,
            filterPercent, seedSize)

        with open (queryOutputLocation, "r") as myfile:
            result = myfile.read()

        expected = ""

        self.assertEquals(result, expected)

        os.remove(queryOutputLocation)

if __name__ == '__main__':
    
    unittest.main()
