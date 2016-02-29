#!/usr/bin/env python

"""
# =============================================================================

Copyright Government of Canada 2015-2016

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

from neptune.Signature import *
from neptune.Utility import *

import unittest

"""
# =============================================================================

READ SIGNATURES

# =============================================================================
"""
class TestReadSignatures(unittest.TestCase):

    """ 
    # =============================================================================

    test_simple

    PURPOSE:
        Tests a simple read signatures example.

    INPUT:

        >long1 score=0.0000 in=0.0000 ex=0.0000 len=84 ref=reference1 pos=0
        ACTGAACCTTGGAAACCCTTTGGGAAAACCCCTTTTGGGGAAAAACCCCCTTTTTGGGGGAAAAAACCCCCCTTTTTTGGGGGG
        >long2 score=0.0000 in=0.0000 ex=0.0000 len=84 ref=reference3 pos=100
        ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT

    EXPECTED:

        signature["long1"] -> LONG1 SIGNATURE
        signature["long2"] -> LONG2 SIGNATURE

    # =============================================================================
    """
    def test_simple(self):

        fileLocation = getPath("tests/data/signature/multiple.fasta")
        signatures = readSignatures(fileLocation)

        self.assertEquals(signatures["long1"].ID, "long1")
        self.assertEquals(signatures["long1"].reference, "reference1")

        self.assertEquals(signatures["long2"].ID, "long2")
        self.assertEquals(signatures["long2"].reference, "reference3")

"""
# =============================================================================

WRITE SIGNATURE

# =============================================================================
"""
class TestWriteSignatures(unittest.TestCase):

    """ 
    # =============================================================================

    test_simple

    PURPOSE:
        Tests a simple write signatures example.

    INPUT:

        signature = Signature("0", 0.0, 0.0, 0.0, "ACGTACGT", "ref", "20")

    EXPECTED:

        >0 score=0.0000 in=0.0000 ex=0.0000 len=8 ref=ref pos=20
        ACGTACGT

    # =============================================================================
    """
    def test_simple(self):

        destination = StringIO.StringIO()
        signature = Signature("0", 0.0, 0.0, 0.0, "ACGTACGT", "ref", "20")

        writeSignature(signature, destination)
        result = destination.getvalue()

        expected = ">0 score=0.0000 in=0.0000 ex=0.0000 len=8 ref=ref pos=20\nACGTACGT\n"
        self.assertEquals(result, expected)

"""
# =============================================================================

SORT SIGNATURES

# =============================================================================
"""
class TestSortSignatures(unittest.TestCase):

    """ 
    # =============================================================================

    test_simple

    PURPOSE:
        Tests a simple signature sorting example.

    INPUT:

        signature1 = Signature("0", 0.1, 0.1, 0.0, "ACGTACGT", "ref1", "20")
        signature2 = Signature("0", 0.9, 0.9, 0.0, "ACGTACGT", "ref2", "20")
        signature3 = Signature("0", 0.5, 0.5, 0.0, "ACGTACGT", "ref3", "20")

    EXPECTED:

        [signature2, signature3, signature1]

    # =============================================================================
    """
    def test_simple(self):

        signature1 = Signature("1", 0.1, 0.1, 0.0, "ACGTACGT", "ref1", "20")
        signature2 = Signature("2", 0.9, 0.9, 0.0, "ACGTACGT", "ref2", "20")
        signature3 = Signature("3", 0.5, 0.5, 0.0, "ACGTACGT", "ref3", "20")

        signatures = {}
        signatures[signature1.ID] = signature1
        signatures[signature2.ID] = signature2
        signatures[signature3.ID] = signature3

        sortedSignatures = sortSignatures(signatures)

        self.assertEquals(sortedSignatures[0], signature2)
        self.assertEquals(sortedSignatures[1], signature3)
        self.assertEquals(sortedSignatures[2], signature1)

if __name__ == '__main__':
    
    unittest.main()
