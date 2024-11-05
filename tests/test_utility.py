#!/usr/bin/env python

"""
# =============================================================================

Copyright Government of Canada 2015-2024

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
import io

from tests.TestingUtility import *

from neptune.Utility import *

import unittest

""" 
# =============================================================================

GET AGGREGATION TAGS

# =============================================================================
"""
class TestGetAggregationTags(unittest.TestCase):

    """ 
    # =============================================================================

    test_zero

    PURPOSE:
        Test a degree of parallelization of zero.

    INPUT:
        0: parallelization = 0

    EXPECTED:
        0: tags = ["", AGGREGATE_OTHER]

    # =============================================================================
    """
    def test_zero(self):

        # 0: parallelization = 0
        expected = ["", AGGREGATE_OTHER]
        result = getAggregationTags(int(0))
        self.assertSequenceEqual(result, expected)

    """ 
    # =============================================================================

    test_one

    PURPOSE:
        Test a degree of parallelization of one.

    INPUT:
        0: parallelization = 1

    EXPECTED:
        0: tags = ["A", "C", "G", "T", AGGREGATE_OTHER]

    # =============================================================================
    """
    def test_one(self):

        # 0: parallelization = 1
        expected = ["A", "C", "G", "T", AGGREGATE_OTHER]
        result = getAggregationTags(int(1))
        self.assertSequenceEqual(result, expected)

    """ 
    # =============================================================================

    test_two

    PURPOSE:
        Test a degree of parallelization of two.

    INPUT:
        0: parallelization = 2

    EXPECTED:
        0: tags = [
            "AA", "AC", "AG", "AT",
            "CA", "CC", "CG", "CT",
            "GA", "GC", "GG", "GT",
            "TA", "TC", "TG", "TT",
            AGGREGATE_OTHER]

    # =============================================================================
    """
    def test_two(self):

        # 0: parallelization = 2
        expected = [
            "AA", "AC", "AG", "AT",
            "CA", "CC", "CG", "CT",
            "GA", "GC", "GG", "GT",
            "TA", "TC", "TG", "TT",
            AGGREGATE_OTHER]
        result = getAggregationTags(int(2))
        self.assertSequenceEqual(result, expected)

    """ 
    # =============================================================================

    test_negative

    PURPOSE:
        Test a negative degree of parallelization.

    INPUT:
        0: parallelization = -1

    EXPECTED:
        0: RuntimeError

    # =============================================================================
    """
    def test_negative(self):

        # 0: parallelization = -1
        with self.assertRaises(RuntimeError):
            result = getAggregationTags(int(-1))

""" 
# =============================================================================

GENERATE SEQUENCE

# =============================================================================
"""
class TestGenerateSequence(unittest.TestCase):

    """ 
    # =============================================================================

    test_trivial

    PURPOSE:
        Test some trivial use cases.

    INPUT:
        0:
        integer = 0
        length = 1

        1:
        integer = 1
        length = 1

        2:
        integer = 2
        length = 1

        3:
        integer = 3
        length = 1

    EXPECTED:
        0: sequence = "A"
        1: sequence = "C"
        2: sequence = "G"
        3: sequence = "T"

    # =============================================================================
    """
    def test_trivial(self):

        # 0: integer = 0, length = 1
        result = generateSequence(0, 1)
        expected = "A"
        self.assertEqual(result, expected)

        # 1: integer = 1, length = 1
        result = generateSequence(1, 1)
        expected = "C"
        self.assertEqual(result, expected)

        # 2: integer = 2, length = 1
        result = generateSequence(2, 1)
        expected = "G"
        self.assertEqual(result, expected)

        # 3: integer = 3, length = 1
        result = generateSequence(3, 1)
        expected = "T"
        self.assertEqual(result, expected)

    """ 
    # =============================================================================

    test_large

    PURPOSE:
        Test some larger use cases.

    INPUT:
        0:
        integer = 0
        length = 3

        1:
        integer = 32
        length = 3

        2:
        integer = 0
        length = 10

    EXPECTED:
        0: sequence = "AAA"
        1: sequence = "GAA"
        2: sequence = "AAAAAAAAAA"

    # =============================================================================
    """
    def test_large(self):

        # 0: integer = 0, length = 3
        result = generateSequence(0, 3)
        expected = "AAA"
        self.assertEqual(result, expected)   

        # 1: integer = 32, length = 3
        result = generateSequence(32, 3)
        expected = "GAA"
        self.assertEqual(result, expected)

        # 2: integer = 0, length = 10
        result = generateSequence(0, 10)
        expected = "AAAAAAAAAA"
        self.assertEqual(result, expected)

    """ 
    # =============================================================================

    test_bounds

    PURPOSE:
        Test the bounds.

    INPUT:
        0:
        integer = 0
        length = 0

        1:
        integer = 0
        length = -1

        2:
        integer = -1
        length = 0

    EXPECTED:
        0: sequence = ""
        1: sequence = RuntimeError
        2: sequence = RuntimeError

    # =============================================================================
    """
    def test_bounds(self):

        # 0: integer = 0, length = 0
        result = generateSequence(0, 0)
        expected = ""
        self.assertEqual(result, expected)

        # 1: integer = 0, length = -1
        with self.assertRaises(RuntimeError):
            result = generateSequence(0, -1)

        # 2: integer = -1, length = 0
        with self.assertRaises(RuntimeError):
            result = generateSequence(-1, 0)

""" 
# =============================================================================

REVERSE COMPLEMENT

# =============================================================================
"""
class TestReverseComplement(unittest.TestCase):

    """ 
    # =============================================================================

    test_trivial

    PURPOSE:
        Test a trivial use case.

    INPUT:
        0: sequence = "A"

    EXPECTED:
        0: reverse = "T"

    # =============================================================================
    """
    def test_trivial(self):

        # 0: sequence = "A"
        result = reverseComplement("A")
        expected = "T"
        self.assertEqual(result, expected)

    """ 
    # =============================================================================

    test_large

    PURPOSE:
        Test a larger use case.

    INPUT:
        0: sequence = "ACGTACGTACGTA"

    EXPECTED:
        0: reverse = "TACGTACGTACGT"

    # =============================================================================
    """
    def test_large(self):

        # 0: sequence = "ACGTACGTACGTA"
        result = reverseComplement("ACGTACGTACGTA")
        expected = "TACGTACGTACGT"
        self.assertEqual(result, expected)  

    """ 
    # =============================================================================

    test_palindrome

    PURPOSE:
        Test the reverse complement of a palindrome.

    INPUT:
        0: sequence = "GATATAG"

    EXPECTED:
        0: reverse = "CTATATC"

    # =============================================================================
    """
    def test_palindrome(self):

        # 0: sequence = "GATATAG"
        result = reverseComplement("GATATAG")
        expected = "CTATATC"
        self.assertEqual(result, expected)

""" 
# =============================================================================

BUILD REFERENCES

# =============================================================================
"""
class TestBuildReferences(unittest.TestCase):

    """ 
    # =============================================================================

    test_simple

    PURPOSE:
        Test a simple use case.

    INPUT:
        0: referenceFile = 
            >0
            ACGTACGTACGT

    EXPECTED:
        0: references = {'0': "ACGTACGTACGT"}

    # =============================================================================
    """
    def test_simple(self):

        # 0:
        referenceFile = io.StringIO()
        referenceFile.write(">0\n")
        referenceFile.write("ACGTACGTACGT")

        referenceFile.seek(0)

        result = buildReferences(referenceFile)

        expected = {}
        expected["0"] = "ACGTACGTACGT"

        self.assertDictEqual(result, expected)

    """ 
    # =============================================================================

    test_multiline_fasta

    PURPOSE:
        Test a multiple-line fasta file.

    INPUT:
        0: referenceFile = 
            >0
            ACGTACGTACGT
            TTTTTTTTTTTT

    EXPECTED:
        0: references = {'0' : 'ACGTACGTACGTTTTTTTTTTTTT'}

    # =============================================================================
    """
    def test_multiline_fasta(self):

        # 0:
        referenceFile = io.StringIO()
        referenceFile.write(">0\n")
        referenceFile.write("ACGTACGTACGT\n")
        referenceFile.write("TTTTTTTTTTTT")

        referenceFile.seek(0)

        result = buildReferences(referenceFile)

        expected = {}
        expected["0"] = "ACGTACGTACGTTTTTTTTTTTTT"

        self.assertDictEqual(result, expected)

    """ 
    # =============================================================================

    test_multirecord_fasta

    PURPOSE:
        Test a multi-record fasta file.

    INPUT:
        0: referenceFile = 
            >0
            ACGTACGTACGT
            >1
            TTTTTTTTTTTT

    EXPECTED:
        0: references = {'0' : 'ACGTACGTACGT', '1' : 'TTTTTTTTTTTT'}

    # =============================================================================
    """
    def test_multirecord_fasta(self):

        # 0:
        referenceFile = io.StringIO()
        referenceFile.write(">0\n")
        referenceFile.write("ACGTACGTACGT\n")
        referenceFile.write(">1\n")
        referenceFile.write("TTTTTTTTTTTT\n")

        referenceFile.seek(0)

        result = buildReferences(referenceFile)

        expected = {}
        expected["0"] = "ACGTACGTACGT"
        expected["1"] = "TTTTTTTTTTTT"

        self.assertDictEqual(result, expected)

""" 
# =============================================================================

ESTIMATE REFERENCE PARAMETERS

# =============================================================================
"""
class TestEstimateReferenceParameters(unittest.TestCase):

    """ 
    # =============================================================================

    test_simple

    PURPOSE:
        Test a simple use case.

    INPUT:
        0: references = {'1': "ACGTACGTACGT"}

    EXPECTED:
        0: size = 12, gcContent = 0.5

    # =============================================================================
    """
    def test_simple(self):

        # 0: references = {'1': "ACGTACGTACGT"}
        references = {'1': "ACGTACGTACGT"}

        size, gcContent = estimateReferenceParameters(references)

        expectedSize = 12
        expectedGCContent = 0.5

        self.assertEqual(size, expectedSize)
        self.assertAlmostEqual(gcContent, expectedGCContent, 5)

    """ 
    # =============================================================================

    test_AT

    PURPOSE:
        Test all AT.

    INPUT:
        0: references = {'1': "ATATATATATAT"}

    EXPECTED:
        0: size = 12, gcContent = 0

    # =============================================================================
    """
    def test_AT(self):

        # 0: references = {'1': "ATATATATATAT"}
        references = {'1': "ATATATATATAT"}

        size, gcContent = estimateReferenceParameters(references)

        expectedSize = 12
        expectedGCContent = 0

        self.assertEqual(size, expectedSize)
        self.assertAlmostEqual(gcContent, expectedGCContent, 5)

    """ 
    # =============================================================================

    test_GC

    PURPOSE:
        Test all GC.

    INPUT:
        0: references = {'1': "GCGCGCGCGCGC"}

    EXPECTED:
        0: size = 12, gcContent = 1.0

    # =============================================================================
    """
    def test_GC(self):

        # 0: references = {'1': "GCGCGCGCGCGC"}
        references = {'1': "GCGCGCGCGCGC"}

        size, gcContent = estimateReferenceParameters(references)

        expectedSize = 12
        expectedGCContent = 1.0

        self.assertEqual(size, expectedSize)
        self.assertAlmostEqual(gcContent, expectedGCContent, 5)

    """ 
    # =============================================================================

    test_multiple_references

    PURPOSE:
        Test multiple references.

    INPUT:
        0: references = {'1': "ACGTACGTACGT", '2': "ACGT"}

    EXPECTED:
        0: size = 12, gcContent = 0.5

    # =============================================================================
    """
    def test_multiple_references(self):

        # 0: references = {'1': "ACGTACGTACGT", '2': "ACGT"}
        references = {'1': "ACGTACGTACGT", '2': "ACGT"}

        size, gcContent = estimateReferenceParameters(references)

        expectedSize = 16
        expectedGCContent = 0.5

        self.assertEqual(size, expectedSize)
        self.assertAlmostEqual(gcContent, expectedGCContent, 5)

    """ 
    # =============================================================================

    test_no_references

    PURPOSE:
        Test no references.

    INPUT:
        0: references = {}

    EXPECTED:
        0: RuntimeError

    # =============================================================================
    """
    def test_no_references(self):

        # 0: references = {}
        references = {}

        with self.assertRaises(RuntimeError):
            size, gcContent = estimateReferenceParameters(references)

if __name__ == '__main__':
    
    unittest.main()   
