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

from neptune.FilterSignatures import *
from neptune.Utility import *
from neptune.Database import *

import unittest

class DefaultFilterSignatures():

    def __init__(self):

        candidatesLocation = getPath("tests/data/filter/long.fasta")
        filteredLocation = getPath("tests/output/filter/temp.out")
        sortedLocation = getPath("tests/output/filter/temp.out")
        totalInclusion = 1
        totalExclusion = 1
        filterLength = 0.5

        self.default = FilterSignatures(
            candidatesLocation, filteredLocation, sortedLocation,
            totalInclusion, totalExclusion, filterLength)

"""
# =============================================================================

UPDATE HIT OVERALL DICTIONARY

# =============================================================================
"""
class TestUpdateExclusionOverallDictionary(unittest.TestCase):

    """ 
    # =============================================================================

    test_simple

    PURPOSE:
        Tests a simple example where entries are continually added to the
        dictionary.

    INPUT:

        hit1 = "query1 10 subject 20 0.50 1"
        hit2 = "query2 10 subject 20 0.50 1"
        hit3 = "query1 13 subject 20 0.60 3"
        hit4 = "query2 13 subject 20 0.60 3"

        0: hit1 = "query1 10 subject 20 0.50 1"
        1: hit2 = "query2 10 subject 20 0.50 1"
        2: hit3 = "query1 13 subject 20 0.60 3"
        3: hit4 = "query2 13 subject 20 0.60 3"

    EXPECTED:

        0:

        bestOverall:
            query1 -> hit1

        1:

        bestOverall:
            query1 -> hit1
            query2 -> hit2

        2:

        bestOverall:
            query1 -> hit3
            query2 -> hit2

        3:

        bestOverall:
            query1 -> hit3
            query2 -> hit4

    # =============================================================================
    """
    def test_simple(self):

        filterSignatures = DefaultFilterSignatures().default

        query1 = "query1"
        query2 = "query2"

        # empty
        self.assertDictEqual(filterSignatures.exclusionOverallDictionary, {})

        #0:
        hit1 = Hit(str(query1) + " 10 subject 20 0.50 1")
        filterSignatures.updateExclusionOverallDictionary(hit1)
        expected = {}
        expected[query1] = hit1
        self.assertDictEqual(filterSignatures.exclusionOverallDictionary, expected)

        #1:
        hit2 = Hit(str(query2) + " 10 subject 20 0.50 1")
        filterSignatures.updateExclusionOverallDictionary(hit2)
        expected = {}
        expected[query1] = hit1
        expected[query2] = hit2
        self.assertDictEqual(filterSignatures.exclusionOverallDictionary, expected)

        #2:
        hit3 = Hit(str(query1) + " 13 subject 20 0.60 3")
        filterSignatures.updateExclusionOverallDictionary(hit3)
        expected = {}
        expected[query1] = hit3
        expected[query2] = hit2
        self.assertDictEqual(filterSignatures.exclusionOverallDictionary, expected)

        #3:
        hit4 = Hit(str(query2) + " 13 subject 20 0.60 3")
        filterSignatures.updateExclusionOverallDictionary(hit4)
        expected = {}
        expected[query1] = hit3
        expected[query2] = hit4
        self.assertDictEqual(filterSignatures.exclusionOverallDictionary, expected)

"""
# =============================================================================

UPDATE BEST PAIRS

# =============================================================================
"""
class TestUpdatePairDictionary(unittest.TestCase):

    """ 
    # =============================================================================

    test_simple

    PURPOSE:
        Tests a simple example where items are continually added to the dictionary.

    INPUT:

        0: hit1 = Hit("query1 10 subject1 20 0.50 1")
        1: hit2 = Hit("query1 10 subject2 20 0.50 1")
        2: hit3 = Hit("query2 10 subject1 20 0.50 1")
        3: hit4 = Hit("query2 10 subject2 20 0.50 1")
        4: hit5 = Hit("query1 13 subject1 20 0.60 3")
        5: hit6 = Hit("query1 13 subject2 20 0.60 3")
        6: hit7 = Hit("query2 13 subject1 20 0.60 3")
        7: hit8 = Hit("query2 13 subject2 20 0.60 3")


    EXPECTED:

        0: bestPairs = {hit1}
        1: bestPairs = {hit1, hit2}
        2: bestPairs = {hit1, hit2, hit3}
        3: bestPairs = {hit1, hit2, hit3, hit4}
        4: bestPairs = {hit5, hit2, hit3, hit4}
        5: bestPairs = {hit5, hit6, hit3, hit4}
        6: bestPairs = {hit5, hit6, hit7, hit4}
        7: bestPairs = {hit5, hit6, hit7, hit8}


    # =============================================================================
    """
    def test_simple(self):

        filterSignatures = DefaultFilterSignatures().default

        pair1 = ("query1", "subject1")
        pair2 = ("query1", "subject2")
        pair3 = ("query2", "subject1")
        pair4 = ("query2", "subject2")

        hit1 = Hit("query1 10 subject1 20 0.50 1")
        hit2 = Hit("query1 10 subject2 20 0.50 1")
        hit3 = Hit("query2 10 subject1 20 0.50 1")
        hit4 = Hit("query2 10 subject2 20 0.50 1")
        hit5 = Hit("query1 13 subject1 20 0.60 3")
        hit6 = Hit("query1 13 subject2 20 0.60 3")
        hit7 = Hit("query2 13 subject1 20 0.60 3")
        hit8 = Hit("query2 13 subject2 20 0.60 3")

        # empty
        self.assertDictEqual(filterSignatures.inclusionPairDictionary, {})

        #0:
        filterSignatures.updatePairDictionary(filterSignatures.inclusionPairDictionary, hit1)
        expected = {}
        expected[pair1] = hit1
        self.assertDictEqual(filterSignatures.inclusionPairDictionary, expected)

        #1:
        filterSignatures.updatePairDictionary(filterSignatures.inclusionPairDictionary, hit2)
        expected = {}
        expected[pair1] = hit1
        expected[pair2] = hit2
        self.assertDictEqual(filterSignatures.inclusionPairDictionary, expected)

        #2:
        filterSignatures.updatePairDictionary(filterSignatures.inclusionPairDictionary, hit3)
        expected = {}
        expected[pair1] = hit1
        expected[pair2] = hit2
        expected[pair3] = hit3
        self.assertDictEqual(filterSignatures.inclusionPairDictionary, expected)

        #3:
        filterSignatures.updatePairDictionary(filterSignatures.inclusionPairDictionary, hit4)
        expected = {}
        expected[pair1] = hit1
        expected[pair2] = hit2
        expected[pair3] = hit3
        expected[pair4] = hit4
        self.assertDictEqual(filterSignatures.inclusionPairDictionary, expected)

        #4:
        filterSignatures.updatePairDictionary(filterSignatures.inclusionPairDictionary, hit5)
        expected = {}
        expected[pair1] = hit5
        expected[pair2] = hit2
        expected[pair3] = hit3
        expected[pair4] = hit4
        self.assertDictEqual(filterSignatures.inclusionPairDictionary, expected)

        #5:
        filterSignatures.updatePairDictionary(filterSignatures.inclusionPairDictionary, hit6)
        expected = {}
        expected[pair1] = hit5
        expected[pair2] = hit6
        expected[pair3] = hit3
        expected[pair4] = hit4
        self.assertDictEqual(filterSignatures.inclusionPairDictionary, expected)

        #6:
        filterSignatures.updatePairDictionary(filterSignatures.inclusionPairDictionary, hit7)
        expected = {}
        expected[pair1] = hit5
        expected[pair2] = hit6
        expected[pair3] = hit7
        expected[pair4] = hit4
        self.assertDictEqual(filterSignatures.inclusionPairDictionary, expected)

        #7:
        filterSignatures.updatePairDictionary(filterSignatures.inclusionPairDictionary, hit8)
        expected = {}
        expected[pair1] = hit5
        expected[pair2] = hit6
        expected[pair3] = hit7
        expected[pair4] = hit8
        self.assertDictEqual(filterSignatures.inclusionPairDictionary, expected)

"""
# =============================================================================

UPDATE EXCLUSION SCORES

# =============================================================================
"""
class TestUpdateExclusionScores(unittest.TestCase):

    """ 
    # =============================================================================

    test_simple

    PURPOSE:
        Tests a simple example.

    INPUT:
        0:

        hit1 = Hit("query1 20 subject1 10 20 1")
        hit2 = Hit("query1 20 subject2 10 40 1")
        hit3 = Hit("query2 20 subject1 10 60 1")
        hit4 = Hit("query2 20 subject2 10 80 1")

    EXPECTED:
        0:

        query1 -> -0.15
        query2 -> -0.35

    # =============================================================================
    """
    def test_simple(self):

        filterSignatures = DefaultFilterSignatures().default
        filterSignatures.totalExclusion = 2

        hit1 = Hit("query1 20 subject1 10 20 1")
        hit2 = Hit("query1 20 subject2 10 40 1")
        hit3 = Hit("query2 20 subject1 10 60 1")
        hit4 = Hit("query2 20 subject2 10 80 1")

        filterSignatures.updatePairDictionary(filterSignatures.exclusionPairDictionary, hit1)
        filterSignatures.updatePairDictionary(filterSignatures.exclusionPairDictionary, hit2)
        filterSignatures.updatePairDictionary(filterSignatures.exclusionPairDictionary, hit3)
        filterSignatures.updatePairDictionary(filterSignatures.exclusionPairDictionary, hit4)

        # empty
        self.assertDictEqual(filterSignatures.overallScore, {})
        self.assertDictEqual(filterSignatures.exclusionScore, {})

        filterSignatures.updateExclusionScores()
  
        self.assertAlmostEqual(filterSignatures.overallScore["query1"], -0.15)
        self.assertAlmostEqual(filterSignatures.overallScore["query2"], -0.35)

        self.assertAlmostEqual(filterSignatures.exclusionScore["query1"], -0.15)
        self.assertAlmostEqual(filterSignatures.exclusionScore["query2"], -0.35)

    """ 
    # =============================================================================

    test_longer_alignment

    PURPOSE:
        Tests an example where the alignment is longer than the query length.

    INPUT:
        0:

        hit1 = Hit("query1 10 subject1 20 10 1")
        hit2 = Hit("query1 10 subject2 20 20 1")
        hit3 = Hit("query2 10 subject1 20 30 1")
        hit4 = Hit("query2 10 subject2 20 40 1")

    EXPECTED:
        0:

        query1 -> -0.3
        query2 -> -0.7

    # =============================================================================
    """
    def test_longer_alignment(self):

        filterSignatures = DefaultFilterSignatures().default
        filterSignatures.totalExclusion = 2

        hit1 = Hit("query1 10 subject1 20 10 1")
        hit2 = Hit("query1 10 subject2 20 20 1")
        hit3 = Hit("query2 10 subject1 20 30 1")
        hit4 = Hit("query2 10 subject2 20 40 1")

        filterSignatures.updatePairDictionary(filterSignatures.exclusionPairDictionary, hit1)
        filterSignatures.updatePairDictionary(filterSignatures.exclusionPairDictionary, hit2)
        filterSignatures.updatePairDictionary(filterSignatures.exclusionPairDictionary, hit3)
        filterSignatures.updatePairDictionary(filterSignatures.exclusionPairDictionary, hit4)

        # empty
        self.assertDictEqual(filterSignatures.overallScore, {})
        self.assertDictEqual(filterSignatures.exclusionScore, {})

        filterSignatures.updateExclusionScores()
  
        self.assertAlmostEqual(filterSignatures.overallScore["query1"], -0.30)
        self.assertAlmostEqual(filterSignatures.overallScore["query2"], -0.70)

        self.assertAlmostEqual(filterSignatures.exclusionScore["query1"], -0.30)
        self.assertAlmostEqual(filterSignatures.exclusionScore["query2"], -0.70)

"""
# =============================================================================

UPDATE INCLUSION SCORES

# =============================================================================
"""
class TestUpdateInclusionScores(unittest.TestCase):

    """ 
    # =============================================================================

    test_simple

    PURPOSE:
        Tests a simple example.

    INPUT:
        0:

        hit1 = Hit("query1 20 subject1 10 20 1")
        hit2 = Hit("query1 20 subject2 10 40 1")
        hit3 = Hit("query2 20 subject1 10 60 1")
        hit4 = Hit("query2 20 subject2 10 80 1")

    EXPECTED:
        0:

        query1 -> 0.15
        query2 -> 0.35

    # =============================================================================
    """
    def test_simple(self):

        filterSignatures = DefaultFilterSignatures().default
        filterSignatures.totalInclusion = 2

        hit1 = Hit("query1 20 subject1 10 20 1")
        hit2 = Hit("query1 20 subject2 10 40 1")
        hit3 = Hit("query2 20 subject1 10 60 1")
        hit4 = Hit("query2 20 subject2 10 80 1")

        filterSignatures.updatePairDictionary(filterSignatures.inclusionPairDictionary, hit1)
        filterSignatures.updatePairDictionary(filterSignatures.inclusionPairDictionary, hit2)
        filterSignatures.updatePairDictionary(filterSignatures.inclusionPairDictionary, hit3)
        filterSignatures.updatePairDictionary(filterSignatures.inclusionPairDictionary, hit4)

        # empty
        self.assertDictEqual(filterSignatures.overallScore, {})
        self.assertDictEqual(filterSignatures.inclusionScore, {})

        filterSignatures.updateInclusionScores()
  
        self.assertAlmostEqual(filterSignatures.overallScore["query1"], 0.15)
        self.assertAlmostEqual(filterSignatures.overallScore["query2"], 0.35)

        self.assertAlmostEqual(filterSignatures.inclusionScore["query1"], 0.15)
        self.assertAlmostEqual(filterSignatures.inclusionScore["query2"], 0.35)

    """ 
    # =============================================================================

    test_longer_alignment

    PURPOSE:
        Tests an example where the alignment is longer than the query length.

    INPUT:
        0:

        hit1 = Hit("query1 10 subject1 20 10 1")
        hit2 = Hit("query1 10 subject2 20 20 1")
        hit3 = Hit("query2 10 subject1 20 30 1")
        hit4 = Hit("query2 10 subject2 20 40 1")

    EXPECTED:
        0:

        query1 -> 0.3
        query2 -> 0.7

    # =============================================================================
    """
    def test_longer_alignment(self):

        filterSignatures = DefaultFilterSignatures().default
        filterSignatures.totalInclusion = 2

        hit1 = Hit("query1 10 subject1 20 10 1")
        hit2 = Hit("query1 10 subject2 20 20 1")
        hit3 = Hit("query2 10 subject1 20 30 1")
        hit4 = Hit("query2 10 subject2 20 40 1")

        filterSignatures.updatePairDictionary(filterSignatures.inclusionPairDictionary, hit1)
        filterSignatures.updatePairDictionary(filterSignatures.inclusionPairDictionary, hit2)
        filterSignatures.updatePairDictionary(filterSignatures.inclusionPairDictionary, hit3)
        filterSignatures.updatePairDictionary(filterSignatures.inclusionPairDictionary, hit4)

        # empty
        self.assertDictEqual(filterSignatures.overallScore, {})
        self.assertDictEqual(filterSignatures.inclusionScore, {})

        filterSignatures.updateInclusionScores()
  
        self.assertAlmostEqual(filterSignatures.overallScore["query1"], 0.30)
        self.assertAlmostEqual(filterSignatures.overallScore["query2"], 0.70)

        self.assertAlmostEqual(filterSignatures.inclusionScore["query1"], 0.30)
        self.assertAlmostEqual(filterSignatures.inclusionScore["query2"], 0.70)

"""
# =============================================================================

REPORT CANDIDATES

# =============================================================================
"""
class TestReportCandidates(unittest.TestCase):

    """ 
    # =============================================================================

    test_simple

    PURPOSE:
        Tests a simple example.

    INPUT:
        0:

    EXPECTED:
        >long score=0.0000 in=0.0000 ex=0.0000 len=84 ref=reference pos=0
        ACTGAACCTTGGAAACCCTTTGGGAAAACCCCTTTTGGGGAAAAACCCCCTTTTTGGGGGAAAAAACCCCCCTTTTTTGGGGGG

    # =============================================================================
    """
    def test_simple(self):

        filterSignatures = DefaultFilterSignatures().default
        
        filterSignatures.candidatesLocation = getPath("tests/data/filter/long.fasta")
        filterSignatures.filteredLocation = getPath("tests/output/filter/temp.out")
        filterSignatures.filterLength = 0.5

        hit1 = Hit("long 84 subject 40 100 10")
        filterSignatures.updateExclusionOverallDictionary(hit1)

        filterSignatures.reportFilteredCandidates()

        with open (filterSignatures.filteredLocation, "r") as myfile:

            result = myfile.read()
            expected = ">long score=0.0000 in=0.0000 ex=0.0000 len=84 ref=reference pos=0\nACTGAACCTTGGAAACCCTTTGGGAAAACCCCTTTTGGGGAAAAACCCCCTTTTTGGGGGAAAAAACCCCCCTTTTTTGGGGGG\n"
            self.assertEqual(result, expected)

        os.remove(filterSignatures.filteredLocation)

    """ 
    # =============================================================================

    test_too_large

    PURPOSE:
        Tests a candidate that's too large.

    INPUT:
        0:

    EXPECTED:
        0:

    # =============================================================================
    """
    def test_too_large(self):

        filterSignatures = DefaultFilterSignatures().default
        
        filterSignatures.candidatesLocation = getPath("tests/data/filter/long.fasta")
        filterSignatures.filteredLocation = getPath("tests/output/filter/temp.out")
        filterSignatures.filterLength = 0.5

        hit1 = Hit("long 84 subject 50 100 10")
        filterSignatures.updateExclusionOverallDictionary(hit1)

        filterSignatures.reportFilteredCandidates()

        with open (filterSignatures.filteredLocation, "r") as myfile:

            result = myfile.read()
            expected = ""
            self.assertEqual(result, expected)

        os.remove(filterSignatures.filteredLocation)

"""
# =============================================================================

REPORT SORTED

# =============================================================================
"""
class TestReportSorted(unittest.TestCase):

    """ 
    # =============================================================================

    test_simple

    PURPOSE:
        Tests a simple example.

    INPUT:
        0:

    EXPECTED:
        >long score=0.1500 in=0.1500 ex=0.0000 len=84 ref=reference pos=0
        ACTGAACCTTGGAAACCCTTTGGGAAAACCCCTTTTGGGGAAAAACCCCCTTTTTGGGGGAAAAAACCCCCCTTTTTTGGGGGG

    # =============================================================================
    """
    def test_simple(self):

        filterSignatures = DefaultFilterSignatures().default

        filterSignatures.filteredLocation = getPath("tests/data/filter/long.fasta")
        filterSignatures.sortedLocation = getPath("tests/output/filter/temp.out")
        filterSignatures.totalInclusion = 2

        hit1 = Hit("long 84 subject1 42 20 1")
        hit2 = Hit("long 84 subject2 42 40 1")

        filterSignatures.updatePairDictionary(filterSignatures.inclusionPairDictionary, hit1)
        filterSignatures.updatePairDictionary(filterSignatures.inclusionPairDictionary, hit2)

        # empty
        self.assertDictEqual(filterSignatures.overallScore, {})
        self.assertDictEqual(filterSignatures.inclusionScore, {})

        filterSignatures.updateInclusionScores()
  
        self.assertAlmostEqual(filterSignatures.overallScore["long"], 0.15)
        self.assertAlmostEqual(filterSignatures.inclusionScore["long"], 0.15)

        sortedSignatureIDs = [ID for (ID, score) in sorted(
            filterSignatures.overallScore.items(), key=operator.itemgetter(1), reverse=True)]

        filterSignatures.reportSorted(sortedSignatureIDs)

        with open (filterSignatures.sortedLocation, "r") as myfile:

            result = myfile.read()

            expected = (
                ">long score=0.1500 in=0.1500 ex=0.0000 len=84 ref=reference pos=0\n"
                + "ACTGAACCTTGGAAACCCTTTGGGAAAACCCCTTTTGGGGAAAAACCCCCTTTTTGGGGGAAAAAACCCCCCTTTTTTGGGGGG\n")

            self.assertEqual(result, expected)

        os.remove(filterSignatures.sortedLocation)

    """ 
    # =============================================================================

    test_multiple_signatures

    PURPOSE:
        Tests an example with multiple signatures with different locations on
        their references.

    INPUT:
        0:

    EXPECTED:
        >long1 score=0.45 in=0.45 ex=0.00 len=84 ref=reference1 pos=0
        ACTGAACCTTGGAAACCCTTTGGGAAAACCCCTTTTGGGGAAAAACCCCCTTTTTGGGGGAAAAAACCCCCCTTTTTTGGGGGG
        >long2 score=0.35 in=0.35 ex=0.00 len=84 ref=reference3 pos=100
        ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT

    # =============================================================================
    """
    def test_multiple_signatures(self):

        filterSignatures = DefaultFilterSignatures().default

        filterSignatures.filteredLocation = getPath("tests/data/filter/multiple.fasta")
        filterSignatures.sortedLocation = getPath("tests/output/filter/temp.out")
        filterSignatures.totalInclusion = 4

        hit1 = Hit("long1 84 reference1 84 80 60")
        hit2 = Hit("long1 84 reference2 84 100 84")

        hit3 = Hit("long2 84 reference3 84 60 40")
        hit4 = Hit("long2 84 reference4 84 80 60")

        filterSignatures.updatePairDictionary(filterSignatures.inclusionPairDictionary, hit1)
        filterSignatures.updatePairDictionary(filterSignatures.inclusionPairDictionary, hit2)
        filterSignatures.updatePairDictionary(filterSignatures.inclusionPairDictionary, hit3)
        filterSignatures.updatePairDictionary(filterSignatures.inclusionPairDictionary, hit4)

        # empty
        self.assertDictEqual(filterSignatures.overallScore, {})
        self.assertDictEqual(filterSignatures.inclusionScore, {})

        filterSignatures.updateInclusionScores()
  
        self.assertAlmostEqual(filterSignatures.overallScore["long1"], 0.45)
        self.assertAlmostEqual(filterSignatures.inclusionScore["long1"], 0.45)

        self.assertAlmostEqual(filterSignatures.overallScore["long2"], 0.35)
        self.assertAlmostEqual(filterSignatures.inclusionScore["long2"], 0.35)

        sortedSignatureIDs = [ID for (ID, score) in sorted(
            filterSignatures.overallScore.items(), key=operator.itemgetter(1), reverse=True)]

        filterSignatures.reportSorted(sortedSignatureIDs)

        with open (filterSignatures.sortedLocation, "r") as myfile:

            result = myfile.read()
            expected = (
                ">long1 score=0.4500 in=0.4500 ex=0.0000 len=84 ref=reference1 pos=0\n"
                + "ACTGAACCTTGGAAACCCTTTGGGAAAACCCCTTTTGGGGAAAAACCCCCTTTTTGGGGGAAAAAACCCCCCTTTTTTGGGGGG\n"
                + ">long2 score=0.3500 in=0.3500 ex=0.0000 len=84 ref=reference3 pos=100\n"
                + "ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT\n")
            self.assertEqual(result, expected)

        os.remove(filterSignatures.sortedLocation)

if __name__ == '__main__':
    
    unittest.main()
