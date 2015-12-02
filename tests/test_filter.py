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

from neptune.FilterSignatures import *
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
        0:

        >long1 84 reference1 0
        ACTGAACCTTGGAAACCCTTTGGGAAAACCCCTTTTGGGGAAAAACCCCCTTTTTGGGGGAAAAAACCCCCCTTTTTTGGGGGG
        >long2 84 reference3 100
        ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT

    EXPECTED:
        0:

        signature["long1"] -> LONG1 SIGNATURE
        signature["long2"] -> LONG2 SIGNATURE

    # =============================================================================
    """
    def test_simple(self):

        fileLocation = getPath("tests/data/filter/multiple.fasta")
        signatures = Signature.readSignatures(fileLocation)

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
        0:

        signature = Signature("ACGTACGT", "0", "ref", "20")

    EXPECTED:
        0:

        >0 8 ref 20
        ACGTACGT

    # =============================================================================
    """
    def test_simple(self):

        destination = StringIO.StringIO()
        signature = Signature("0", "ACGTACGT", "ref", "20")

        Signature.writeSignature(signature, destination)
        result = destination.getvalue()

        expected = ">0 8 ref 20\nACGTACGT\n"
        self.assertEquals(result, expected)

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

"""
# =============================================================================

UPDATE HIT OVERALL DICTIONARY

# =============================================================================
"""
class TestUpdateHitOverallDictionary(unittest.TestCase):

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

        bestOverall = {}

        query1 = "query1"
        query2 = "query2"

        # empty
        self.assertDictEqual(bestOverall, {})

        #0:
        hit1 = Hit(str(query1) + " 10 subject 20 0.50 1")
        updateHitOverallDictionary(hit1, bestOverall)
        expected = {}
        expected[query1] = hit1
        self.assertDictEqual(bestOverall, expected)

        #1:
        hit2 = Hit(str(query2) + " 10 subject 20 0.50 1")
        updateHitOverallDictionary(hit2, bestOverall)
        expected = {}
        expected[query1] = hit1
        expected[query2] = hit2
        self.assertDictEqual(bestOverall, expected)

        #2:
        hit3 = Hit(str(query1) + " 13 subject 20 0.60 3")
        updateHitOverallDictionary(hit3, bestOverall)
        expected = {}
        expected[query1] = hit3
        expected[query2] = hit2
        self.assertDictEqual(bestOverall, expected)

        #3:
        hit4 = Hit(str(query2) + " 13 subject 20 0.60 3")
        updateHitOverallDictionary(hit4, bestOverall)
        expected = {}
        expected[query1] = hit3
        expected[query2] = hit4
        self.assertDictEqual(bestOverall, expected)

"""
# =============================================================================

UPDATE BEST PAIRS

# =============================================================================
"""
class TestUpdateHitPairDictionary(unittest.TestCase):

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

        bestPairs = {}

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
        self.assertDictEqual(bestPairs, {})

        #0:
        updateHitPairDictionary(hit1, bestPairs)
        expected = {}
        expected[pair1] = hit1
        self.assertDictEqual(bestPairs, expected)

        #1:
        updateHitPairDictionary(hit2, bestPairs)
        expected = {}
        expected[pair1] = hit1
        expected[pair2] = hit2
        self.assertDictEqual(bestPairs, expected)

        #2:
        updateHitPairDictionary(hit3, bestPairs)
        expected = {}
        expected[pair1] = hit1
        expected[pair2] = hit2
        expected[pair3] = hit3
        self.assertDictEqual(bestPairs, expected)

        #3:
        updateHitPairDictionary(hit4, bestPairs)
        expected = {}
        expected[pair1] = hit1
        expected[pair2] = hit2
        expected[pair3] = hit3
        expected[pair4] = hit4
        self.assertDictEqual(bestPairs, expected)

        #4:
        updateHitPairDictionary(hit5, bestPairs)
        expected = {}
        expected[pair1] = hit5
        expected[pair2] = hit2
        expected[pair3] = hit3
        expected[pair4] = hit4
        self.assertDictEqual(bestPairs, expected)

        #5:
        updateHitPairDictionary(hit6, bestPairs)
        expected = {}
        expected[pair1] = hit5
        expected[pair2] = hit6
        expected[pair3] = hit3
        expected[pair4] = hit4
        self.assertDictEqual(bestPairs, expected)

        #6:
        updateHitPairDictionary(hit7, bestPairs)
        expected = {}
        expected[pair1] = hit5
        expected[pair2] = hit6
        expected[pair3] = hit7
        expected[pair4] = hit4
        self.assertDictEqual(bestPairs, expected)

        #7:
        updateHitPairDictionary(hit8, bestPairs)
        expected = {}
        expected[pair1] = hit5
        expected[pair2] = hit6
        expected[pair3] = hit7
        expected[pair4] = hit8
        self.assertDictEqual(bestPairs, expected)

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

        overallScore.clear()
        exclusionScore.clear()

        hit1 = Hit("query1 20 subject1 10 20 1")
        hit2 = Hit("query1 20 subject2 10 40 1")
        hit3 = Hit("query2 20 subject1 10 60 1")
        hit4 = Hit("query2 20 subject2 10 80 1")

        bestPairs = {}
        updateHitPairDictionary(hit1, bestPairs)
        updateHitPairDictionary(hit2, bestPairs)
        updateHitPairDictionary(hit3, bestPairs)
        updateHitPairDictionary(hit4, bestPairs)

        # empty
        self.assertDictEqual(overallScore, {})
        self.assertDictEqual(exclusionScore, {})

        updateExclusionScores(bestPairs, 2)
  
        self.assertAlmostEqual(overallScore["query1"], -0.15)
        self.assertAlmostEqual(overallScore["query2"], -0.35)

        self.assertAlmostEqual(exclusionScore["query1"], -0.15)
        self.assertAlmostEqual(exclusionScore["query2"], -0.35)

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

        overallScore.clear()
        exclusionScore.clear()

        hit1 = Hit("query1 10 subject1 20 10 1")
        hit2 = Hit("query1 10 subject2 20 20 1")
        hit3 = Hit("query2 10 subject1 20 30 1")
        hit4 = Hit("query2 10 subject2 20 40 1")

        bestPairs = {}
        updateHitPairDictionary(hit1, bestPairs)
        updateHitPairDictionary(hit2, bestPairs)
        updateHitPairDictionary(hit3, bestPairs)
        updateHitPairDictionary(hit4, bestPairs)

        # empty
        self.assertDictEqual(overallScore, {})
        self.assertDictEqual(exclusionScore, {})

        updateExclusionScores(bestPairs, 2)
  
        self.assertAlmostEqual(overallScore["query1"], -0.30)
        self.assertAlmostEqual(overallScore["query2"], -0.70)

        self.assertAlmostEqual(exclusionScore["query1"], -0.30)
        self.assertAlmostEqual(exclusionScore["query2"], -0.70)

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

        overallScore.clear()
        inclusionScore.clear()

        hit1 = Hit("query1 20 subject1 10 20 1")
        hit2 = Hit("query1 20 subject2 10 40 1")
        hit3 = Hit("query2 20 subject1 10 60 1")
        hit4 = Hit("query2 20 subject2 10 80 1")

        bestPairs = {}
        updateHitPairDictionary(hit1, bestPairs)
        updateHitPairDictionary(hit2, bestPairs)
        updateHitPairDictionary(hit3, bestPairs)
        updateHitPairDictionary(hit4, bestPairs)

        # empty
        self.assertDictEqual(overallScore, {})
        self.assertDictEqual(inclusionScore, {})

        updateInclusionScores(bestPairs, 2)
  
        self.assertAlmostEqual(overallScore["query1"], 0.15)
        self.assertAlmostEqual(overallScore["query2"], 0.35)

        self.assertAlmostEqual(inclusionScore["query1"], 0.15)
        self.assertAlmostEqual(inclusionScore["query2"], 0.35)

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

        overallScore.clear()
        inclusionScore.clear()

        hit1 = Hit("query1 10 subject1 20 10 1")
        hit2 = Hit("query1 10 subject2 20 20 1")
        hit3 = Hit("query2 10 subject1 20 30 1")
        hit4 = Hit("query2 10 subject2 20 40 1")

        bestPairs = {}
        updateHitPairDictionary(hit1, bestPairs)
        updateHitPairDictionary(hit2, bestPairs)
        updateHitPairDictionary(hit3, bestPairs)
        updateHitPairDictionary(hit4, bestPairs)

        # empty
        self.assertDictEqual(overallScore, {})
        self.assertDictEqual(inclusionScore, {})

        updateInclusionScores(bestPairs, 2)
  
        self.assertAlmostEqual(overallScore["query1"], 0.30)
        self.assertAlmostEqual(overallScore["query2"], 0.70)

        self.assertAlmostEqual(inclusionScore["query1"], 0.30)
        self.assertAlmostEqual(inclusionScore["query2"], 0.70)

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
        0:

    # =============================================================================
    """
    def test_simple(self):
        
        candidatesLocation = getPath("tests/data/filter/long.fasta")
        outputLocation = getPath("tests/output/filter/temp.out")

        hitOverallDictionary = {}

        filterLength = 0.5

        hit1 = Hit("long 84 subject 40 100 10")
        updateHitOverallDictionary(hit1, hitOverallDictionary)

        reportCandidates(candidatesLocation, outputLocation, hitOverallDictionary, filterLength)

        with open (outputLocation, "r") as myfile:

            result = myfile.read()
            expected = ">long 84 reference 0\nACTGAACCTTGGAAACCCTTTGGGAAAACCCCTTTTGGGGAAAAACCCCCTTTTTGGGGGAAAAAACCCCCCTTTTTTGGGGGG\n"
            self.assertEquals(result, expected)

        os.remove(outputLocation)

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
        
        candidatesLocation = getPath("tests/data/filter/long.fasta")
        outputLocation = getPath("tests/output/filter/temp.out")

        bestOverall = {}

        filterLength = 0.5

        hit1 = Hit("long 84 subject 50 100 10")
        updateHitOverallDictionary(hit1, bestOverall)

        reportCandidates(candidatesLocation, outputLocation, bestOverall, filterLength)

        with open (outputLocation, "r") as myfile:

            result = myfile.read()
            expected = ""
            self.assertEquals(result, expected)

        os.remove(outputLocation)

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
        0:

    # =============================================================================
    """
    def test_simple(self):

        filteredLocation = getPath("tests/data/filter/long.fasta")
        outputLocation = getPath("tests/output/filter/temp.out")

        overallScore.clear()
        inclusionScore.clear()
        exclusionScore.clear()

        hit1 = Hit("long 84 subject1 42 20 1")
        hit2 = Hit("long 84 subject2 42 40 1")

        bestPairs = {}
        updateHitPairDictionary(hit1, bestPairs)
        updateHitPairDictionary(hit2, bestPairs)

        # empty
        self.assertDictEqual(overallScore, {})
        self.assertDictEqual(inclusionScore, {})

        updateInclusionScores(bestPairs, 2)
  
        self.assertAlmostEqual(overallScore["long"], 0.15)
        self.assertAlmostEqual(inclusionScore["long"], 0.15)

        sortedSignatureIDs = [ID for (ID, score) in sorted(
            overallScore.items(), key=operator.itemgetter(1), reverse=True)]

        reportSorted(filteredLocation, outputLocation, sortedSignatureIDs)

        with open (outputLocation, "r") as myfile:

            result = myfile.read()

            expected = (
                ">long score=0.15 in=0.15 ex=0.00 len=84 ref=reference pos=0\n"
                + "ACTGAACCTTGGAAACCCTTTGGGAAAACCCCTTTTGGGGAAAAACCCCCTTTTTGGGGGAAAAAACCCCCCTTTTTTGGGGGG\n")

            self.assertEquals(result, expected)

        os.remove(outputLocation)

    """ 
    # =============================================================================

    test_multiple_signatures

    PURPOSE:
        Tests an example with multiple signatures with different locations on
        their references.

    INPUT:
        0:

    EXPECTED:
        0:

    # =============================================================================
    """
    def test_multiple_signatures(self):

        filteredLocation = getPath("tests/data/filter/multiple.fasta")
        outputLocation = getPath("tests/output/filter/temp.out")

        overallScore.clear()
        inclusionScore.clear()
        exclusionScore.clear()

        hit1 = Hit("long1 84 reference1 84 80 60")
        hit2 = Hit("long1 84 reference2 84 100 84")

        hit3 = Hit("long2 84 reference3 84 60 40")
        hit4 = Hit("long2 84 reference4 84 80 60")

        hitPairDictionary = {}
        updateHitPairDictionary(hit1, hitPairDictionary)
        updateHitPairDictionary(hit2, hitPairDictionary)
        updateHitPairDictionary(hit3, hitPairDictionary)
        updateHitPairDictionary(hit4, hitPairDictionary)

        # empty
        self.assertDictEqual(overallScore, {})
        self.assertDictEqual(inclusionScore, {})

        updateInclusionScores(hitPairDictionary, 4)
  
        self.assertAlmostEqual(overallScore["long1"], 0.45)
        self.assertAlmostEqual(inclusionScore["long1"], 0.45)

        self.assertAlmostEqual(overallScore["long2"], 0.35)
        self.assertAlmostEqual(inclusionScore["long2"], 0.35)

        sortedSignatureIDs = [ID for (ID, score) in sorted(
            overallScore.items(), key=operator.itemgetter(1), reverse=True)]

        reportSorted(filteredLocation, outputLocation, sortedSignatureIDs)

        with open (outputLocation, "r") as myfile:

            result = myfile.read()
            expected = (
                ">long1 score=0.45 in=0.45 ex=0.00 len=84 ref=reference1 pos=0\n"
                + "ACTGAACCTTGGAAACCCTTTGGGAAAACCCCTTTTGGGGAAAAACCCCCTTTTTGGGGGAAAAAACCCCCCTTTTTTGGGGGG\n"
                + ">long2 score=0.35 in=0.35 ex=0.00 len=84 ref=reference3 pos=100\n"
                + "ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT\n")
            self.assertEquals(result, expected)

        os.remove(outputLocation)

if __name__ == '__main__':
    
    unittest.main()
