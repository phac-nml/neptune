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
import io

from tests.TestingUtility import *

from neptune.CountKMers import *
from neptune.Utility import *

import unittest

""" 
# =============================================================================

WRITE SINGLE FILE

# =============================================================================
"""
class TestWriteSingleFile(unittest.TestCase):

    """ 
    # =============================================================================

    test_simple

    PURPOSE:
        Tests a simple use case of writing a single k-mer file.

    INPUT:
        kmers = [["AAA", 1], ["CAA", 2], ["GAA", 1], ["TAA", 3]]

    EXPECTED:

        OUTPUT:
        AAA
        CAA
        GAA
        TAA

    # =============================================================================
    """
    def test_simple(self):

        kmers = [["AAA", 1], ["CAA", 2], ["GAA", 1], ["TAA", 3]]

        buff = io.StringIO()
        writeSingleFile(kmers, buff)

        result = buff.getvalue()
        expected = "AAA 1\nCAA 2\nGAA 1\nTAA 3\n"

        self.assertEqual(result, expected)

""" 
# =============================================================================

COUNT

# =============================================================================
"""
class TestCount(unittest.TestCase):

    """ 
    # =============================================================================

    test_simple

    PURPOSE:
        Tests a simple use case of counting k-mers, while writing to one file.

    INPUT:

        ount1.fasta:
        >0
        ACGTACGTACGT
        
         k = 7

    EXPECTED:

        count1.kmers:
        ACGTACG 4
        GTACGTA 2

    # =============================================================================
    """
    def test_simple(self):

        inputLocation = "tests/data/count/count1.fasta"
        outputLocation = getPath("tests/output/count/count1.kmers")
        k = 7
        parallelization = 0

        count(inputLocation, outputLocation, k, parallelization)

        with open (outputLocation, "r") as myfile:

            result = myfile.read()
            expected = "ACGTACG 4\nGTACGTA 2\n"
            self.assertEqual(result, expected)

        os.remove(outputLocation)

    """ 
    # =============================================================================

    test_simple_multiple

    PURPOSE:
        Tests a simple use case of counting k-mers, while writing to multuple files.

    INPUT:

        count1.fasta:
        >0
        ACGTACGTACGT
        
         k = 7

    EXPECTED:

        count1.kmers.A:
        ACGTACG 4

        count1.kmers.G
        GTACGTA 2

    # =============================================================================
    """
    def test_simple_multiple(self):

        inputLocation = "tests/data/count/count1.fasta"
        outputLocation = getPath("tests/output/count/count1.kmers")
        k = 7
        parallelization = 1

        count(inputLocation, outputLocation, k, parallelization)

        with open (getPath("tests/output/count/count1.kmers.A"), "r") as myfile:

            result = myfile.read()
            expected = "ACGTACG 4\n"
            self.assertEqual(result, expected)

        with open (getPath("tests/output/count/count1.kmers.G"), "r") as myfile:

            result = myfile.read()
            expected = "GTACGTA 2\n"
            self.assertEqual(result, expected)

        tags = Utility.getAggregationTags(parallelization)

        for tag in tags:

            outputName = outputLocation + "." + tag
            os.remove(outputName)

""" 
# =============================================================================

WRITE MULTIPLE FILES

# =============================================================================
"""
class TestWriteMultiple(unittest.TestCase):

    """ 
    # =============================================================================

    test_simple

    PURPOSE:
        Tests a simple use case of writing to multiple files.

    INPUT:
        kmers = [["AAA", 1], ["CAA", 2], ["ACA", 1], ["CCA", 3], ["RCA", 3]]

    EXPECTED:

        kmers.out.A:
        AAA 1
        ACA 1

        kmers.out.C:
        CAA 2
        CCA 3

        kmers.out.G:

        kmers.out.T:

        kmers.out.__OTHER__:
        RCA 3

    # =============================================================================
    """
    def test_simple(self):

        kmers = [["AAA", 1], ["CAA", 2], ["ACA", 1], ["CCA", 3], ["RCA", 3]]

        outputLocation = getPath("tests/output/count/kmers.out")
        parallelization = 1

        # run function
        writeMultipleFiles(kmers, outputLocation, parallelization)

        # verify
        outputFiles = {}
        tags = Utility.getAggregationTags(parallelization)

        for tag in tags:

            outputName = outputLocation + "." + tag
            outputFiles[tag] = open(outputName, 'r')
            self.assertTrue(outputFiles[tag])

        # check output
        result = outputFiles["A"].read()
        expected = "AAA 1\nACA 1\n"
        self.assertEqual(result, expected)

        result = outputFiles["C"].read()
        expected = "CAA 2\nCCA 3\n"
        self.assertEqual(result, expected)

        result = outputFiles["G"].read()
        expected = ""
        self.assertEqual(result, expected)

        result = outputFiles["T"].read()
        expected = ""
        self.assertEqual(result, expected)        

        result = outputFiles[Utility.AGGREGATE_OTHER].read()
        expected = "RCA 3\n"
        self.assertEqual(result, expected) 

        for tag in tags:

            outputName = outputLocation + "." + tag
            os.remove(outputName)

""" 
# =============================================================================

MAIN

# =============================================================================
"""
class TestMain(unittest.TestCase):

    """ 
    # =============================================================================

    test_simple

    PURPOSE:
        A simple test of the main function and entire execution of CountKMers.

    INPUT:

        count1.kmers:
        >0
        ACGTACGTACGT

        output: "tests/data/count/count1.fasta"

        k = 7

        p = 0

    EXPECTED:

        count1.kmers:
        ACGTACG 4
        GTACGTA 2

    # =============================================================================
    """
    def test_simple(self):

        outputLocation = getPath("tests/output/count/count1.kmers")

        sys.argv[1:] = ["-i", "tests/data/count/count1.fasta", "-o", outputLocation, KMER_LONG, "7", ORGANIZATION_LONG, "0"]
        main()

        with open (outputLocation, "r") as myfile:

            result = myfile.read()
            expected = "ACGTACG 4\nGTACGTA 2\n"
            self.assertEqual(result, expected)

        os.remove(outputLocation)

    """ 
    # =============================================================================

    test_missing_input

    PURPOSE:
        Tests behaviour when the input is missing.

    INPUT:

        input: tests/data/count/this_file_does_not_exist.fake
        

    EXPECTED:
        [RuntimeError]

    # =============================================================================
    """
    def test_missing_input(self):

        sys.argv[1:] = ["-i", "tests/data/this_file_does_not_exist.fake", "-o",
            getPath("tests/output/kmers.out"), KMER_LONG, "7", ORGANIZATION_LONG, "0"]
        
        with self.assertRaises(RuntimeError):
            main()

if __name__ == '__main__':
    
    unittest.main()
