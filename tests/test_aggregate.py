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

from neptune.AggregateKMers import *

import unittest

"""
# =============================================================================

FIND SMALLEST

# =============================================================================
"""
class TestFindSmallest(unittest.TestCase):

    """ 
    # =============================================================================

    test_simple

    PURPOSE:
        Tests a simple use case of finding the smallest element in a list.

    INPUT:
        strings = ["a", "b", "c"], SENTINEL = "~"

    EXPECTED:
        "a"

    # =============================================================================
    """
    def test_simple(self):

        strings = ["a", "b", "c"]
        SENTINEL = "~"

        result = findSmallest(strings, SENTINEL)
        expected = "a"
        self.assertEquals(result, expected)

    """ 
    # =============================================================================

    test_empty_strings

    PURPOSE:
        Tests finding the SENTINEL when all strings in the list are empty.

    INPUT:
        strings = ["", "", ""], SENTINEL = "~"

    EXPECTED:
        "~"

    # =============================================================================
    """
    def test_empty_strings(self):

        strings = ["", "", ""]
        SENTINEL = "~"

        result = findSmallest(strings, SENTINEL)
        expected = SENTINEL
        self.assertEquals(result, expected)

    """ 
    # =============================================================================

    test_mix

    PURPOSE:
        Tests finding the smallest element in a list containing a mix of empty and
        non-empty strings.

    INPUT:
        strings = ["c", "", "a"], SENTINEL = "~"

    EXPECTED:
        "a"

    # =============================================================================
    """
    def test_mix(self):

        strings = ["c", "", "a"]
        SENTINEL = "~"

        result = findSmallest(strings, SENTINEL)
        expected = "a"
        self.assertEquals(result, expected)

"""
# =============================================================================

AGGREGATE

# =============================================================================
"""
class TestAggregate(unittest.TestCase):

    """ 
    # =============================================================================

    test_simple

    PURPOSE:
        Tests a simple k-mer aggregation use case.

    INPUT:

        IN1: aggregate1.kmers:
        AAA
        CAA
        GAA
        TAA

        IN2: aggregate2.kmers:
        AAA
        TAA

        EX1: aggregate3.kmers:
        CAA
        GAA

        EX2: aggregate4.kmers:
        TAA

    EXPECTED:

        OUT:
        AAA 2 0
        CAA 1 1
        GAA 1 1
        TAA 2 1

    # =============================================================================
    """
    def test_simple(self):

        aggregateLocation1 = "tests/data/aggregate/aggregate1.kmers"
        aggregateLocation2 = "tests/data/aggregate/aggregate2.kmers"
        aggregateLocation3 = "tests/data/aggregate/aggregate3.kmers"
        aggregateLocation4 = "tests/data/aggregate/aggregate4.kmers"

        inclusionLocations = [aggregateLocation1, aggregateLocation2]
        exclusionLocations = [aggregateLocation3, aggregateLocation4]
        outputLocation = "tests/output/aggregate/kmers.out"
        delete = False

        aggregate(inclusionLocations, exclusionLocations, outputLocation, delete)

        # check the output
        with open(outputLocation, "r") as myfile:

            result = myfile.read()
            expected = "AAA 2 0\nCAA 1 1\nGAA 1 1\nTAA 2 1\n"
            self.assertEquals(result, expected)

        os.remove(outputLocation)

    """ 
    # =============================================================================

    test_all_inclusion

    PURPOSE:
        Tests aggregation when there are inclusion files, but no exclusion files.

    INPUT:

        IN1: aggregate1.kmers:
        AAA
        CAA
        GAA
        TAA

        IN2: aggregate2.kmers:
        AAA
        TAA

    EXPECTED:

        OUT:
        AAA 2 0
        CAA 1 0
        GAA 1 0
        TAA 2 0

    # =============================================================================
    """
    def test_all_inclusion(self):

        aggregateLocation1 = "tests/data/aggregate/aggregate1.kmers"
        aggregateLocation2 = "tests/data/aggregate/aggregate2.kmers"

        inclusionLocations = [aggregateLocation1, aggregateLocation2]
        exclusionLocations = []
        outputLocation = "tests/output/aggregate/kmers.out"
        delete = False

        aggregate(inclusionLocations, exclusionLocations, outputLocation, delete)

        # check the output
        with open(outputLocation, "r") as myfile:

            result = myfile.read()
            expected = "AAA 2 0\nCAA 1 0\nGAA 1 0\nTAA 2 0\n"
            self.assertEquals(result, expected)

        os.remove(outputLocation)

    """ 
    # =============================================================================

    test_all_exclusion

    PURPOSE:
        Tests aggregation when there are exclusion files, but not inclusion files.

    INPUT:

        EX1: aggregate3.kmers:
        CAA
        GAA

        EX2: aggregate4.kmers:
        TAA

    EXPECTED:

        OUT:
        CAA 0 1
        GAA 0 1
        TAA 0 1

    # =============================================================================
    """
    def test_all_exclusion(self):

        aggregateLocation3 = "tests/data/aggregate/aggregate3.kmers"
        aggregateLocation4 = "tests/data/aggregate/aggregate4.kmers"

        inclusionLocations = []
        exclusionLocations = [aggregateLocation3, aggregateLocation4]
        outputLocation = "tests/output/aggregate/kmers.out"
        delete = False

        aggregate(inclusionLocations, exclusionLocations, outputLocation, delete)

        # check the output
        with open(outputLocation, "r") as myfile:

            result = myfile.read()
            expected = "CAA 0 1\nGAA 0 1\nTAA 0 1\n"
            self.assertEquals(result, expected)

        os.remove(outputLocation)
        

    """ 
    # =============================================================================

    test_nothing

    PURPOSE:
        Tests aggregation when there are neither inclusion or exclusion input
        files.

    INPUT:
        [NONE]

    EXPECTED:
        [EMPTY STRING] ("")

    # =============================================================================
    """
    def test_nothing(self):

        inclusionLocations = []
        exclusionLocations = []
        outputLocation = "tests/output/aggregate/kmers.out"
        delete = False

        aggregate(inclusionLocations, exclusionLocations, outputLocation, delete)

        # check the output
        with open(outputLocation, "r") as myfile:

            result = myfile.read()
            expected = ""
            self.assertEquals(result, expected)

        os.remove(outputLocation)

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
        Tests a simple use case of the entire AggregateKMers execution.

    INPUT:

        IN1: aggregate1.kmers:
        AAA
        CAA
        GAA
        TAA

        IN2: aggregate2.kmers:
        AAA
        TAA

        EX1: aggregate3.kmers:
        CAA
        GAA

        EX2: aggregate4.kmers:
        TAA

    EXPECTED:

        OUT:
        AAA 2 0
        CAA 1 1
        GAA 1 1
        TAA 2 1

    # =============================================================================
    """
    def test_simple(self):

        aggregateLocation1 = "tests/data/aggregate/aggregate1.kmers"
        aggregateLocation2 = "tests/data/aggregate/aggregate2.kmers"
        aggregateLocation3 = "tests/data/aggregate/aggregate3.kmers"
        aggregateLocation4 = "tests/data/aggregate/aggregate4.kmers"

        outputLocation = getPath("tests/output/aggregate/kmers.out")

        sys.argv[1:] = [INCLUSION_LONG, aggregateLocation1, aggregateLocation2,
            EXCLUSION_LONG, aggregateLocation3, aggregateLocation4,
            OUTPUT_LONG, outputLocation]

        main()

        # check the output
        with open(outputLocation, "r") as myfile:

            result = myfile.read()
            expected = "AAA 2 0\nCAA 1 1\nGAA 1 1\nTAA 2 1\n"
            self.assertEquals(result, expected)

        os.remove(outputLocation)

    """ 
    # =============================================================================

    test_inclusion_file_missing

    PURPOSE:
        Tests the entire AggregateKMers execution when an inclusion file is
        specified, but no such file exists on the file system.

    INPUT:

        IN1: aggregate1.kmers:
        AAA
        CAA
        GAA
        TAA

        IN2: aggregate2.kmers:
        AAA
        TAA

        IN3: MISSING.FAKE
        [DOES NOT EXIST]

        EX1: aggregate3.kmers:
        CAA
        GAA

        EX2: aggregate4.kmers:
        TAA

    EXPECTED:
        [RuntimeError]

    # =============================================================================
    """
    def test_inclusion_file_missing(self):

        aggregateLocation1 = "tests/data/aggregate/aggregate1.kmers"
        aggregateLocation2 = "tests/data/aggregate/aggregate2.kmers"
        aggregateLocation3 = "tests/data/aggregate/aggregate3.kmers"
        aggregateLocation4 = "tests/data/aggregate/aggregate4.kmers"
        MISSING = "tests/data/aggregate/MISSING.FAKE"

        outputLocation = getPath("tests/output/aggregate/kmers.out")

        sys.argv[1:] = [INCLUSION_LONG, aggregateLocation1, aggregateLocation2, MISSING,
            EXCLUSION_LONG, aggregateLocation3, aggregateLocation4,
            OUTPUT_LONG, outputLocation]

        with self.assertRaises(RuntimeError):
            main()

    """ 
    # =============================================================================

    test_exclusion_file_missing

    PURPOSE:
        Tests the entire AggregateKMers execution when an exclusion file is
        specified, but no such file exists on the file system.

    INPUT:

        IN1: aggregate1.kmers:
        AAA
        CAA
        GAA
        TAA

        IN2: aggregate2.kmers:
        AAA
        TAA

        EX1: aggregate3.kmers:
        CAA
        GAA

        EX2: aggregate4.kmers:
        TAA

        EX3: MISSING.FAKE
        [DOES NOT EXIST]

    EXPECTED:
        [RuntimeError]

    # =============================================================================
    """
    def test_exclusion_file_missing(self):

        aggregateLocation1 = "tests/data/aggregate/aggregate1.kmers"
        aggregateLocation2 = "tests/data/aggregate/aggregate2.kmers"
        aggregateLocation3 = "tests/data/aggregate/aggregate3.kmers"
        aggregateLocation4 = "tests/data/aggregate/aggregate4.kmers"
        MISSING = "tests/data/aggregate/MISSING.FAKE"

        outputLocation = getPath("tests/output/aggregate/kmers.out")

        sys.argv[1:] = [INCLUSION_LONG, aggregateLocation1, aggregateLocation2,
            EXCLUSION_LONG, aggregateLocation3, aggregateLocation4, MISSING,
            OUTPUT_LONG, outputLocation]

        with self.assertRaises(RuntimeError):
            main()

    """ 
    # =============================================================================

    test_delete_input

    PURPOSE:
        Tests the entire AggregateKMers execution with file deletion after
        execution.

    INPUT:

        IN1: tempLocation1.kmers
        AAA
        CAA

        EX1: tempLocation2.kmers
        CAA
        GAA

    EXPECTED:
        (files tempLocation1.kmers and tempLocation2.kmers do not exist)

    # =============================================================================
    """
    def test_delete_input(self):

        tempLocation1 = getPath("tests/output/aggregate/temp1.kmers")
        tempLocation2 = getPath("tests/output/aggregate/temp2.kmers")

        outputLocation = getPath("tests/output/aggregate/kmers.out")

        tempFile1 = open(tempLocation1, 'w')
        tempFile2 = open(tempLocation2, 'w')

        tempFile1.write("AAA\nCAA")
        tempFile1.close()

        tempFile2.write("CAA\nGAA")
        tempFile2.close()

        # ensure they exist before deleting
        self.assertTrue(os.path.isfile(tempLocation1))
        self.assertTrue(os.path.isfile(tempLocation2))

        sys.argv[1:] = [INCLUSION_LONG, tempLocation1,
            EXCLUSION_LONG, tempLocation2,
            OUTPUT_LONG, outputLocation,
            DELETE_LONG]

        main()

        self.assertFalse(os.path.isfile(tempLocation1))
        self.assertFalse(os.path.isfile(tempLocation2))

        os.remove(outputLocation)

if __name__ == '__main__':
    
    unittest.main()
