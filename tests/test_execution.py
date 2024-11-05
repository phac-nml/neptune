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

import unittest

import os
import sys
import io
import shutil

from tests.TestingUtility import *

from neptune.Execution import *

import neptune.JobManager
import neptune.JobManagerParallel

import neptune.Neptune
import neptune.CountKMers
import neptune.ExtractSignatures
import neptune.FilterSignatures
import neptune.Utility

class DefaultArgs():

        def __init__(self):

            parameters = {}

            parameters[CountKMers.KMER] = 5
            parameters[CountKMers.ORGANIZATION] = 0

            parameters[ExtractSignatures.RATE] = 0.01
            parameters[ExtractSignatures.INHITS] = 1
            parameters[ExtractSignatures.EXHITS] = 2
            parameters[ExtractSignatures.GAP] = 5
            parameters[ExtractSignatures.SIZE] = 5
            parameters[ExtractSignatures.GC_CONTENT] = 0.5
            parameters[ExtractSignatures.CONFIDENCE] = 0.95
            parameters[ExtractSignatures.INCLUSION] = [getPath("tests/data/execution/simple.fasta")]
            parameters[ExtractSignatures.EXCLUSION] = [getPath("tests/data/execution/alternative.fasta")]
            parameters[ExtractSignatures.REFERENCE] = ["tests/data/execution/simple.fasta"]
            parameters[ExtractSignatures.REFERENCE_SIZE] = 12

            parameters[FilterSignatures.FILTER_LENGTH] = 0.5
            parameters[FilterSignatures.FILTER_PERCENT] = 0.5
            parameters[FilterSignatures.SEED_SIZE] = 9

            parameters[Neptune.OUTPUT] = getPath("tests/output/execution/temp.dir")

            self.parameters = parameters

def buildParallelJobManager():

    outputDirectoryLocation = "tests/output/execution/temp.dir"
    logDirectoryLocation = "tests/output/execution/temp.dir"

    jobManager = neptune.JobManagerParallel.JobManagerParallel(
        outputDirectoryLocation, logDirectoryLocation)

    return jobManager

""" 
# =============================================================================

TEST EXECUTION CONSTRUCTOR

# =============================================================================
"""
class TestExecutionConstructor(unittest.TestCase):

    """ 
    # =============================================================================

    test_simple

    PURPOSE:
        Tests a simple case of Execution construction.

    INPUT:
        ParallelJobManager
        DefaultArgs().parameters

    EXPECTED:
        The execution does not encounter any errors.

    # =============================================================================
    """
    def test_simple(self):

        jobManager = buildParallelJobManager()

        parameters = DefaultArgs().parameters

        execution = Execution(jobManager, parameters)

    """ 
    # =============================================================================

    test_no_inclusion

    PURPOSE:
        Tests Execution construction when the inclusion targets are missing.

    INPUT:
        ParallelJobManager
        DefaultArgs().parameters

        parameters[ExtractSignatures.INCLUSION] = None

    EXPECTED:
        A RuntimeError occurs.

    # =============================================================================
    """
    def test_no_inclusion(self):

        jobManager = buildParallelJobManager()   

        parameters = DefaultArgs().parameters
        parameters[ExtractSignatures.INCLUSION] = None

        with self.assertRaises(RuntimeError):
            execution = Execution(jobManager, parameters)

    """ 
    # =============================================================================

    test_no_exclusion

    PURPOSE:
        Tests Execution construction when the exclusion targets are missing.

    INPUT:
        ParallelJobManager
        DefaultArgs().parameters

        parameters[ExtractSignatures.EXCLUSION] = None

    EXPECTED:
        A RuntimeError occurs.

    # =============================================================================
    """
    def test_no_exclusion(self):

        jobManager = buildParallelJobManager()

        parameters = DefaultArgs().parameters
        parameters[ExtractSignatures.EXCLUSION] = None

        with self.assertRaises(RuntimeError):
            execution = Execution(jobManager, parameters)

    """ 
    # =============================================================================

    test_no_output

    PURPOSE:
        Tests Execution construction when the output is None.

    INPUT:
        ParallelJobManager
        DefaultArgs().parameters

        parameters[Neptune.OUTPUT] = None

    EXPECTED:
        A RuntimeError occurs.

    # =============================================================================
    """
    def test_no_output(self):

        jobManager = buildParallelJobManager() 

        parameters = DefaultArgs().parameters
        parameters[Neptune.OUTPUT] = None

        with self.assertRaises(RuntimeError):
            execution = Execution(jobManager, parameters)


""" 
# =============================================================================

CALCULATE EXPECTED K-MER HITS

# =============================================================================
"""
class CalculateExpectedKMerHits(unittest.TestCase):

    """ 
    # =============================================================================

    test_simple

    PURPOSE:
        Tests a simple use case of calculating the expected number of k-mer hits.

    INPUT:

        gc-content = 0.50
        length = 10000
        k-mer size = 11

    EXPECTED:
        expected = 11.8959

        ((w -k + 1) choose (2)) * (2* ((1 - L) / 2)^2 + 2 * (L / 2)^2)^k 

    # =============================================================================
    """
    def test_simple(self):

        result = Execution.calculateExpectedKMerHits(0.50, 10000, 11)
        expected = 11.8959

        self.assertAlmostEqual(result, expected, 4)

    """ 
    # =============================================================================

    test_low_gc

    PURPOSE:
        Tests when the GC-content is low.

    INPUT:

        gc-content = 0.01
        length = 10000
        k-mer size = 11

    EXPECTED:
        expected = 19551.9

        ((w -k + 1) choose (2)) * (2* ((1 - L) / 2)^2 + 2 * (L / 2)^2)^k 

    # =============================================================================
    """
    def test_low_gc(self):

        result = Execution.calculateExpectedKMerHits(0.01, 10000, 11)
        expected = 19551.9

        self.assertAlmostEqual(result, expected, 1)

    """ 
    # =============================================================================

    test_high_gc

    PURPOSE:
        Tests when the GC-content is high.

    INPUT:

        gc-content = 0.99
        length = 10000
        k-mer size = 11

    EXPECTED:
        expected = 19551.9

        ((w -k + 1) choose (2)) * (2* ((1 - L) / 2)^2 + 2 * (L / 2)^2)^k 

    # =============================================================================
    """
    def test_high_gc(self):

        result = Execution.calculateExpectedKMerHits(0.99, 10000, 11)
        expected = 19551.9

        self.assertAlmostEqual(result, expected, 1)

    """ 
    # =============================================================================

    test_short_genome_length

    PURPOSE:
        Tests when the genome size is short.

    INPUT:

        gc-content = 0.50
        length = 12
        k-mer size = 11

    EXPECTED:
        expected = 0.000000238419

        ((w -k + 1) choose (2)) * (2* ((1 - L) / 2)^2 + 2 * (L / 2)^2)^k 

    # =============================================================================
    """
    def test_short_genome_length(self):

        result = Execution.calculateExpectedKMerHits(0.50, 12, 11)
        expected = 0.000000238419

        self.assertAlmostEqual(result, expected, 12)

    """ 
    # =============================================================================

    test_long_genome_length

    PURPOSE:
        Tests when the genome size is long.

    INPUT:

        gc-content = 0.50
        length = 1000000000
        k-mer size = 11

    EXPECTED:
        expected = 119209000000

        ((w -k + 1) choose (2)) * (2* ((1 - L) / 2)^2 + 2 * (L / 2)^2)^k 

    # =============================================================================
    """
    def test_long_genome_length(self):

        result = Execution.calculateExpectedKMerHits(0.50, 1000000000, 11)
        expected = 119209287047.3861825466

        self.assertAlmostEqual(result / math.pow(10, 10), expected / math.pow(10, 10), 4)

    """ 
    # =============================================================================

    test_small_kmer_size

    PURPOSE:
        Tests when the k-mer size is small.

    INPUT:

        gc-content = 0.50
        length = 10000
        k-mer size = 2

    EXPECTED:
        expected = 3124060

        ((w -k + 1) choose (2)) * (2* ((1 - L) / 2)^2 + 2 * (L / 2)^2)^k 

    # =============================================================================
    """
    def test_small_kmer_size(self):

        result = Execution.calculateExpectedKMerHits(0.50, 10000, 2)
        expected = 3124062.6

        self.assertAlmostEqual(result, expected, 1)

    """ 
    # =============================================================================

    test_large_kmer_size

    PURPOSE:
        Tests when the k-mer size is large.

    INPUT:

        gc-content = 0.50
        length = 10000
        k-mer size = 101

    EXPECTED:
        expected = 0

        ((w -k + 1) choose (2)) * (2* ((1 - L) / 2)^2 + 2 * (L / 2)^2)^k 

    # =============================================================================
    """
    def test_large_kmer_size(self):

        result = Execution.calculateExpectedKMerHits(0.50, 10000, 101)
        expected = 0

        self.assertAlmostEqual(result, expected, 2)


if __name__ == '__main__':
    
    unittest.main()
