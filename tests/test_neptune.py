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
import shutil

from TestingUtility import *
prepareSystemPath()

from neptune.Neptune import *

import unittest

class DefaultArgs():

        def __init__(self):

            parameters = {}

            parameters[CountKMers.KMER] = 5
            parameters[CountKMers.PARALLEL] = 0

            parameters[ExtractSignatures.RATE] = 0.01
            parameters[ExtractSignatures.INHITS] = 1
            parameters[ExtractSignatures.EXHITS] = 2
            parameters[ExtractSignatures.GAP] = 5
            parameters[ExtractSignatures.SIZE] = 5
            parameters[ExtractSignatures.GC_CONTENT] = 0.5
            parameters[ExtractSignatures.CONFIDENCE] = 0.95
            parameters[ExtractSignatures.INCLUSION] = getPath("tests/data/neptune/simple.fasta")
            parameters[ExtractSignatures.EXCLUSION] = getPath("tests/data/neptune/alternative.fasta")
            parameters[ExtractSignatures.REFERENCE] = ["tests/data/neptune/simple.fasta"]
            parameters[ExtractSignatures.REFERENCE_SIZE] = 12

            parameters[FilterSignatures.FILTER_LENGTH] = 0.5
            parameters[FilterSignatures.FILTER_PERCENT] = 0.5
            parameters[FilterSignatures.SEED_SIZE] = 9

            parameters[OUTPUT] = getPath("tests/output/neptune/temp.dir")

            parameters[DEFAULT_SPECIFICATION] = None
            parameters[COUNT_SPECIFICATION] = None
            parameters[AGGREGATE_SPECIFICATION] = None
            parameters[EXTRACT_SPECIFICATION] = None
            parameters[DATABASE_SPECIFICATION] = None
            parameters[FILTER_SPECIFICATION] = None
            parameters[CONSOLIDATE_SPECIFICATION] = None

            self.parameters = parameters

def setCommandLineArguments(parameters):

        sys.argv[1:] = [
            CountKMers.KMER_LONG, str(parameters[CountKMers.KMER]),
            ExtractSignatures.RATE_LONG, str(parameters[ExtractSignatures.RATE]),
            ExtractSignatures.INHITS_LONG, str(parameters[ExtractSignatures.INHITS]),
            ExtractSignatures.EXHITS_LONG, str(parameters[ExtractSignatures.EXHITS]),
            ExtractSignatures.GAP_LONG, str(parameters[ExtractSignatures.GAP]),
            ExtractSignatures.SIZE_LONG, str(parameters[ExtractSignatures.SIZE]),
            ExtractSignatures.GC_LONG, str(parameters[ExtractSignatures.GC_CONTENT]),
            FilterSignatures.FILTER_LENGTH_LONG, str(parameters[FilterSignatures.FILTER_LENGTH]),
            FilterSignatures.FILTER_PERCENT_LONG, str(parameters[FilterSignatures.FILTER_PERCENT]),
            CountKMers.PARALLEL_LONG, str(parameters[CountKMers.PARALLEL]),
            ExtractSignatures.INCLUSION_LONG, str(parameters[ExtractSignatures.INCLUSION]),
            ExtractSignatures.EXCLUSION_LONG, str(parameters[ExtractSignatures.EXCLUSION]),
            ExtractSignatures.REFERENCE_SIZE_LONG, str(parameters[ExtractSignatures.REFERENCE_SIZE]),
            OUTPUT_LONG, str(parameters[OUTPUT])]

class TestMain(unittest.TestCase):

    def test_simple(self):

        parameters = DefaultArgs().parameters

        setCommandLineArguments(parameters)
        main()

        sortedOutputDirectory = os.path.join(parameters[OUTPUT], "sorted")
        with open (os.path.join(sortedOutputDirectory, "simple.fasta"), "r") as myfile:

            myfile.readline()
            result = myfile.readline()
            expected = "ACGTACGTACGTACGTACGTACGTACGT\n"
            self.assertEquals(result, expected)

        shutil.rmtree(parameters[OUTPUT])

    def test_simple_input_directories(self):

        parameters = DefaultArgs().parameters
        parameters[ExtractSignatures.INCLUSION] = getPath("tests/data/neptune/simple_inclusion")
        parameters[ExtractSignatures.EXCLUSION] = getPath("tests/data/neptune/simple_exclusion")

        setCommandLineArguments(parameters)
        main()

        sortedOutputDirectory = os.path.join(parameters[OUTPUT], "sorted")
        with open (os.path.join(sortedOutputDirectory, "inclusion1.fasta"), "r") as myfile:

            myfile.readline()
            result = myfile.readline()
            expected = "ACGTACGTACGTACGTACGTACGTACGT\n"
            self.assertEquals(result, expected)

        shutil.rmtree(parameters[OUTPUT])

    def test_parallel(self):

        parameters = DefaultArgs().parameters
        parameters[CountKMers.PARALLEL] = 1

        setCommandLineArguments(parameters)
        main()

        sortedOutputDirectory = os.path.join(parameters[OUTPUT], "sorted")
        with open (os.path.join(sortedOutputDirectory, "simple.fasta"), "r") as myfile:

            myfile.readline()
            result = myfile.readline()
            expected = "ACGTACGTACGTACGTACGTACGTACGT\n"
            self.assertEquals(result, expected)

        shutil.rmtree(parameters[OUTPUT])

    def test_specify_reference(self):

        parameters = DefaultArgs().parameters
        parameters[ExtractSignatures.REFERENCE] = getPath("tests/data/neptune/simple.fasta")

        setCommandLineArguments(parameters)
        main()

        sortedOutputDirectory = os.path.join(parameters[OUTPUT], "sorted")
        with open (os.path.join(sortedOutputDirectory, "simple.fasta"), "r") as myfile:

            myfile.readline()
            result = myfile.readline()
            expected = "ACGTACGTACGTACGTACGTACGTACGT\n"
            self.assertEquals(result, expected)

        shutil.rmtree(parameters[OUTPUT])

    def test_invalid_k(self):

        parameters = DefaultArgs().parameters
        parameters[CountKMers.KMER] = -1

        setCommandLineArguments(parameters)

        with self.assertRaises(RuntimeError):
            main()

    def test_invalid_rate(self):

        parameters = DefaultArgs().parameters
        parameters[ExtractSignatures.RATE] = -0.01

        setCommandLineArguments(parameters)

        with self.assertRaises(RuntimeError):
            main()

        parameters = DefaultArgs().parameters
        parameters[ExtractSignatures.RATE] = -0.01

        setCommandLineArguments(parameters)

        with self.assertRaises(RuntimeError):
            main()

    def test_invalid_inhits(self):

        parameters = DefaultArgs().parameters
        parameters[ExtractSignatures.INHITS] = -1

        setCommandLineArguments(parameters)

        with self.assertRaises(RuntimeError):
            main()

    def test_invalid_exhits(self):

        parameters = DefaultArgs().parameters
        parameters[ExtractSignatures.EXHITS] = -1

        setCommandLineArguments(parameters)

        with self.assertRaises(RuntimeError):
            main()

    def test_invalid_gap(self):

        parameters = DefaultArgs().parameters
        parameters[ExtractSignatures.GAP] = -1

        setCommandLineArguments(parameters)

        with self.assertRaises(RuntimeError):
            main()

    def test_invalid_size(self):

        parameters = DefaultArgs().parameters
        parameters[ExtractSignatures.SIZE] = -1

        setCommandLineArguments(parameters)

        with self.assertRaises(RuntimeError):
            main()

    def test_invalid_gc(self):

        parameters = DefaultArgs().parameters
        parameters[ExtractSignatures.GC_CONTENT] = -0.01

        setCommandLineArguments(parameters)

        with self.assertRaises(RuntimeError):
            main()

        parameters = DefaultArgs().parameters
        parameters[ExtractSignatures.GC_CONTENT] = 1.01

        setCommandLineArguments(parameters)

        with self.assertRaises(RuntimeError):
            main()

    def test_invalid_filter_length(self):

        parameters = DefaultArgs().parameters
        parameters[FilterSignatures.FILTER_LENGTH] = -0.01

        setCommandLineArguments(parameters)

        with self.assertRaises(RuntimeError):
            main()

        parameters = DefaultArgs().parameters
        parameters[FilterSignatures.FILTER_LENGTH] = 1.01

        setCommandLineArguments(parameters)

        with self.assertRaises(RuntimeError):
            main()

    def test_invalid_filter_percent(self):

        parameters = DefaultArgs().parameters
        parameters[FilterSignatures.FILTER_PERCENT] = -0.01

        setCommandLineArguments(parameters)

        with self.assertRaises(RuntimeError):
            main()

        parameters = DefaultArgs().parameters
        parameters[FilterSignatures.FILTER_PERCENT] = 1.01

        setCommandLineArguments(parameters)

        with self.assertRaises(RuntimeError):
            main()

    def test_invalid_parallelization(self):

        parameters = DefaultArgs().parameters
        parameters[CountKMers.PARALLEL] = -1

        setCommandLineArguments(parameters)

        with self.assertRaises(RuntimeError):
            main()

    # designed to catch possible parallel race-condition failure with many files
    def test_multiple_files(self):

        inclusion = getPath("tests/data/neptune/large_inclusion")
        exclusion = getPath("tests/data/neptune/large_exclusion")
        outputDirectoryLocation = getPath("tests/output/neptune/temp.dir")
        sortedDirectoryLocation = os.path.join(outputDirectoryLocation, "sorted")

        sys.argv[1:] = [
            ExtractSignatures.INCLUSION_LONG, str(inclusion),
            ExtractSignatures.EXCLUSION_LONG, str(exclusion),
            OUTPUT_LONG, str(outputDirectoryLocation)]

        main()

        with open (os.path.join(sortedDirectoryLocation, "in0.fasta"), "r") as myfile:

            result = myfile.readline()
            expected = ">0 score=0.8889 in=1.0000 ex=-0.1111 len=99 ref=inclusion0 pos=999\n"
            self.assertEquals(result, expected)

            result = myfile.readline()
            expected = "CGGTTTCTTCATATATAACCCCGTCGGCGCTTCAGAAAACAGGGATGTATAGAATCTCTGCGTCAGAACGGCATCTAAAATCAAAACGGTATTGATGAC\n"
            self.assertEquals(result, expected)

        shutil.rmtree(outputDirectoryLocation)

if __name__ == '__main__':
    
    unittest.main()
