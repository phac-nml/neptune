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

class TestMain(unittest.TestCase):

    def test_simple(self):

        kmer = 5
        rate = 0.01
        inhits = 1
        exhits = 2
        gap = 5
        size = 4
        gcContent = 0.5
        filterLength = 0.5
        filterPercent = 0.5
        parallelization = 0
        inclusion = getPath("tests/data/neptune/simple.fasta")
        exclusion = getPath("tests/data/neptune/alternative.fasta")
        referenceSize = 12
        outputDirectoryLocation = getPath("tests/output/neptune/temp.dir")
        sortedDirectoryLocation = os.path.join(outputDirectoryLocation, "sorted")

        sys.argv[1:] = [
            CountKMers.KMER_LONG, str(kmer),
            ExtractSignatures.RATE_LONG, str(rate),
            ExtractSignatures.INHITS_LONG, str(inhits),
            ExtractSignatures.EXHITS_LONG, str(exhits),
            ExtractSignatures.GAP_LONG, str(gap),
            ExtractSignatures.SIZE_LONG, str(size),
            ExtractSignatures.GC_LONG, str(gcContent),
            FilterSignatures.FILTER_LENGTH_LONG, str(filterLength),
            FilterSignatures.FILTER_PERCENT_LONG, str(filterPercent),
            CountKMers.PARALLEL_LONG, str(parallelization),
            ExtractSignatures.INCLUSION_LONG, str(inclusion),
            ExtractSignatures.EXCLUSION_LONG, str(exclusion),
            ExtractSignatures.REFERENCE_SIZE_LONG, str(referenceSize),
            OUTPUT_LONG, str(outputDirectoryLocation)]

        main()

        with open (os.path.join(sortedDirectoryLocation, "simple.fasta"), "r") as myfile:

            myfile.readline()
            result = myfile.readline()
            expected = "ACGTACGTACGTACGTACGTACGTACGT\n"
            self.assertEquals(result, expected)

        shutil.rmtree(outputDirectoryLocation)

    def test_simple_directories(self):

        kmer = 5
        rate = 0.01
        inhits = 1
        exhits = 2
        gap = 5
        size = 4
        gcContent = 0.5
        filterLength = 0.5
        filterPercent = 0.5
        parallelization = 0
        inclusion = getPath("tests/data/neptune/simple_inclusion")
        exclusion = getPath("tests/data/neptune/simple_exclusion")
        referenceSize = 12
        outputDirectoryLocation = getPath("tests/output/neptune/temp.dir")
        sortedDirectoryLocation = os.path.join(outputDirectoryLocation, "sorted")

        sys.argv[1:] = [
            CountKMers.KMER_LONG, str(kmer),
            ExtractSignatures.RATE_LONG, str(rate),
            ExtractSignatures.INHITS_LONG, str(inhits),
            ExtractSignatures.EXHITS_LONG, str(exhits),
            ExtractSignatures.GAP_LONG, str(gap),
            ExtractSignatures.SIZE_LONG, str(size),
            ExtractSignatures.GC_LONG, str(gcContent),
            FilterSignatures.FILTER_LENGTH_LONG, str(filterLength),
            FilterSignatures.FILTER_PERCENT_LONG, str(filterPercent),
            CountKMers.PARALLEL_LONG, str(parallelization),
            ExtractSignatures.INCLUSION_LONG, str(inclusion),
            ExtractSignatures.EXCLUSION_LONG, str(exclusion),
            ExtractSignatures.REFERENCE_SIZE_LONG, str(referenceSize),
            OUTPUT_LONG, str(outputDirectoryLocation)]

        main()

        with open (os.path.join(sortedDirectoryLocation, "inclusion1.fasta"), "r") as myfile:

            myfile.readline()
            result = myfile.readline()
            expected = "ACGTACGTACGTACGTACGTACGTACGT\n"
            self.assertEquals(result, expected)

        shutil.rmtree(outputDirectoryLocation)

    def test_parallel(self):

        kmer = 5
        rate = 0.01
        inhits = 1
        exhits = 2
        gap = 5
        size = 4
        gcContent = 0.5
        filterLength = 0.5
        filterPercent = 0.5
        parallelization = 1
        inclusion = getPath("tests/data/neptune/simple.fasta")
        exclusion = getPath("tests/data/neptune/alternative.fasta")
        referenceSize = 12
        outputDirectoryLocation = getPath("tests/output/neptune/temp.dir")
        sortedDirectoryLocation = os.path.join(outputDirectoryLocation, "sorted")

        sys.argv[1:] = [
            CountKMers.KMER_LONG, str(kmer),
            ExtractSignatures.RATE_LONG, str(rate),
            ExtractSignatures.INHITS_LONG, str(inhits),
            ExtractSignatures.EXHITS_LONG, str(exhits),
            ExtractSignatures.GAP_LONG, str(gap),
            ExtractSignatures.SIZE_LONG, str(size),
            ExtractSignatures.GC_LONG, str(gcContent),
            FilterSignatures.FILTER_LENGTH_LONG, str(filterLength),
            FilterSignatures.FILTER_PERCENT_LONG, str(filterPercent),
            CountKMers.PARALLEL_LONG, str(parallelization),
            ExtractSignatures.INCLUSION_LONG, str(inclusion),
            ExtractSignatures.EXCLUSION_LONG, str(exclusion),
            ExtractSignatures.REFERENCE_SIZE_LONG, str(referenceSize),
            OUTPUT_LONG, str(outputDirectoryLocation)]

        main()

        with open (os.path.join(sortedDirectoryLocation, "simple.fasta"), "r") as myfile:

            myfile.readline()
            result = myfile.readline()
            expected = "ACGTACGTACGTACGTACGTACGTACGT\n"
            self.assertEquals(result, expected)

        shutil.rmtree(outputDirectoryLocation)

    def test_specify_reference(self):

        kmer = 5
        rate = 0.01
        inhits = 1
        exhits = 2
        gap = 5
        size = 4
        gcContent = 0.5
        filterLength = 0.5
        filterPercent = 0.5
        parallelization = 0
        inclusion = getPath("tests/data/neptune/simple.fasta")
        exclusion = getPath("tests/data/neptune/alternative.fasta")
        referenceSize = 12
        outputDirectoryLocation = getPath("tests/output/neptune/temp.dir")
        sortedDirectoryLocation = os.path.join(outputDirectoryLocation, "sorted")
        reference = getPath("tests/data/neptune/simple.fasta")

        sys.argv[1:] = [
            ExtractSignatures.REFERENCE_LONG, str(reference),
            CountKMers.KMER_LONG, str(kmer),
            ExtractSignatures.RATE_LONG, str(rate),
            ExtractSignatures.INHITS_LONG, str(inhits),
            ExtractSignatures.EXHITS_LONG, str(exhits),
            ExtractSignatures.GAP_LONG, str(gap),
            ExtractSignatures.SIZE_LONG, str(size),
            ExtractSignatures.GC_LONG, str(gcContent),
            FilterSignatures.FILTER_LENGTH_LONG, str(filterLength),
            FilterSignatures.FILTER_PERCENT_LONG, str(filterPercent),
            CountKMers.PARALLEL_LONG, str(parallelization),
            ExtractSignatures.INCLUSION_LONG, str(inclusion),
            ExtractSignatures.EXCLUSION_LONG, str(exclusion),
            ExtractSignatures.REFERENCE_SIZE_LONG, str(referenceSize),
            OUTPUT_LONG, str(outputDirectoryLocation)]

        main()

        with open (os.path.join(sortedDirectoryLocation, "simple.fasta"), "r") as myfile:

            myfile.readline()
            result = myfile.readline()
            expected = "ACGTACGTACGTACGTACGTACGTACGT\n"
            self.assertEquals(result, expected)

        shutil.rmtree(outputDirectoryLocation)

    def test_invalid_k(self):

        kmer = -1
        rate = 0.01
        inhits = 1
        exhits = 2
        gap = 5
        size = 7
        gcContent = 0.5
        filterLength = 0.5
        filterPercent = 0.5
        parallelization = 0
        inclusion = getPath("tests/data/neptune/simple.fasta")
        exclusion = getPath("tests/data/neptune/alternative.fasta")
        referenceSize = 12
        outputDirectoryLocation = getPath("tests/output/neptune/temp.dir")

        sys.argv[1:] = [
            CountKMers.KMER_LONG, str(kmer),
            ExtractSignatures.RATE_LONG, str(rate),
            ExtractSignatures.INHITS_LONG, str(inhits),
            ExtractSignatures.EXHITS_LONG, str(exhits),
            ExtractSignatures.GAP_LONG, str(gap),
            ExtractSignatures.SIZE_LONG, str(size),
            ExtractSignatures.GC_LONG, str(gcContent),
            FilterSignatures.FILTER_LENGTH_LONG, str(filterLength),
            FilterSignatures.FILTER_PERCENT_LONG, str(filterPercent),
            CountKMers.PARALLEL_LONG, str(parallelization),
            ExtractSignatures.INCLUSION_LONG, str(inclusion),
            ExtractSignatures.EXCLUSION_LONG, str(exclusion),
            ExtractSignatures.REFERENCE_SIZE_LONG, str(referenceSize),
            OUTPUT_LONG, str(outputDirectoryLocation)]

        with self.assertRaises(RuntimeError):
            main()

    def test_invalid_rate(self):

        kmer = 5
        rate = -0.01
        inhits = 1
        exhits = 2
        gap = 5
        size = 7
        gcContent = 0.5
        filterLength = 0.5
        filterPercent = 0.5
        parallelization = 0
        inclusion = getPath("tests/data/simple.fasta")
        exclusion = getPath("tests/data/alternative.fasta")
        referenceSize = 12
        outputDirectoryLocation = getPath("tests/output/temp.dir")

        sys.argv[1:] = [
            CountKMers.KMER_LONG, str(kmer),
            ExtractSignatures.RATE_LONG, str(rate),
            ExtractSignatures.INHITS_LONG, str(inhits),
            ExtractSignatures.EXHITS_LONG, str(exhits),
            ExtractSignatures.GAP_LONG, str(gap),
            ExtractSignatures.SIZE_LONG, str(size),
            ExtractSignatures.GC_LONG, str(gcContent),
            FilterSignatures.FILTER_LENGTH_LONG, str(filterLength),
            FilterSignatures.FILTER_PERCENT_LONG, str(filterPercent),
            CountKMers.PARALLEL_LONG, str(parallelization),
            ExtractSignatures.INCLUSION_LONG, str(inclusion),
            ExtractSignatures.EXCLUSION_LONG, str(exclusion),
            ExtractSignatures.REFERENCE_SIZE_LONG, str(referenceSize),
            OUTPUT_LONG, str(outputDirectoryLocation)]

        with self.assertRaises(RuntimeError):
            main()

        rate = 1.01

        sys.argv[1:] = [
            CountKMers.KMER_LONG, str(kmer),
            ExtractSignatures.RATE_LONG, str(rate),
            ExtractSignatures.INHITS_LONG, str(inhits),
            ExtractSignatures.EXHITS_LONG, str(exhits),
            ExtractSignatures.GAP_LONG, str(gap),
            ExtractSignatures.SIZE_LONG, str(size),
            ExtractSignatures.GC_LONG, str(gcContent),
            FilterSignatures.FILTER_LENGTH_LONG, str(filterLength),
            FilterSignatures.FILTER_PERCENT_LONG, str(filterPercent),
            CountKMers.PARALLEL_LONG, str(parallelization),
            ExtractSignatures.INCLUSION_LONG, str(inclusion),
            ExtractSignatures.EXCLUSION_LONG, str(exclusion),
            ExtractSignatures.REFERENCE_SIZE_LONG, str(referenceSize),
            OUTPUT_LONG, str(outputDirectoryLocation)]

        with self.assertRaises(RuntimeError):
            main()

    def test_invalid_inhits(self):

        kmer = 5
        rate = 0.01
        inhits = -1
        exhits = 2
        gap = 5
        size = 7
        gcContent = 0.5
        filterLength = 0.5
        filterPercent = 0.5
        parallelization = 0
        inclusion = getPath("tests/data/simple.fasta")
        exclusion = getPath("tests/data/alternative.fasta")
        referenceSize = 12
        outputDirectoryLocation = getPath("tests/output/temp.dir")

        sys.argv[1:] = [
            CountKMers.KMER_LONG, str(kmer),
            ExtractSignatures.RATE_LONG, str(rate),
            ExtractSignatures.INHITS_LONG, str(inhits),
            ExtractSignatures.EXHITS_LONG, str(exhits),
            ExtractSignatures.GAP_LONG, str(gap),
            ExtractSignatures.SIZE_LONG, str(size),
            ExtractSignatures.GC_LONG, str(gcContent),
            FilterSignatures.FILTER_LENGTH_LONG, str(filterLength),
            FilterSignatures.FILTER_PERCENT_LONG, str(filterPercent),
            CountKMers.PARALLEL_LONG, str(parallelization),
            ExtractSignatures.INCLUSION_LONG, str(inclusion),
            ExtractSignatures.EXCLUSION_LONG, str(exclusion),
            ExtractSignatures.REFERENCE_SIZE_LONG, str(referenceSize),
            OUTPUT_LONG, str(outputDirectoryLocation)]

        with self.assertRaises(RuntimeError):
            main()

    def test_invalid_exhits(self):

        kmer = 5
        rate = 0.01
        inhits = 1
        exhits = -2
        gap = 5
        size = 7
        gcContent = 0.5
        filterLength = 0.5
        filterPercent = 0.5
        parallelization = 0
        inclusion = getPath("tests/data/simple.fasta")
        exclusion = getPath("tests/data/alternative.fasta")
        referenceSize = 12
        outputDirectoryLocation = getPath("tests/output/temp.dir")

        sys.argv[1:] = [
            CountKMers.KMER_LONG, str(kmer),
            ExtractSignatures.RATE_LONG, str(rate),
            ExtractSignatures.INHITS_LONG, str(inhits),
            ExtractSignatures.EXHITS_LONG, str(exhits),
            ExtractSignatures.GAP_LONG, str(gap),
            ExtractSignatures.SIZE_LONG, str(size),
            ExtractSignatures.GC_LONG, str(gcContent),
            FilterSignatures.FILTER_LENGTH_LONG, str(filterLength),
            FilterSignatures.FILTER_PERCENT_LONG, str(filterPercent),
            CountKMers.PARALLEL_LONG, str(parallelization),
            ExtractSignatures.INCLUSION_LONG, str(inclusion),
            ExtractSignatures.EXCLUSION_LONG, str(exclusion),
            ExtractSignatures.REFERENCE_SIZE_LONG, str(referenceSize),
            OUTPUT_LONG, str(outputDirectoryLocation)]

        with self.assertRaises(RuntimeError):
            main()

    def test_invalid_gap(self):

        kmer = 5
        rate = 0.01
        inhits = 1
        exhits = 2
        gap = -1
        size = 7
        gcContent = 0.5
        filterLength = 0.5
        filterPercent = 0.5
        parallelization = 0
        inclusion = getPath("tests/data/simple.fasta")
        exclusion = getPath("tests/data/alternative.fasta")
        referenceSize = 12
        outputDirectoryLocation = getPath("tests/output/temp.dir")

        sys.argv[1:] = [
            CountKMers.KMER_LONG, str(kmer),
            ExtractSignatures.RATE_LONG, str(rate),
            ExtractSignatures.INHITS_LONG, str(inhits),
            ExtractSignatures.EXHITS_LONG, str(exhits),
            ExtractSignatures.GAP_LONG, str(gap),
            ExtractSignatures.SIZE_LONG, str(size),
            ExtractSignatures.GC_LONG, str(gcContent),
            FilterSignatures.FILTER_LENGTH_LONG, str(filterLength),
            FilterSignatures.FILTER_PERCENT_LONG, str(filterPercent),
            CountKMers.PARALLEL_LONG, str(parallelization),
            ExtractSignatures.INCLUSION_LONG, str(inclusion),
            ExtractSignatures.EXCLUSION_LONG, str(exclusion),
            ExtractSignatures.REFERENCE_SIZE_LONG, str(referenceSize),
            OUTPUT_LONG, str(outputDirectoryLocation)]

        with self.assertRaises(RuntimeError):
            main()

    def test_invalid_size(self):

        kmer = 5
        rate = 0.01
        inhits = 1
        exhits = 2
        gap = 5
        size = -1
        gcContent = 0.5
        filterLength = 0.5
        filterPercent = 0.5
        parallelization = 0
        inclusion = getPath("tests/data/simple.fasta")
        exclusion = getPath("tests/data/alternative.fasta")
        referenceSize = 12
        outputDirectoryLocation = getPath("tests/output/temp.dir")

        sys.argv[1:] = [
            CountKMers.KMER_LONG, str(kmer),
            ExtractSignatures.RATE_LONG, str(rate),
            ExtractSignatures.INHITS_LONG, str(inhits),
            ExtractSignatures.EXHITS_LONG, str(exhits),
            ExtractSignatures.GAP_LONG, str(gap),
            ExtractSignatures.SIZE_LONG, str(size),
            ExtractSignatures.GC_LONG, str(gcContent),
            FilterSignatures.FILTER_LENGTH_LONG, str(filterLength),
            FilterSignatures.FILTER_PERCENT_LONG, str(filterPercent),
            CountKMers.PARALLEL_LONG, str(parallelization),
            ExtractSignatures.INCLUSION_LONG, str(inclusion),
            ExtractSignatures.EXCLUSION_LONG, str(exclusion),
            ExtractSignatures.REFERENCE_SIZE_LONG, str(referenceSize),
            OUTPUT_LONG, str(outputDirectoryLocation)]

        with self.assertRaises(RuntimeError):
            main()

    def test_invalid_gc(self):

        kmer = 5
        rate = 0.01
        inhits = 1
        exhits = 2
        gap = 5
        size = 7
        gcContent = -0.01
        filterLength = 0.5
        filterPercent = 0.5
        parallelization = 0
        inclusion = getPath("tests/data/simple.fasta")
        exclusion = getPath("tests/data/alternative.fasta")
        referenceSize = 12
        outputDirectoryLocation = getPath("tests/output/temp.dir")

        sys.argv[1:] = [
            CountKMers.KMER_LONG, str(kmer),
            ExtractSignatures.RATE_LONG, str(rate),
            ExtractSignatures.INHITS_LONG, str(inhits),
            ExtractSignatures.EXHITS_LONG, str(exhits),
            ExtractSignatures.GAP_LONG, str(gap),
            ExtractSignatures.SIZE_LONG, str(size),
            ExtractSignatures.GC_LONG, str(gcContent),
            FilterSignatures.FILTER_LENGTH_LONG, str(filterLength),
            FilterSignatures.FILTER_PERCENT_LONG, str(filterPercent),
            CountKMers.PARALLEL_LONG, str(parallelization),
            ExtractSignatures.INCLUSION_LONG, str(inclusion),
            ExtractSignatures.EXCLUSION_LONG, str(exclusion),
            ExtractSignatures.REFERENCE_SIZE_LONG, str(referenceSize),
            OUTPUT_LONG, str(outputDirectoryLocation)]

        with self.assertRaises(RuntimeError):
            main()

        gcContent = 1.01

        sys.argv[1:] = [
            CountKMers.KMER_LONG, str(kmer),
            ExtractSignatures.RATE_LONG, str(rate),
            ExtractSignatures.INHITS_LONG, str(inhits),
            ExtractSignatures.EXHITS_LONG, str(exhits),
            ExtractSignatures.GAP_LONG, str(gap),
            ExtractSignatures.SIZE_LONG, str(size),
            ExtractSignatures.GC_LONG, str(gcContent),
            FilterSignatures.FILTER_LENGTH_LONG, str(filterLength),
            FilterSignatures.FILTER_PERCENT_LONG, str(filterPercent),
            CountKMers.PARALLEL_LONG, str(parallelization),
            ExtractSignatures.INCLUSION_LONG, str(inclusion),
            ExtractSignatures.EXCLUSION_LONG, str(exclusion),
            ExtractSignatures.REFERENCE_SIZE_LONG, str(referenceSize),
            OUTPUT_LONG, str(outputDirectoryLocation)]

        with self.assertRaises(RuntimeError):
            main()

    def test_invalid_filter_length(self):

        kmer = 5
        rate = 0.01
        inhits = 1
        exhits = 2
        gap = 5
        size = 7
        gcContent = 0.05
        filterLength = -0.01
        filterPercent = 0.5
        parallelization = 0
        inclusion = getPath("tests/data/simple.fasta")
        exclusion = getPath("tests/data/alternative.fasta")
        referenceSize = 12
        outputDirectoryLocation = getPath("tests/output/temp.dir")

        sys.argv[1:] = [
            CountKMers.KMER_LONG, str(kmer),
            ExtractSignatures.RATE_LONG, str(rate),
            ExtractSignatures.INHITS_LONG, str(inhits),
            ExtractSignatures.EXHITS_LONG, str(exhits),
            ExtractSignatures.GAP_LONG, str(gap),
            ExtractSignatures.SIZE_LONG, str(size),
            ExtractSignatures.GC_LONG, str(gcContent),
            FilterSignatures.FILTER_LENGTH_LONG, str(filterLength),
            FilterSignatures.FILTER_PERCENT_LONG, str(filterPercent),
            CountKMers.PARALLEL_LONG, str(parallelization),
            ExtractSignatures.INCLUSION_LONG, str(inclusion),
            ExtractSignatures.EXCLUSION_LONG, str(exclusion),
            ExtractSignatures.REFERENCE_SIZE_LONG, str(referenceSize),
            OUTPUT_LONG, str(outputDirectoryLocation)]

        with self.assertRaises(RuntimeError):
            main()

        filterLength = 1.01

        sys.argv[1:] = [
            CountKMers.KMER_LONG, str(kmer),
            ExtractSignatures.RATE_LONG, str(rate),
            ExtractSignatures.INHITS_LONG, str(inhits),
            ExtractSignatures.EXHITS_LONG, str(exhits),
            ExtractSignatures.GAP_LONG, str(gap),
            ExtractSignatures.SIZE_LONG, str(size),
            ExtractSignatures.GC_LONG, str(gcContent),
            FilterSignatures.FILTER_LENGTH_LONG, str(filterLength),
            FilterSignatures.FILTER_PERCENT_LONG, str(filterPercent),
            CountKMers.PARALLEL_LONG, str(parallelization),
            ExtractSignatures.INCLUSION_LONG, str(inclusion),
            ExtractSignatures.EXCLUSION_LONG, str(exclusion),
            ExtractSignatures.REFERENCE_SIZE_LONG, str(referenceSize),
            OUTPUT_LONG, str(outputDirectoryLocation)]

        with self.assertRaises(RuntimeError):
            main()

    def test_invalid_filter_percent(self):

        kmer = 5
        rate = 0.01
        inhits = 1
        exhits = 2
        gap = 5
        size = 7
        gcContent = 0.05
        filterLength = 0.5
        filterPercent = -0.01
        parallelization = 0
        inclusion = getPath("tests/data/simple.fasta")
        exclusion = getPath("tests/data/alternative.fasta")
        referenceSize = 12
        outputDirectoryLocation = getPath("tests/output/temp.dir")

        sys.argv[1:] = [
            CountKMers.KMER_LONG, str(kmer),
            ExtractSignatures.RATE_LONG, str(rate),
            ExtractSignatures.INHITS_LONG, str(inhits),
            ExtractSignatures.EXHITS_LONG, str(exhits),
            ExtractSignatures.GAP_LONG, str(gap),
            ExtractSignatures.SIZE_LONG, str(size),
            ExtractSignatures.GC_LONG, str(gcContent),
            FilterSignatures.FILTER_LENGTH_LONG, str(filterLength),
            FilterSignatures.FILTER_PERCENT_LONG, str(filterPercent),
            CountKMers.PARALLEL_LONG, str(parallelization),
            ExtractSignatures.INCLUSION_LONG, str(inclusion),
            ExtractSignatures.EXCLUSION_LONG, str(exclusion),
            ExtractSignatures.REFERENCE_SIZE_LONG, str(referenceSize),
            OUTPUT_LONG, str(outputDirectoryLocation)]

        with self.assertRaises(RuntimeError):
            main()

        filterPercent = 1.01

        sys.argv[1:] = [
            CountKMers.KMER_LONG, str(kmer),
            ExtractSignatures.RATE_LONG, str(rate),
            ExtractSignatures.INHITS_LONG, str(inhits),
            ExtractSignatures.EXHITS_LONG, str(exhits),
            ExtractSignatures.GAP_LONG, str(gap),
            ExtractSignatures.SIZE_LONG, str(size),
            ExtractSignatures.GC_LONG, str(gcContent),
            FilterSignatures.FILTER_LENGTH_LONG, str(filterLength),
            FilterSignatures.FILTER_PERCENT_LONG, str(filterPercent),
            CountKMers.PARALLEL_LONG, str(parallelization),
            ExtractSignatures.INCLUSION_LONG, str(inclusion),
            ExtractSignatures.EXCLUSION_LONG, str(exclusion),
            ExtractSignatures.REFERENCE_SIZE_LONG, str(referenceSize),
            OUTPUT_LONG, str(outputDirectoryLocation)]

        with self.assertRaises(RuntimeError):
            main()

    def test_invalid_parallelization(self):

        kmer = 5
        rate = 0.01
        inhits = 1
        exhits = 2
        gap = 5
        size = 5
        gcContent = 0.5
        filterLength = 0.5
        filterPercent = 0.5
        parallelization = -1
        inclusion = getPath("tests/data/simple.fasta")
        exclusion = getPath("tests/data/alternative.fasta")
        referenceSize = 12
        outputDirectoryLocation = getPath("tests/output/temp.dir")

        sys.argv[1:] = [
            CountKMers.KMER_LONG, str(kmer),
            ExtractSignatures.RATE_LONG, str(rate),
            ExtractSignatures.INHITS_LONG, str(inhits),
            ExtractSignatures.EXHITS_LONG, str(exhits),
            ExtractSignatures.GAP_LONG, str(gap),
            ExtractSignatures.SIZE_LONG, str(size),
            ExtractSignatures.GC_LONG, str(gcContent),
            FilterSignatures.FILTER_LENGTH_LONG, str(filterLength),
            FilterSignatures.FILTER_PERCENT_LONG, str(filterPercent),
            CountKMers.PARALLEL_LONG, str(parallelization),
            ExtractSignatures.INCLUSION_LONG, str(inclusion),
            ExtractSignatures.EXCLUSION_LONG, str(exclusion),
            ExtractSignatures.REFERENCE_SIZE_LONG, str(referenceSize),
            OUTPUT_LONG, str(outputDirectoryLocation)]

        with self.assertRaises(RuntimeError):
            main()

    # designed to catch possible parallel race-condition failure with many files
    def test_multiple_files(self):

        kmer = 5
        rate = 0.01
        inhits = 1
        exhits = 2
        gap = 5
        size = 4
        gcContent = 0.5
        filterLength = 0.5
        filterPercent = 0.5
        parallelization = 0
        inclusion = getPath("tests/data/neptune/large_inclusion")
        exclusion = getPath("tests/data/neptune/large_exclusion")
        referenceSize = 12
        outputDirectoryLocation = getPath("tests/output/neptune/temp.dir")
        sortedDirectoryLocation = os.path.join(outputDirectoryLocation, "sorted")

        sys.argv[1:] = [
            CountKMers.KMER_LONG, str(kmer),
            ExtractSignatures.RATE_LONG, str(rate),
            ExtractSignatures.INHITS_LONG, str(inhits),
            ExtractSignatures.EXHITS_LONG, str(exhits),
            ExtractSignatures.GAP_LONG, str(gap),
            ExtractSignatures.SIZE_LONG, str(size),
            ExtractSignatures.GC_LONG, str(gcContent),
            FilterSignatures.FILTER_LENGTH_LONG, str(filterLength),
            FilterSignatures.FILTER_PERCENT_LONG, str(filterPercent),
            CountKMers.PARALLEL_LONG, str(parallelization),
            ExtractSignatures.INCLUSION_LONG, str(inclusion),
            ExtractSignatures.EXCLUSION_LONG, str(exclusion),
            ExtractSignatures.REFERENCE_SIZE_LONG, str(referenceSize),
            OUTPUT_LONG, str(outputDirectoryLocation)]

        main()

        with open (os.path.join(sortedDirectoryLocation, "in0.fasta"), "r") as myfile:

            result = myfile.readline()
            expected = ">0 score=0.8889 in=1.0000 ex=-0.1111 len=99 ref=inclusion0 pos=999"
            self.assertEquals(result, expected)

            result = myfile.readline()
            expected = "CGGTTTCTTCATATATAACCCCGTCGGCGCTTCAGAAAACAGGGATGTATAGAATCTCTGCGTCAGAACGGCATCTAAAATCAAAACGGTATTGATGAC\n"
            self.assertEquals(result, expected)

        shutil.rmtree(outputDirectoryLocation)

if __name__ == '__main__':
    
    unittest.main()
