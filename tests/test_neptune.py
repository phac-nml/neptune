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
import shutil

from tests.TestingUtility import *

from neptune.Neptune import *

import unittest

"""
# =============================================================================

DEFAULT ARGS

# =============================================================================
"""
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
            parameters[ExtractSignatures.INCLUSION] = getPath("tests/data/neptune/simple.fasta")
            parameters[ExtractSignatures.EXCLUSION] = getPath("tests/data/neptune/alternative.fasta")
            parameters[ExtractSignatures.REFERENCE] = ["tests/data/neptune/simple.fasta"]
            parameters[ExtractSignatures.REFERENCE_SIZE] = 12

            parameters[FilterSignatures.FILTER_LENGTH] = 0.5
            parameters[FilterSignatures.FILTER_PERCENT] = 0.5
            parameters[FilterSignatures.SEED_SIZE] = 9

            parameters[OUTPUT] = getPath("tests/output/neptune/temp.dir")

            self.parameters = parameters

"""
# =============================================================================

SET COMMAND LINE ARGUMENTS

PURPOSE:
    Sets the Python command line arguments for Neptune using the provided
    [parameters].

INPUT:
    [(STRING -> STRING) DICTIONARY] [parameters] - The Neptune command line
    parameters. These parameters are in a dictionary mapping the argument to
    it's parameter value.

RETURN:
    [NONE]

POST:
    The [sys.argv] command line arguments will be set.

# =============================================================================
"""
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
            CountKMers.ORGANIZATION_LONG, str(parameters[CountKMers.ORGANIZATION]),
            ExtractSignatures.INCLUSION_LONG, str(parameters[ExtractSignatures.INCLUSION]),
            ExtractSignatures.EXCLUSION_LONG, str(parameters[ExtractSignatures.EXCLUSION]),
            ExtractSignatures.REFERENCE_SIZE_LONG, str(parameters[ExtractSignatures.REFERENCE_SIZE]),
            OUTPUT_LONG, str(parameters[OUTPUT])]

"""
# =============================================================================

TEST MAIN

# =============================================================================
"""
class TestMain(unittest.TestCase):

    """
    # =========================================================================

    test_simple

    PURPOSE:
        Tests a simple execution of Neptune from the main method.

    INPUT:
        [DefaultArgs()]

    EXPECT:
        Signature: "ACGTACGTACGTACGTACGTACGTACGT\n"

    # =========================================================================
    """
    def test_simple(self):

        parameters = DefaultArgs().parameters

        setCommandLineArguments(parameters)
        main()

        sortedOutputDirectory = os.path.join(parameters[OUTPUT], "sorted")
        with open (os.path.join(sortedOutputDirectory, "simple.fasta"), "r") as myfile:

            myfile.readline()
            result = myfile.readline()
            expected = "ACGTACGTACGTACGTACGTACGTACGT\n"
            self.assertEqual(result, expected)

        shutil.rmtree(parameters[OUTPUT])

    """
    # =========================================================================

    test_simple_input_directories

    PURPOSE:
        Tests a simple execution of Neptune with inclusion and exclusion input
        directories instead of single files.

    INPUT:
        [DefaultArgs()]
        parameters[ExtractSignatures.INCLUSION] = getPath("tests/data/neptune/simple_inclusion")
        parameters[ExtractSignatures.EXCLUSION] = getPath("tests/data/neptune/simple_exclusion")

    EXPECTED:
        Signature: "ACGTACGTACGTACGTACGTACGTACGT\n"

    # =========================================================================
    """
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
            self.assertEqual(result, expected)

        shutil.rmtree(parameters[OUTPUT])

    """
    # =========================================================================

    Test a Genome with Only Ns (Simple)

    PURPOSE:
        Tests an execution of Neptune from the main method where an input
        contains a genome of entirely Ns. This only specifies mandatory
        parameters.

    INPUT:
        parameters[ExtractSignatures.INCLUSION] = getPath("tests/data/neptune/n.fasta")
        parameters[ExtractSignatures.EXCLUSION] = ["tests/data/neptune/alternative.fasta"]
        parameters[OUTPUT] = getPath("tests/output/neptune/temp.dir")

    EXPECT:
        RuntimeError when estimating the GC content of the n.fasta file.

    # =========================================================================
    """
    def test_n_genome_simple(self):

        parameters = {}
        parameters[ExtractSignatures.INCLUSION] = getPath("tests/data/neptune/n.fasta")
        parameters[ExtractSignatures.EXCLUSION] = ["tests/data/neptune/alternative.fasta"]
        parameters[OUTPUT] = getPath("tests/output/neptune/temp.dir")

        sys.argv[1:] = [
            ExtractSignatures.INCLUSION_LONG, str(parameters[ExtractSignatures.INCLUSION]),
            ExtractSignatures.EXCLUSION_LONG, str(parameters[ExtractSignatures.EXCLUSION]),
            OUTPUT_LONG, str(parameters[OUTPUT])]

        with self.assertRaises(RuntimeError):
            main()

    """
    # =========================================================================

    Test a Genome with Only Ns (Specifying GC and k)

    PURPOSE:
        Tests an execution of Neptune from the main method where an input
        contains a genome of entirely Ns. This function specifies the GC
        content and k-mer size, so as to circumvent the RuntimeError that would
        be thrown when trying to calculate these parameters for an input
        containing only Ns.

    INPUT:
        parameters[ExtractSignatures.INCLUSION] = getPath("tests/data/neptune/n.fasta")
        parameters[ExtractSignatures.EXCLUSION] = ["tests/data/neptune/alternative.fasta"]
        parameters[OUTPUT] = getPath("tests/output/neptune/temp.dir")

    EXPECT:
        RuntimeError when estimating the GC content of the n.fasta file.

    # =========================================================================
    """
    def test_n_genome_simple_specify_gc_k(self):

        parameters = {}
        parameters[ExtractSignatures.INCLUSION] = getPath("tests/data/neptune/n.fasta")
        parameters[ExtractSignatures.EXCLUSION] = ["tests/data/neptune/alternative.fasta"]
        parameters[OUTPUT] = getPath("tests/output/neptune/temp.dir")
        parameters[ExtractSignatures.GC_CONTENT] = 0.5
        parameters[CountKMers.KMER] = 5

        sys.argv[1:] = [
            ExtractSignatures.INCLUSION_LONG, str(parameters[ExtractSignatures.INCLUSION]),
            ExtractSignatures.EXCLUSION_LONG, str(parameters[ExtractSignatures.EXCLUSION]),
            OUTPUT_LONG, str(parameters[OUTPUT]),
            ExtractSignatures.GC_LONG, str(parameters[ExtractSignatures.GC_CONTENT]),
            CountKMers.KMER_LONG, str(parameters[CountKMers.KMER])]

        print(sys.argv[1:])

        with self.assertRaises(RuntimeError):
            main()

    """
    # =========================================================================

    test_organization

    PURPOSE:
        Tests a CountKMers multiple file organization execution when running
        Neptune.

    INPUT:
        [DefaultArgs()]
        parameters[CountKMers.ORGANIZATION] = 1

    EXPECTED:
        Signature: "ACGTACGTACGTACGTACGTACGTACGT\n"

    # =========================================================================
    """
    def test_organization(self):

        parameters = DefaultArgs().parameters
        parameters[CountKMers.ORGANIZATION] = 1

        setCommandLineArguments(parameters)
        main()

        sortedOutputDirectory = os.path.join(parameters[OUTPUT], "sorted")
        with open (os.path.join(sortedOutputDirectory, "simple.fasta"), "r") as myfile:

            myfile.readline()
            result = myfile.readline()
            expected = "ACGTACGTACGTACGTACGTACGTACGT\n"
            self.assertEqual(result, expected)

        shutil.rmtree(parameters[OUTPUT])

    """
    # =========================================================================

    test_specify_reference

    PURPOSE:
        Tests a Neptune execution when the reference is specified.

    INPUT:
        [DefaultArgs()]
        parameters[ExtractSignatures.REFERENCE] = getPath("tests/data/neptune/simple.fasta")

    EXPECTED:
        Signature: "ACGTACGTACGTACGTACGTACGTACGT\n"

    # =========================================================================
    """
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
            self.assertEqual(result, expected)

        shutil.rmtree(parameters[OUTPUT])

    """
    # =========================================================================

    test_invalid_k

    PURPOSE:
        Tests a Neptune execution when the k-mer size is invalid.

    INPUT:
        [DefaultArgs()]
        parameters[CountKMers.KMER] = -1

    EXPECTED:
        [RuntimeError]

    # =========================================================================
    """
    def test_invalid_k(self):

        parameters = DefaultArgs().parameters
        parameters[CountKMers.KMER] = -1

        setCommandLineArguments(parameters)

        with self.assertRaises(RuntimeError):
            main()

    """
    # =========================================================================

    test_invalid_rate

    PURPOSE:
        Tests a Neptune execution when the SNV rate is invalid.

    INPUT:
        [DefaultArgs()]
        parameters[ExtractSignatures.RATE] = -0.01

        AND

        [DefaultArgs()]
        parameters[ExtractSignatures.RATE] = 1.01

    EXPECTED:
        [RuntimeError]

    # =========================================================================
    """
    def test_invalid_rate(self):

        parameters = DefaultArgs().parameters
        parameters[ExtractSignatures.RATE] = -0.01

        setCommandLineArguments(parameters)

        with self.assertRaises(RuntimeError):
            main()

        parameters = DefaultArgs().parameters
        parameters[ExtractSignatures.RATE] = 1.01

        setCommandLineArguments(parameters)

        with self.assertRaises(RuntimeError):
            main()

    """
    # =========================================================================

    test_invalid_inhits

    PURPOSE:
        Tests a Neptune execution when the minimum inclusion hits is invalid.

    INPUT:
        [DefaultArgs()]
        parameters[ExtractSignatures.INHITS] = -1

    EXPECTED:
        [RuntimeError]

    # =========================================================================
    """
    def test_invalid_inhits(self):

        parameters = DefaultArgs().parameters
        parameters[ExtractSignatures.INHITS] = -1

        setCommandLineArguments(parameters)

        with self.assertRaises(RuntimeError):
            main()

    """
    # =========================================================================

    test_invalid_exhits

    PURPOSE:
        Tests a Neptune execution when the minimum exclusion hits is invalid.

    INPUT:
        [DefaultArgs()]
        parameters[ExtractSignatures.EXHITS] = -1

    EXPECTED:
        [RuntimeError]

    # =========================================================================
    """
    def test_invalid_exhits(self):

        parameters = DefaultArgs().parameters
        parameters[ExtractSignatures.EXHITS] = -1

        setCommandLineArguments(parameters)

        with self.assertRaises(RuntimeError):
            main()

    """
    # =========================================================================

    test_invalid_gap

    PURPOSE:
        Tests a Neptune execution when the signature gap size is invalid.

    INPUT:
        [DefaultArgs()]
        parameters[ExtractSignatures.GAP] = -1

    EXPECTED:
        [RuntimeError]

    # =========================================================================
    """
    def test_invalid_gap(self):

        parameters = DefaultArgs().parameters
        parameters[ExtractSignatures.GAP] = -1

        setCommandLineArguments(parameters)

        with self.assertRaises(RuntimeError):
            main()

    """
    # =========================================================================

    test_invalid_size

    PURPOSE:
        Tests a Neptune execution when the minimum signature size is invalid.

    INPUT:
        [DefaultArgs()]
        parameters[ExtractSignatures.SIZE] = -1

    EXPECTED:
        [RuntimeError]

    # =========================================================================
    """
    def test_invalid_size(self):

        parameters = DefaultArgs().parameters
        parameters[ExtractSignatures.SIZE] = -1

        setCommandLineArguments(parameters)

        with self.assertRaises(RuntimeError):
            main()

    """
    # =========================================================================

    test_invalid_gc

    PURPOSE:
        Tests a Neptune execution when the GC content is invalid.

    INPUT:
        [DefaultArgs()]
        parameters[ExtractSignatures.GC_CONTENT] = -0.01

        AND

        [DefaultArgs()]
        parameters[ExtractSignatures.GC_CONTENT] = 1.01

    EXPECTED:
        [RuntimeError]

    # =========================================================================
    """
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

    """
    # =========================================================================

    test_invalid_filter_length

    PURPOSE:
        Tests a Neptune execution when the alignment filter length is invalid.

    INPUT:
        [DefaultArgs()]
        parameters[FilterSignatures.FILTER_LENGTH] = -0.01

        AND

        [DefaultArgs()]
        parameters[FilterSignatures.FILTER_LENGTH] = 1.01

    EXPECTED:
        [RuntimeError]

    # =========================================================================
    """
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

    """
    # =========================================================================

    test_invalid_filter_percent

    PURPOSE:
        Tests a Neptune execution when the alignment filter percent is invalid.

    INPUT:
        [DefaultArgs()]
        parameters[FilterSignatures.FILTER_PERCENT] = -0.01

        AND

        [DefaultArgs()]
        parameters[FilterSignatures.FILTER_PERCENT] = 1.01

    EXPECTED:
        [RuntimeError]

    # =========================================================================
    """
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

    """
    # =========================================================================

    test_invalid_organization

    PURPOSE:
        Tests a Neptune execution when the CountKMers organization value is
        invalid.

    INPUT:
        [DefaultArgs()]
        parameters[CountKMers.ORGANIZATION] = -1

    EXPECTED:
        [RuntimeError]

    # =========================================================================
    """
    def test_invalid_organization(self):

        parameters = DefaultArgs().parameters
        parameters[CountKMers.ORGANIZATION] = -1

        setCommandLineArguments(parameters)

        with self.assertRaises(RuntimeError):
            main()

    """
    # =========================================================================

    test_multiple_files

    PURPOSE:
        Tests an execution of Neptune for a possible parallel race-condition
        when run with several input files.

    INPUT:
        inclusion = "tests/data/neptune/large_inclusion" (10 files)
        exclusion = "tests/data/neptune/large_exclusion" (10 files)
        output = "tests/output/neptune/temp.dir"

    EXPECT:
        signature:

        ">0 score=0.8889 in=1.0000 ex=-0.1111 len=99 ref=inclusion0 pos=999\n"
        "CGGTTTCTTCATATATAACCCCGTCGGCGCTTCAGAAAACAGGGATGTATAGAATCTCTGCGTCAGAACGGCATCTAAAATCAAAACGGTATTGATGAC\n"

    # =========================================================================
    """
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
            expected = ">0 score=0.8889 in=1.0000 ex=0.1111 len=99 ref=inclusion0 pos=999\n"
            self.assertEqual(result, expected)

            result = myfile.readline()
            expected = "CGGTTTCTTCATATATAACCCCGTCGGCGCTTCAGAAAACAGGGATGTATAGAATCTCTGCGTCAGAACGGCATCTAAAATCAAAACGGTATTGATGAC\n"
            self.assertEqual(result, expected)

        shutil.rmtree(outputDirectoryLocation)
    
    """
    # =========================================================================

    test_ambiguous_signature

    PURPOSE:
        Tests that Neptune does not crash when attempting to consolidate a
        candidate signature with many ambiguous sequence characters.

        Specifically, the issue is that in previous versions of Neptune
        (<=1.2.5), Neptune makes an assumption during the consolidate
        signatures step that every signature will at least align with itself.
        Normally, this is reasonable because a sequence should align perfectly
        with itself, but when such a sequence contains ambiguous sequence
        characters (N, etc.), BLAST will not (normally) produce an alignment.
        This was causing a crash because Neptune attempted a dictionary look-up
        for signatures alignments, but there were no alignments associated with
        the candidate signature ID. This would look like: "KeyError: '0.0'"

    INPUT:
        inclusion = "tests/data/ambiguous/inclusion.fasta"
        exclusion = "tests/data/ambiguous/exclusion.fasta"
        output = "tests/output/neptune/temp.dir"

    EXPECT:
        >0.1 score=1.0000 in=1.0000 ex=0.0000 len=640 ref=inclusion1 pos=3497
        CGCGGGCGATATTTTCACAGCCATTTCAGGAGTTCAGCCATGAACGCTTATTACATTCAGGATCGTCTT..
        >0.2 score=1.0000 in=1.0000 ex=0.0000 len=102 ref=inclusion1 pos=5205
        CATGGCGAGTTTTGCGAGATGGTGCCGGAGTTCATCGAAAAAATGGACGAGGCACTGCTGAAATTGGTT..
        >0.0 score=0.2039 in=0.2039 ex=0.0000 len=103 ref=inclusion1 pos=99
        TAGTCTCCAGGATTCCCGGGGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN..

    # =========================================================================
    """
    def test_ambiguous_signature(self):

        inclusion = getPath("tests/data/ambiguous/inclusion.fasta")
        exclusion = getPath("tests/data/ambiguous/exclusion.fasta")
        outputDirectoryLocation = getPath("tests/output/neptune/temp.dir")
        consolidatedDirectoryLocation = os.path.join(outputDirectoryLocation, "consolidated")

        sys.argv[1:] = [
            ExtractSignatures.INCLUSION_LONG, str(inclusion),
            ExtractSignatures.EXCLUSION_LONG, str(exclusion),
            OUTPUT_LONG, str(outputDirectoryLocation)]

        main()

        # Normally, opening a consolidated file might cause problems because of
        # race conditions, but with one inclusion and exclusion file, there
        # will be no problem.
        with open (os.path.join(consolidatedDirectoryLocation, "consolidated.fasta"), "r") as myfile:

            contents = myfile.read()

            expected = ">0.1 score=1.0000 in=1.0000 ex=0.0000 len=640 ref=inclusion1 pos=3497"
            self.assertTrue(expected in contents)
            expected = "CGCGGGCGATATTTTCACAGCCATTTCAGGAGTTCAGCCATGAACGCTTATTACATTCAGGATCGTCTTGAGGCTCAGAGCTGGGAGCGTCACTACCAGCAG" # trimming this one
            self.assertTrue(expected in contents)

            expected = ">0.2 score=1.0000 in=1.0000 ex=0.0000 len=102 ref=inclusion1 pos=5205"
            self.assertTrue(expected in contents)
            expected = "CATGGCGAGTTTTGCGAGATGGTGCCGGAGTTCATCGAAAAAATGGACGAGGCACTGCTGAAATTGGTTTTGTATTTGGGGAGCAATGGCGATGAAGCATCC"
            self.assertTrue(expected in contents)

            expected = ">0.0 score=0.2039 in=0.2039 ex=0.0000 len=103 ref=inclusion1 pos=99"
            self.assertTrue(expected in contents)
            expected = "TAGTCTCCAGGATTCCCGGGGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAGCAGTGATGTAAGAAAA"
            self.assertTrue(expected in contents)

        shutil.rmtree(outputDirectoryLocation)

if __name__ == '__main__':
    
    unittest.main()
