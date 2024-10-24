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

from neptune.JobManager import *
from neptune.JobManagerParallel import JobManagerParallel

import neptune.CountKMers as CountKMers
import neptune.AggregateKMers as AggregateKMers
import neptune.ExtractSignatures as ExtractSignatures
import neptune.FilterSignatures as FilterSignatures
import neptune.Utility as Utility

import unittest


class TestParallelConstructor(unittest.TestCase):

    def test_simple(self):

        outputDirectoryLocation = getPath("tests/output/manager")
        logDirectoryLocation = getPath("tests/output/manager/log")

        jobManager = JobManagerParallel(outputDirectoryLocation, logDirectoryLocation)

        self.assertEqual(jobManager.outputDirectoryLocation, outputDirectoryLocation)
        self.assertEqual(jobManager.logDirectoryLocation, logDirectoryLocation)

class TestRunJobs(unittest.TestCase):

    def test_simple(self):

        outputDirectoryLocation = getPath("tests/output/manager/output")
        logDirectoryLocation = getPath("tests/output/manager/log")
        jobManager = JobManagerParallel(outputDirectoryLocation, logDirectoryLocation)

        inputLocation = getPath("tests/data/manager/simple.fasta")
        outputLocation = getPath("tests/output/manager/temp.out")
        k = 7
        organization = 0

        job = jobManager.createCountJob(inputLocation, outputLocation, k, organization)

        jobManager.runJobs([job])

        with open (outputLocation, "r") as myfile:
            result = myfile.read()

        expected = "ACGTACG 4\nGTACGTA 2\n"

        self.assertEqual(result, expected)

        os.remove(outputLocation)


if __name__ == '__main__':
    unittest.main()
