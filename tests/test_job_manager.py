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

import drmaa
import os
import sys
import StringIO

from TestingUtility import *
prepareSystemPath()

from neptune.JobManager import *
from neptune.JobManagerDRMAA import JobManagerDRMAA
from neptune.JobManagerParallel import JobManagerParallel

import neptune.CountKMers as CountKMers
import neptune.AggregateKMers as AggregateKMers
import neptune.ExtractSignatures as ExtractSignatures
import neptune.FilterSignatures as FilterSignatures
import neptune.Utility as Utility

import unittest

class TestDRMAAConstructor(unittest.TestCase):

    def test_no_defaults(self):

        with drmaa.Session() as session:

            outputDirectoryLocation = getPath("tests/output/manager")
            logDirectoryLocation = getPath("tests/output/manager/log")
            defaultSpecification = None

            jobManager = JobManagerDRMAA(
                outputDirectoryLocation, logDirectoryLocation, session, defaultSpecification)

            self.assertEquals(jobManager.session, session)
            self.assertEquals(jobManager.outputDirectoryLocation, outputDirectoryLocation)
            self.assertEquals(jobManager.logDirectoryLocation, logDirectoryLocation)

            self.assertEquals(jobManager.countSpecification, defaultSpecification)
            self.assertEquals(jobManager.aggregateSpecification, defaultSpecification)
            self.assertEquals(jobManager.extractSpecification, defaultSpecification)
            self.assertEquals(jobManager.databaseSpecification, defaultSpecification)
            self.assertEquals(jobManager.filterSpecification, defaultSpecification)

    def test_defaults(self):

        with drmaa.Session() as session:

            outputDirectoryLocation = getPath("tests/output/manager")
            logDirectoryLocation = getPath("tests/output/manager/log")
            defaultSpecification = "-l h_vmem=16G -pe smp 8"

            jobManager = JobManagerDRMAA(
                outputDirectoryLocation, logDirectoryLocation, session, defaultSpecification)

            self.assertEquals(jobManager.session, session)
            self.assertEquals(jobManager.outputDirectoryLocation, outputDirectoryLocation)

            self.assertEquals(jobManager.countSpecification, defaultSpecification)
            self.assertEquals(jobManager.aggregateSpecification, defaultSpecification)
            self.assertEquals(jobManager.extractSpecification, defaultSpecification)
            self.assertEquals(jobManager.databaseSpecification, defaultSpecification)
            self.assertEquals(jobManager.filterSpecification, defaultSpecification)

class TestParallelConstructor(unittest.TestCase):

    def test_simple(self):

        outputDirectoryLocation = getPath("tests/output/manager")
        logDirectoryLocation = getPath("tests/output/manager/log")

        jobManager = JobManagerParallel(outputDirectoryLocation, logDirectoryLocation)

        self.assertEquals(jobManager.outputDirectoryLocation, outputDirectoryLocation)
        self.assertEquals(jobManager.logDirectoryLocation, logDirectoryLocation)

class TestSetCount(unittest.TestCase):

    def test_simple(self):

        specification = "-l h_vmem=16G -pe smp 8"

        with drmaa.Session() as session:

            jobManager = JobManagerDRMAA(getPath("tests/output/manager"), getPath("tests/output/manager/log"), session, None)

            self.assertEquals(jobManager.countSpecification, None)
            jobManager.setCountSpecification(specification)
            self.assertEquals(jobManager.countSpecification, specification)

class TestSetAggregate(unittest.TestCase):

    def test_simple(self):

        specification = "-l h_vmem=16G -pe smp 8"

        with drmaa.Session() as session:

            jobManager = JobManagerDRMAA(getPath("tests/output/manager"), getPath("tests/output/manager/log"), session, None)

            self.assertEquals(jobManager.aggregateSpecification, None)
            jobManager.setAggregateSpecification(specification)
            self.assertEquals(jobManager.aggregateSpecification, specification)

class TestSetExtract(unittest.TestCase):

    def test_simple(self):

        specification = "-l h_vmem=16G -pe smp 8"

        with drmaa.Session() as session:

            jobManager = JobManagerDRMAA(getPath("tests/output/manager"), getPath("tests/output/manager/log"), session, None)

            self.assertEquals(jobManager.extractSpecification, None)
            jobManager.setExtractSpecification(specification)
            self.assertEquals(jobManager.extractSpecification, specification)

class TestSetDatabase(unittest.TestCase):

    def test_simple(self):

        specification = "-l h_vmem=16G -pe smp 8"

        with drmaa.Session() as session:

            jobManager = JobManagerDRMAA(getPath("tests/output/manager"), getPath("tests/output/manager/log"), session, None)

            self.assertEquals(jobManager.databaseSpecification, None)
            jobManager.setDatabaseSpecification(specification)
            self.assertEquals(jobManager.databaseSpecification, specification)

class TestSetFilter(unittest.TestCase):

    def test_simple(self):

        specification = "-l h_vmem=16G -pe smp 8"

        with drmaa.Session() as session:

            jobManager = JobManagerDRMAA(getPath("tests/output/manager"), getPath("tests/output/manager/log"), session, None)

            self.assertEquals(jobManager.filterSpecification, None)
            jobManager.setFilterSpecification(specification)
            self.assertEquals(jobManager.filterSpecification, specification)

class TestSetConsolidate(unittest.TestCase):

    def test_simple(self):

        specification = "-l h_vmem=16G -pe smp 8"

        with drmaa.Session() as session:

            jobManager = JobManagerDRMAA(getPath("tests/output/manager"), getPath("tests/output/manager/log"), session, None)

            self.assertEquals(jobManager.consolidateSpecification, None)
            jobManager.setConsolidateSpecification(specification)
            self.assertEquals(jobManager.consolidateSpecification, specification)

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

        self.assertEquals(result, expected)

        os.remove(outputLocation)

class TestCreateJob(unittest.TestCase):

    def test_simple(self):

        with drmaa.Session() as session:

            outputDirectoryLocation = getPath("tests/output/manager")
            logDirectoryLocation = getPath("tests/output/manager/log")
            defaultSpecification = "-l h_vmem=2G -pe smp 1"

            jobManager = JobManagerDRMAA(outputDirectoryLocation, logDirectoryLocation, session, defaultSpecification)

            job = jobManager.createJob()
            self.assertTrue(job)

class TestCreateCountJob(unittest.TestCase):

    def test_simple(self):

        with drmaa.Session() as session:

            outputDirectoryLocation = getPath("tests/output")
            logDirectoryLocation = getPath("tests/output/log")
            specification = "-l h_vmem=2G -pe smp 1"

            jobManager = JobManagerDRMAA(outputDirectoryLocation, logDirectoryLocation, session, None)
            jobManager.setCountSpecification(specification)

            inputLocation = "tests/data/manager/simple.fasta"
            outputLocation = getPath("tests/output/manager/temp.out")
            k = 7
            organization = 0

            job = jobManager.createCountJob(inputLocation, outputLocation, k, organization)

            args = [
			    CountKMers.INPUT_LONG, str(inputLocation), 
			    CountKMers.OUTPUT_LONG, str(outputLocation), 
			    CountKMers.KMER_LONG, str(k), 
			    CountKMers.ORGANIZATION_LONG, str(organization)]

            self.assertEquals(job.outputPath, ":" + os.path.join(logDirectoryLocation, "Neptune-CountKMers1.o"))
            self.assertEquals(job.errorPath, ":" + os.path.join(logDirectoryLocation, "Neptune-CountKMers1.e"))
            self.assertEquals(job.args[1:], args)
            self.assertEquals(job.nativeSpecification, specification)

class TestCreateAggregateJob(unittest.TestCase):

    def test_simple(self):

        with drmaa.Session() as session:

            outputDirectoryLocation = getPath("tests/output/manager")
            logDirectoryLocation = getPath("tests/output/manager/log")
            specification = "-l h_vmem=2G -pe smp 1"

            jobManager = JobManagerDRMAA(outputDirectoryLocation, logDirectoryLocation, session, None)
            jobManager.setAggregateSpecification(specification)

            inclusionLocations = ["tests/data/manager/simple.fasta", "tests/data/manager/alternative.fasta"]
            exclusionLocations = ["tests/data/manager/simple.fasta", "tests/data/manager/alternative.fasta"]
            outputLocation = getPath("tests/output/manager/temp.out")
            tag = "A"

            job = jobManager.createAggregateJob(inclusionLocations, exclusionLocations, outputLocation, tag)

            args = [
                AggregateKMers.INCLUSION_LONG, "tests/data/manager/simple.fasta" + "." + tag, "tests/data/manager/alternative.fasta" + "." + tag,
                AggregateKMers.EXCLUSION_LONG, "tests/data/manager/simple.fasta" + "." + tag, "tests/data/manager/alternative.fasta" + "." + tag,
                AggregateKMers.OUTPUT_LONG, outputLocation,
                AggregateKMers.DELETE_LONG]

            self.assertEquals(job.outputPath, ":" + os.path.join(logDirectoryLocation, "Neptune-AggregateKMers1.o"))
            self.assertEquals(job.errorPath, ":" + os.path.join(logDirectoryLocation, "Neptune-AggregateKMers1.e"))
            self.assertEquals(job.args[1:], args)
            self.assertEquals(job.nativeSpecification, specification)

    def test_no_tag(self):

        with drmaa.Session() as session:

            outputDirectoryLocation = getPath("tests/output/manager")
            logDirectoryLocation = getPath("tests/output/manager/log")
            specification = "-l h_vmem=2G -pe smp 1"

            jobManager = JobManagerDRMAA(outputDirectoryLocation, logDirectoryLocation, session, None)
            jobManager.setAggregateSpecification(specification)

            inclusionLocations = ["tests/data/manager/simple.fasta", "tests/data/manager/alternative.fasta"]
            exclusionLocations = ["tests/data/manager/simple.fasta", "tests/data/manager/alternative.fasta"]
            outputLocation = getPath("tests/output/manager/temp.out")
            tag = None

            job = jobManager.createAggregateJob(inclusionLocations, exclusionLocations, outputLocation, tag)

            args = [
                AggregateKMers.INCLUSION_LONG, "tests/data/manager/simple.fasta", "tests/data/manager/alternative.fasta",
                AggregateKMers.EXCLUSION_LONG, "tests/data/manager/simple.fasta", "tests/data/manager/alternative.fasta",
                AggregateKMers.OUTPUT_LONG, outputLocation,
                AggregateKMers.DELETE_LONG]

            self.assertEquals(job.outputPath, ":" + os.path.join(logDirectoryLocation, "Neptune-AggregateKMers1.o"))
            self.assertEquals(job.errorPath, ":" + os.path.join(logDirectoryLocation, "Neptune-AggregateKMers1.e"))
            self.assertEquals(job.args[1:], args)
            self.assertEquals(job.nativeSpecification, specification)

class TestCreateExtractJob(unittest.TestCase):

    def test_simple(self):

        with drmaa.Session() as session:

            outputDirectoryLocation = getPath("tests/output/manager")
            logDirectoryLocation = getPath("tests/output/manager/log")
            specification = "-l h_vmem=2G -pe smp 1"

            jobManager = JobManagerDRMAA(outputDirectoryLocation, logDirectoryLocation, session, None)
            jobManager.setExtractSpecification(specification)

            referenceLocation = "tests/data/manager/simple.fasta"
            referenceSize = 12
            rate = 0.01
            inclusion = ["tests/data/manager/simple.fasta", "tests/data/manager/alternative.fasta"]
            inhits = 2
            exclusion = ["tests/data/manager/simple.fasta", "tests/data/manager/alternative.fasta"]
            exhits = 2
            gap = 3
            size = 5
            GC = 0.5
            confidence = 0.95
            aggregateLocation = "tests/data/manager/simple.kmers"
            outputLocation = getPath("tests/output/manager/temp.out")

            job = jobManager.createExtractJob(referenceLocation, referenceSize, rate, inclusion, inhits, 
		        exclusion, exhits, gap, size, GC, confidence, aggregateLocation, outputLocation)

            args = [
			    ExtractSignatures.REFERENCE_LONG, str(referenceLocation), 
			    ExtractSignatures.REFERENCE_SIZE_LONG, str(referenceSize), 
			    ExtractSignatures.RATE_LONG, str(rate), 
			    ExtractSignatures.INCLUSION_LONG, "tests/data/manager/simple.fasta", "tests/data/manager/alternative.fasta",
                ExtractSignatures.INHITS_LONG, str(inhits),
                ExtractSignatures.EXCLUSION_LONG, "tests/data/manager/simple.fasta", "tests/data/manager/alternative.fasta",
                ExtractSignatures.EXHITS_LONG, str(exhits),
                ExtractSignatures.GAP_LONG, str(gap),
                ExtractSignatures.SIZE_LONG, str(size),
                ExtractSignatures.GC_LONG, str(GC),
                ExtractSignatures.CONFIDENCE_LONG, str(confidence),
                ExtractSignatures.KMERS_LONG, aggregateLocation,
                ExtractSignatures.OUTPUT_LONG, outputLocation]

            self.assertEquals(job.outputPath, ":" + os.path.join(logDirectoryLocation, "Neptune-ExtractSignatures1.o"))
            self.assertEquals(job.errorPath, ":" + os.path.join(logDirectoryLocation, "Neptune-ExtractSignatures1.e"))
            self.assertEquals(job.args[1:], args)
            self.assertEquals(job.nativeSpecification, specification)

class TestCreateDatabaseJob(unittest.TestCase):

    def test_simple(self):

        with drmaa.Session() as session:

            outputDirectoryLocation = getPath("tests/output/manager")
            logDirectoryLocation = getPath("tests/output/manager/log")
            specification = "-l h_vmem=2G -pe smp 1"

            jobManager = JobManagerDRMAA(outputDirectoryLocation, logDirectoryLocation, session, None)
            jobManager.setDatabaseSpecification(specification)

            inputLocations = ["tests/data/manager/simple.fasta", "tests/data/manager/alternative.fasta"]
            aggregatedLocation = getPath("tests/output/manager/aggregated.out")
            outputLocation = getPath("tests/output/manager/temp.out")

            job = jobManager.createDatabaseJob(inputLocations, aggregatedLocation, outputLocation)

            args = [
			    "-dbtype", "nucl",
                "-in", aggregatedLocation,
                "-title", "DATABASE",
                "-out", outputLocation]

            self.assertEquals(job.outputPath, ":" + os.path.join(logDirectoryLocation, "Neptune-CreateDatabase1.o"))
            self.assertEquals(job.errorPath, ":" + os.path.join(logDirectoryLocation, "Neptune-CreateDatabase1.e"))
            self.assertEquals(job.args, args)
            self.assertEquals(job.nativeSpecification, specification)

class TestCreateFilterJob(unittest.TestCase):

    def test_simple(self):

        with drmaa.Session() as session:

            outputDirectoryLocation = getPath("tests/output/manager")
            logDirectoryLocation = getPath("tests/output/manager/log")
            specification = "-l h_vmem=2G -pe smp 1"

            jobManager = JobManagerDRMAA(outputDirectoryLocation, logDirectoryLocation, session, None)
            jobManager.setFilterSpecification(specification)

            inclusionDatabaseLocation = "tests/data/manager/FAKE_IN_DB.FAKE"
            exclusionDatabaseLocation = "tests/data/manager/FAKE_EX_DB.FAKE"
            inclusion = ["tests/data/manager/simple.fasta", "tests/data/manager/alternative.fasta"]
            exclusion = ["tests/data/manager/simple.fasta", "tests/data/manager/alternative.fasta"]
            inputLocation = "tests/data/manager/simple.fasta"
            filteredOutputLocation = getPath("tests/output/manager/simple.filtered")
            sortedOutputLocation = getPath("tests/output/manager/simple.sorted")
            filterLength = 0.5
            filterPercent = 0.5
            seedSize = 11

            job = jobManager.createFilterJob(inclusionDatabaseLocation, exclusionDatabaseLocation, 
		        inclusion, exclusion, inputLocation, filteredOutputLocation, sortedOutputLocation,
                filterLength, filterPercent, seedSize)

            args = [
			    FilterSignatures.INCLUSION_DATABASE_LONG, str(inclusionDatabaseLocation), 
                FilterSignatures.EXCLUSION_DATABASE_LONG, str(exclusionDatabaseLocation), 
                FilterSignatures.INCLUSION_LONG, "tests/data/manager/simple.fasta", "tests/data/manager/alternative.fasta", 
                FilterSignatures.EXCLUSION_LONG, "tests/data/manager/simple.fasta", "tests/data/manager/alternative.fasta", 
                FilterSignatures.INPUT_LONG, str(inputLocation), 
                FilterSignatures.FILTERED_OUTPUT_LONG, str(filteredOutputLocation), 
                FilterSignatures.SORTED_OUTPUT_LONG, str(sortedOutputLocation), 
                FilterSignatures.FILTER_LENGTH_LONG, str(filterLength), 
                FilterSignatures.FILTER_PERCENT_LONG, str(filterPercent),
                FilterSignatures.SEED_SIZE_LONG, str(seedSize)]

            self.assertEquals(job.outputPath, ":" + os.path.join(logDirectoryLocation, "Neptune-FilterSignatures1.o"))
            self.assertEquals(job.errorPath, ":" + os.path.join(logDirectoryLocation, "Neptune-FilterSignatures1.e"))
            self.assertEquals(job.args[1:], args)
            self.assertEquals(job.nativeSpecification, specification)

if __name__ == '__main__':
    unittest.main()
