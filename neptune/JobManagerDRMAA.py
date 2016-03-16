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

"""
# =============================================================================

This file contains the JobManager class. JobManager is responsible for
managing the creation and execution of DRMAA jobs.

# =============================================================================
"""

import drmaa
import os
import sys
import inspect

import JobManager
import CountKMers
import AggregateKMers
import ExtractSignatures
import FilterSignatures
import ConsolidateSignatures

"""
# =============================================================================

JOB MANAGER

# =============================================================================
"""

class JobManagerDRMAA(JobManager.JobManager):

    """
    # =========================================================================

    CONSTRUCTOR

    INPUT:
        [FILE LOCATION] [outputDirectoryLocation] - The directory location to
            write DRMAA output.
        [FILE LOCATION] [logDirectoryLocation] - The directory location to
            write DRMAA output logs and error logs.
        [DRMAA SESSION] [session] - The DRMAA session.
        [STRING - OPTIONAL] [defaultSpecification] - Implementation specific
            command line options.

    # =========================================================================
    """
    def __init__(
            self, outputDirectoryLocation, logDirectoryLocation,
            session, defaultSpecification):

        # JobManager Parent Constructor
        JobManager.JobManager.__init__(
            self, outputDirectoryLocation, logDirectoryLocation)

        self.session = session

        if defaultSpecification:

            defaultSpecification = defaultSpecification.strip()

            self.countSpecification = defaultSpecification
            self.aggregateSpecification = defaultSpecification
            self.extractSpecification = defaultSpecification
            self.databaseSpecification = defaultSpecification
            self.filterSpecification = defaultSpecification
            self.consolidateSpecification = defaultSpecification

        else:

            self.countSpecification = None
            self.aggregateSpecification = None
            self.extractSpecification = None
            self.databaseSpecification = None
            self.filterSpecification = None
            self.consolidateSpecification = None

    """
    # =========================================================================

    SET COUNT SPECIFICATION

    # =========================================================================
    """
    def setCountSpecification(self, specification):
        if specification:
            self.countSpecification = specification.strip()

    """
    # =========================================================================

    SET AGGREGATE SPECIFICATION

    # =========================================================================
    """
    def setAggregateSpecification(self, specification):
        if specification:
            self.aggregateSpecification = specification.strip()

    """
    # =========================================================================

    SET EXTRACT SPECIFICATION

    # =========================================================================
    """
    def setExtractSpecification(self, specification):
        if specification:
            self.extractSpecification = specification.strip()

    """
    # =========================================================================

    SET DATABASE SPECIFICATION

    # =========================================================================
    """
    def setDatabaseSpecification(self, specification):
        if specification:
            self.databaseSpecification = specification.strip()

    """
    # =========================================================================

    SET FILTER SPECIFICATION

    # =========================================================================
    """
    def setFilterSpecification(self, specification):
        if specification:
            self.filterSpecification = specification.strip()

    """
    # =========================================================================

    SET CONSOLIDATE SPECIFICATION

    # =========================================================================
    """
    def setConsolidateSpecification(self, specification):
        if specification:
            self.consolidateSpecification = specification.strip()

    """
    # =========================================================================

    RUN JOBS

    PURPOSE:
        Runs all the DRMAA jobs provided to the function. The jobs are
        synchronized and execution problems are reported if possible.

    INPUT:
        [DRMAA JOB TEMPLATE ITERATOR] [jobs] - The DRMAA jobs to run in
            parallel.

    RETURN:
        [NONE]

    POST:
        Job submission message and errors are printed to standard output.

    # =========================================================================
    """
    def runJobs(self, jobs):

        jobIDs = []

        for job in jobs:

            jobID = self.session.runJob(job)
            jobIDs.append(jobID)

        print "Submitted " + str(len(jobIDs)) + " jobs."
        self.synchronize(jobIDs)

    """
    # =========================================================================

    SYNCHRONIZE

    PURPOSE:
        Synchronizes all the jobs associated with the passed IDs. Will output
        an error message if a job returned with an error-associated exit
        status.

    INPUT:
        [STRING ITERATOR] [jobIDs] - The unique IDs associated with every job.

    RETURN:
        [NONE]

    POST:
        The function will return when all jobs have completed. Any error
        messages will be written if there was an error-associated exit-status
        with a job.

    # =========================================================================
    """
    def synchronize(self, jobIDs):

        # Synchronize:
        self.session.synchronize(
            jobIDs, drmaa.Session.TIMEOUT_WAIT_FOREVER, False)

        # Check status and clean-up:
        for jobID in jobIDs:

            jobInfo = self.session.wait(
                jobID, drmaa.Session.TIMEOUT_WAIT_FOREVER)
            status = jobInfo.exitStatus

            # Problem:
            if status:

                print "ERROR: Job did not complete successfully."
                print "JOB: " + str(jobID)
                print "STATUS: " + str(status)

                exit(1)

    """
    # =========================================================================

    CREATE JOB

    PURPOSE:
        Creates and returns a generic DRMAA job.

    INPUT:
        [NONE]

    RETURN:
        [DRMAA JOB TEMPLATE] [job] - A distributable DRMAA job which may be
            run.

    # =========================================================================
    """
    def createJob(self):

        job = self.session.createJobTemplate()

        job.joinFiles = False
        job.jobName = self.NEPTUNE_JOB
        job.jobEnvironment = os.environ
        job.remoteCommand = sys.executable

        return job

    """
    # =========================================================================

    CREATE PYTHON JOB

    PURPOSE:
        Creates and returns a generic Python DRMAA job.

    INPUT:
        [NONE]

    RETURN:
        [DRMAA JOB TEMPLATE] [job] - A distributable Python DRMAA job which may
            be run.

    # =========================================================================
    """
    def createPythonJob(self):

        job = self.createJob()
        job.remoteCommand = sys.executable

        return job

    """
    # =========================================================================

    CREATE COUNT JOB

    PURPOSE:
        Creates a CountKMers job.

    INPUT:
        [FILE LOCATION] [inputLocation] - The location of the input file.
        [FILE LOCATION] [outputLocation] - The location of the output file.
        [1 <= INT] [k] - The size of the k-mers.
        [0 <= INT] [organization] - The degree of organization.

    RETURN:
        [JOB] [job] - A CountKMers job that may be passed to RunJobs(...).

    # =========================================================================
    """
    def createCountJob(
            self, inputLocation, outputLocation, k, organization):

        # JOB CREATION
        job = self.createPythonJob()

        job.jobName = self.COUNT_JOB
        ID = self.generateID()
        job.outputPath = ":" + os.path.join(
            self.logDirectoryLocation, self.COUNT_JOB + str(ID) + ".o")
        job.errorPath = ":" + os.path.join(
            self.logDirectoryLocation, self.COUNT_JOB + str(ID) + ".e")

        job.args = [
            os.path.realpath(inspect.getsourcefile(CountKMers)),
            CountKMers.INPUT_LONG, str(inputLocation),
            CountKMers.OUTPUT_LONG, str(outputLocation),
            CountKMers.KMER_LONG, str(k),
            CountKMers.ORGANIZATION_LONG, str(organization)]

        if self.countSpecification:
            job.nativeSpecification = self.countSpecification

        return job

    """
    # =========================================================================

    CREATE AGGREGATE JOB

    PURPOSE:
        Creates an AggregateKMers job.

    INPUT:
        [STRING ITERATOR] [inclusionLocations] - An iterable object of all
            inclusion file locations.
        [STRING ITERATOR] [exclusionLocations] - An iterable object of all
            exclusion file locations.
        [FILE LOCATION] [outputLocation] - The output file location.
        [STRING -- OPTIONAL] [tag] - The organization tag; used to generate
            appropriate file names from the inclusion and exclusion iterators.

    RETURN:
        [JOB] [job] - An AggregateKMers job that may be passed to RunJobs(...).

    # =========================================================================
    """
    def createAggregateJob(
            self, inclusionLocations, exclusionLocations,
            outputLocation, tag):

        # JOB CREATION
        job = self.createPythonJob()

        job.jobName = self.AGGREGATE_JOB
        ID = self.generateID()
        job.outputPath = ":" + os.path.join(
            self.logDirectoryLocation, self.AGGREGATE_JOB + str(ID) + ".o")
        job.errorPath = ":" + os.path.join(
            self.logDirectoryLocation, self.AGGREGATE_JOB + str(ID) + ".e")

        # COMMAND
        args = []
        args.append(os.path.realpath(inspect.getsourcefile(AggregateKMers)))

        # INCLUSION
        args.append(AggregateKMers.INCLUSION_LONG)

        if tag:
            args += (item + "." + tag for item in inclusionLocations)
        else:
            args += inclusionLocations

        # EXCLUSION
        args.append(AggregateKMers.EXCLUSION_LONG)

        if tag:
            args += (item + "." + tag for item in exclusionLocations)
        else:
            args += exclusionLocations

        # OUTPUT
        args.append(AggregateKMers.OUTPUT_LONG)
        args.append(outputLocation)

        # DELETE
        args.append(AggregateKMers.DELETE_LONG)

        job.args = args

        if self.aggregateSpecification:
            job.nativeSpecification = self.aggregateSpecification

        return job

    """
    # =========================================================================

    CREATE EXTRACT JOB

    PURPOSE:
        Creates an ExtractSignatures job.

    INPUT:
        [FILE LOCATION] [referenceLocation] - The location of the reference to
            extract candidates.
        [1 <= INT -- OPTIONAL] [referenceSize] - The size of the reference.
        [0 <= FLOAT <= 1 -- OPTIONAL] [rate] - The SNV rate.
        [1 <= INT -- OPTIONAL] [inclusion] - The number of inclusion genome
            files.
        [0 <= INT -- OPTIONAL] [inhits] - The minimum number of inclusion k-mer
            hits.
        [1 <= INT -- OPTIONAL] [exclusion] - The number of exclusion genome
            files.
        [0 <= INT -- OPTIONAL] [exhits] - The maximum number of exclusion k-mer
            hits.
        [1 <= INT -- OPTIONAL] [gap] - The maximum inclusion k-mer gap size.
        [1 <= INT -- OPTIONAL] [size] - The minimum size of any candidate.
        [0 <= FLOAT <= 1 -- OPTIONAL] [GC] - The GC-content of the environment.
        [0 < FLOAT < 1 -- OPTIONAL] [confidence] - The statistical confidence.
        [FILE LOCATION] [aggregateLocation] - The location of the aggregation
            file.
        [FILE LOCATION] [outputLocation] - The location of the output file.

    RETURN:
        [JOB] [job] - An ExtractSignatures job that may be passed to
            RunJobs(...).

    # =========================================================================
    """
    def createExtractJob(
            self, referenceLocation, referenceSize, rate, inclusion, inhits,
            exclusion, exhits, gap, size, GC, confidence, aggregateLocation,
            outputLocation):

        # JOB CREATION
        job = self.createPythonJob()

        job.jobName = self.EXTRACT_JOB
        ID = self.generateID()
        job.outputPath = ":" + os.path.join(
            self.logDirectoryLocation, self.EXTRACT_JOB + str(ID) + ".o")
        job.errorPath = ":" + os.path.join(
            self.logDirectoryLocation, self.EXTRACT_JOB + str(ID) + ".e")

        # COMMAND
        args = []
        args.append(os.path.realpath(inspect.getsourcefile(ExtractSignatures)))

        # REFERENCE
        args.append(ExtractSignatures.REFERENCE_LONG)
        args.append(str(referenceLocation))

        # REFERENCE SIZE
        if referenceSize:
            args.append(ExtractSignatures.REFERENCE_SIZE_LONG)
            args.append(str(referenceSize))

        # RATE
        if rate:
            args.append(ExtractSignatures.RATE_LONG)
            args.append(str(rate))

        # INCLUSION
        if inclusion:
            args.append(ExtractSignatures.INCLUSION_LONG)
            args += inclusion

        # INHITS
        if inhits:
            args.append(ExtractSignatures.INHITS_LONG)
            args.append(str(inhits))

        # EXCLUSION
        if exclusion:
            args.append(ExtractSignatures.EXCLUSION_LONG)
            args += exclusion

        # EXHITS
        if exhits:
            args.append(ExtractSignatures.EXHITS_LONG)
            args.append(str(exhits))

        # GAP
        if gap:
            args.append(ExtractSignatures.GAP_LONG)
            args.append(str(gap))

        # SIZE
        if size:
            args.append(ExtractSignatures.SIZE_LONG)
            args.append(str(size))

        # GC
        if GC:
            args.append(ExtractSignatures.GC_LONG)
            args.append(str(GC))

        # CONFIDENCE
        if confidence:
            args.append(ExtractSignatures.CONFIDENCE_LONG)
            args.append(str(confidence))

        # AGGREGATED KMERS
        args.append(ExtractSignatures.KMERS_LONG)
        args.append(str(aggregateLocation))

        # OUTPUT
        args.append(ExtractSignatures.OUTPUT_LONG)
        args.append(str(outputLocation))

        job.args = args

        if self.extractSpecification:
            job.nativeSpecification = self.extractSpecification

        return job

    """
    # =========================================================================

    CREATE DATABASE JOB

    PURPOSE:
        Creates a BuildDatabase job.

    INPUT:
        [(FILE LOCATION) ITERATOR] [inputLocations] - The input locations of
            the entries (FASTA) in the database.
        [FILE LOCATION] [aggregatedLocation] - The location to write a single
            database file corresponding to information from the input files.
        [FILE LOCATION] [outputLocation] - The output location of the database.

    RETURN:
        [JOB] [job] - An BuildDatabase job that may be passed to RunJobs(...).

    # =========================================================================
    """
    def createDatabaseJob(
            self, inputLocations, aggregatedLocation, outputLocation):

        # COMMAND LINE
        COMMAND = "makeblastdb"

        TYPE = "-dbtype"
        NUCLEOTIDE = "nucl"

        INPUT = "-in"
        INPUT_LOCATIONS = aggregatedLocation

        TITLE = "-title"
        NAME = "DATABASE"

        OUTPUT = "-out"
        OUTPUT_LOCATION = outputLocation

        aggregatedFile = open(aggregatedLocation, 'w')
        ID = 0

        for inputLocation in inputLocations:

            inputFile = open(inputLocation, 'r')

            for line in inputFile:

                if line[0] is ">":
                    aggregatedFile.write(">" + str(ID) + "\n")

                else:
                    aggregatedFile.write(line)

            inputFile.close()
            ID += 1

        aggregatedFile.close()

        # JOB
        job = self.createJob()
        job.remoteCommand = COMMAND

        job.jobName = self.DATABASE_JOB
        ID = self.generateID()
        job.outputPath = ":" + os.path.join(
            self.logDirectoryLocation, self.DATABASE_JOB + str(ID) + ".o")
        job.errorPath = ":" + os.path.join(
            self.logDirectoryLocation, self.DATABASE_JOB + str(ID) + ".e")

        args = []

        args.append(TYPE)
        args.append(NUCLEOTIDE)

        args.append(INPUT)
        args.append(INPUT_LOCATIONS)

        args.append(TITLE)
        args.append(NAME)

        args.append(OUTPUT)
        args.append(OUTPUT_LOCATION)

        job.args = args

        if self.databaseSpecification:
            job.nativeSpecification = self.databaseSpecification

        return job

    """
    # =========================================================================

    CREATE FILTER JOB

    PURPOSE:
        Creates a FilterSignatures job.

    INPUT:
        [FILE LOCATION] [inclusionDatabaseLocation] - The location of the
            inclusion database to compare signatures against.
        [FILE LOCATION] [exclusionDatabaseLocation] - The location of the
            exclusion database to compare signatures against.
        [(FILE LOCATION) LIST] [inclusion] - The list of inclusion files.
        [(FILE LOCATION) LIST] [exclusion] - The list of exclusion files.
        [FILE LOCATION] [inputLocation] - The candidate signatures to filter.
        [FILE LOCATION] [filteredOutputLocation] - The filtered output
            location.
        [FILE LOCATION] [sortedOutputLocation] - The sorted output location.
        [0 <= FLOAT <= 1] [filterLength] - The maximum percent length of an
            exclusion hit with a candidate.
        [0 <= FLOAT <= 1] [filterPercent] - The maximum percent identity of an
            exclusion hit with a candidate.
        [4 <= INT] [seedSize] - The seed size used in alignments.

    RETURN:
        [JOB] [job] - A FilterSignatures job that may be passed to
            RunJobs(...).

    # =========================================================================
    """
    def createFilterJob(
            self, inclusionDatabaseLocation, exclusionDatabaseLocation,
            inclusion, exclusion, inputLocation, filteredOutputLocation,
            sortedOutputLocation, filterLength, filterPercent, seedSize):

        # JOB CREATION
        job = self.createPythonJob()

        job.jobName = self.FILTER_JOB
        ID = self.generateID()
        job.outputPath = ":" + os.path.join(
            self.logDirectoryLocation, self.FILTER_JOB + str(ID) + ".o")
        job.errorPath = ":" + os.path.join(
            self.logDirectoryLocation, self.FILTER_JOB + str(ID) + ".e")

        # COMMAND
        args = []
        args.append(os.path.realpath(inspect.getsourcefile(FilterSignatures)))

        # INCLUSION DATABASE
        args.append(FilterSignatures.INCLUSION_DATABASE_LONG)
        args.append(str(inclusionDatabaseLocation))

        # EXCLUSION DATABASE
        args.append(FilterSignatures.EXCLUSION_DATABASE_LONG)
        args.append(str(exclusionDatabaseLocation))

        # INCLUSION
        if inclusion:
            args.append(FilterSignatures.INCLUSION_LONG)
            args += inclusion

        # EXCLUSION
        if exclusion:
            args.append(FilterSignatures.EXCLUSION_LONG)
            args += exclusion

        # INPUT
        args.append(FilterSignatures.INPUT_LONG)
        args.append(str(inputLocation))

        # FILTERED OUTPUT
        args.append(FilterSignatures.FILTERED_OUTPUT_LONG)
        args.append(str(filteredOutputLocation))

        # SORTED OUTPUT
        args.append(FilterSignatures.SORTED_OUTPUT_LONG)
        args.append(str(sortedOutputLocation))

        # FILTER LENGTH
        if filterLength:
            args.append(FilterSignatures.FILTER_LENGTH_LONG)
            args.append(str(filterLength))

        # FILTER PERCENT
        if filterPercent:
            args.append(FilterSignatures.FILTER_PERCENT_LONG)
            args.append(str(filterPercent))

        # SEED SIZE
        if seedSize:
            args.append(FilterSignatures.SEED_SIZE_LONG)
            args.append(str(seedSize))

        job.args = args

        if self.filterSpecification:
            job.nativeSpecification = self.filterSpecification

        return job

    """
    # =========================================================================

    CREATE CONSOLIDATE JOB

    PURPOSE:
        Creates a ConsolidateSignatures job.

    INPUT:
        [(FILE LOCATION) LIST] [signatureLocations] - A list of Neptune
            signature file locations corresponding to files to consolidate.
        [4 <= INT] [seedSize] - The seed size used in alignments.
        [(FILE DIRECTORY) LOCATION] [outputDirectoryLocation] - The directory
            to write the output files.

    RETURN:
        [JOB] [job] - A ConsolidateSignatures job that may be passed to
            RunJobs(...).

    # =========================================================================
    """
    def createConsolidateJob(
            self, signatureLocations, seedSize, outputDirectoryLocation):

        # JOB CREATION
        job = self.createPythonJob()

        job.jobName = self.CONSOLIDATE_JOB
        ID = self.generateID()
        job.outputPath = ":" + os.path.join(
            self.logDirectoryLocation, self.CONSOLIDATE_JOB + str(ID) + ".o")
        job.errorPath = ":" + os.path.join(
            self.logDirectoryLocation, self.CONSOLIDATE_JOB + str(ID) + ".e")

        # COMMAND
        args = []
        args.append(os.path.realpath(
            inspect.getsourcefile(ConsolidateSignatures)))

        # SIGNATURE LOCATIONS
        if signatureLocations:
            args.append(ConsolidateSignatures.SIGNATURES_LONG)
            args += signatureLocations

        # SEED SIZE
        if seedSize:
            args.append(ConsolidateSignatures.SEED_SIZE_LONG)
            args.append(str(seedSize))

        # OUTPUT DIRECTORY LOCATION
        if outputDirectoryLocation:
            args.append(ConsolidateSignatures.OUTPUT_LONG)
            args.append(str(outputDirectoryLocation))

        job.args = args

        if self.consolidateSpecification:
            job.nativeSpecification = self.consolidateSpecification

        return job
