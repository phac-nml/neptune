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

"""
# =============================================================================

This file contains the JobManagerParallel class. JobManagerParallel is
responsible for managing the creation and execution of Python multiprocessing
parallel jobs. This class is designed for Neptune execution on a single CPU or
node.

# =============================================================================
"""

import multiprocessing

import JobManager
import CountKMers
import AggregateKMers
import ExtractSignatures
import FilterSignatures
import ConsolidateSignatures
import Database

"""
# =============================================================================

JOB MANAGER

# =============================================================================
"""

class JobManagerParallel(JobManager.JobManager):

    pool = multiprocessing.Pool(processes=8)

    """
    # =========================================================================

    CONSTRUCTOR

    INPUT:
        [FILE LOCATION] [outputDirectoryLocation] - The directory location to
            write Neptune output.
        [FILE LOCATION] [logDirectoryLocation] - The directory location to
            write output logs and error logs.

    # =========================================================================
    """
    def __init__(
            self, outputDirectoryLocation, logDirectoryLocation):

        # JobManager Parent Constructor
        JobManager.JobManager.__init__(
            self, outputDirectoryLocation, logDirectoryLocation)

    """
    # =========================================================================

    RUN JOBS

    PURPOSE:
        Runs all the Neptune jobs provided to the function. The jobs are
        synchronized and execution problems are reported if possible.

    INPUT:
        [JOB] [jobs] - The jobs to run in parallel.

    RETURN:
        [NONE]

    POST:
        Job submission message and errors are printed to standard output.

    # =========================================================================
    """
    def runJobs(self, jobs):

        print "Submitted " + str(len(jobs)) + " jobs."
        self.synchronize(jobs)

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
    def synchronize(self, jobs):

        for job in jobs:
            job.wait()

    """
    # =========================================================================

    CREATE COUNT JOB

    PURPOSE:
        Creates a CountKMers job.

    INPUT:
        [FILE LOCATION] [inputLocation] - The location of the input file.
        [FILE LOCATION] [outputLocation] - The location of the output file.
        [1 <= INT] [k] - The size of the k-mers.
        [0 <= INT] [parallelization] - The degree of parallelization.

    RETURN:
        [JOB] [job] - A CountKMers job that may be passed to RunJobs(...).

    # =========================================================================
    """
    def createCountJob(
            self, inputLocation, outputLocation, k, parallelization):

        parameters = {}

        parameters[CountKMers.INPUT] = inputLocation
        parameters[CountKMers.OUTPUT] = outputLocation
        parameters[CountKMers.KMER] = k
        parameters[CountKMers.PARALLEL] = parallelization

        job = self.pool.apply_async(
            CountKMers.parse,
            args=(parameters,))

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
        [STRING -- OPTIONAL] [tag] - The parallelization tag; used to generate
            appropriate file names from the inclusion and exclusion iterators.

    RETURN:
        [JOB] [job] - An AggregateKMers job that may be passed to RunJobs(...).

    # =========================================================================
    """
    def createAggregateJob(
            self, inclusionLocations, exclusionLocations,
            outputLocation, tag):

        parameters = {}

        inclusion = []
        exclusion = []

        # INCLUSION
        if tag:
            inclusion += (item + "." + tag for item in inclusionLocations)
        else:
            inclusion += inclusionLocations

        parameters[AggregateKMers.INCLUSION] = inclusion

        # EXCLUSION
        if tag:
            exclusion += (item + "." + tag for item in exclusionLocations)
        else:
            exclusion += exclusionLocations

        parameters[AggregateKMers.EXCLUSION] = exclusion

        # OUTPUT
        parameters[AggregateKMers.OUTPUT] = outputLocation

        # DELETE
        parameters[AggregateKMers.DELETE] = True

        job = self.pool.apply_async(
            AggregateKMers.parse,
            args=(parameters,))

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

        parameters = {}

        # REFERENCE
        parameters[ExtractSignatures.REFERENCE] = referenceLocation

        # REFERENCE SIZE
        parameters[ExtractSignatures.REFERENCE_SIZE] = referenceSize \
            if referenceSize else None

        # RATE
        parameters[ExtractSignatures.RATE] = rate \
            if rate else None

        # INCLUSION
        parameters[ExtractSignatures.INCLUSION] = inclusion \
            if inclusion else None

        # INHITS
        parameters[ExtractSignatures.INHITS] = inhits \
            if inhits else None

        # EXCLUSION
        parameters[ExtractSignatures.EXCLUSION] = exclusion \
            if exclusion else None

        # EXHITS
        parameters[ExtractSignatures.EXHITS] = exhits \
            if exhits else None

        # GAP
        parameters[ExtractSignatures.GAP] = gap \
            if gap else None

        # SIZE
        parameters[ExtractSignatures.SIZE] = size \
            if size else None

        # GC
        parameters[ExtractSignatures.GC_CONTENT] = GC \
            if GC else None

        # CONFIDENCE
        parameters[ExtractSignatures.CONFIDENCE] = confidence \
            if confidence else None

        # AGGREGATED KMERS
        parameters[ExtractSignatures.KMERS] = aggregateLocation

        # OUTPUT
        parameters[ExtractSignatures.OUTPUT] = outputLocation

        job = self.pool.apply_async(
            ExtractSignatures.parse,
            args=(parameters,))

        return job

    """
    # =========================================================================

    CREATE DATABASE JOB

    PURPOSE:
        Creates a BuildDatabase job.

    INPUT:
        [(FILE LOCATION) ITERATOR] [inputLocations] - The input locations of
            the entries (FASTA) in the database.
        [FILE LOCATION] [outputLocation] - The output location of the database.

    RETURN:
        [JOB] [job] - An BuildDatabase job that may be passed to RunJobs(...).

    # =========================================================================
    """
    def createDatabaseJob(
            self, inputLocations, aggregatedLocation, outputLocation):

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

        job = self.pool.apply_async(
            Database.createDatabaseJob,
            args=(aggregatedLocation, outputLocation,))

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

        parameters = {}

        # INCLUSION DATABASE
        parameters[FilterSignatures.INCLUSION_DATABASE] = \
            inclusionDatabaseLocation

        # EXCLUSION DATABASE
        parameters[FilterSignatures.EXCLUSION_DATABASE] = \
            exclusionDatabaseLocation

        # INCLUSION
        parameters[FilterSignatures.INCLUSION] = inclusion \
            if inclusion else None

        # EXCLUSION
        parameters[FilterSignatures.EXCLUSION] = exclusion \
            if exclusion else None

        # INPUT
        parameters[FilterSignatures.INPUT] = inputLocation

        # FILTERED OUTPUT
        parameters[FilterSignatures.FILTERED_OUTPUT] = filteredOutputLocation

        # SORTED OUTPUT
        parameters[FilterSignatures.SORTED_OUTPUT] = sortedOutputLocation

        # FILTER LENGTH
        parameters[FilterSignatures.FILTER_LENGTH] = filterLength \
            if filterLength else None

        # FILTER PERCENT
        parameters[FilterSignatures.FILTER_PERCENT] = filterPercent \
            if filterPercent else None

        # SEED SIZE
        parameters[FilterSignatures.SEED_SIZE] = seedSize \
            if seedSize else None

        job = self.pool.apply_async(
            FilterSignatures.parse,
            args=(parameters,))

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

        parameters = {}

        # SIGNATURE LOCATIONS
        parameters[ConsolidateSignatures.SIGNATURES] = signatureLocations \
            if signatureLocations else None

        # SEED SIZE
        parameters[ConsolidateSignatures.SEED_SIZE] = seedSize \
            if seedSize else None

        # OUTPUT DIRECTORY LOCATION
        parameters[ConsolidateSignatures.OUTPUT] = outputDirectoryLocation \
            if outputDirectoryLocation else None

        job = self.pool.apply_async(
            ConsolidateSignatures.parse,
            args=(parameters,))

        return job
