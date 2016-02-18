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

This file contains the JobManager abstract class. JobManager is responsible for
managing the creation and execution of Neptune jobs.

# =============================================================================
"""

import abc

import drmaa
import os
import sys
import inspect

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

class JobManager:

    __metaclass__ = abc.ABCMeta

    NEPTUNE_JOB = "Neptune"
    COUNT_JOB = "Neptune-CountKMers"
    AGGREGATE_JOB = "Neptune-AggregateKMers"
    EXTRACT_JOB = "Neptune-ExtractSignatures"
    DATABASE_JOB = "Neptune-CreateDatabase"
    FILTER_JOB = "Neptune-FilterSignatures"
    CONSOLIDATE_JOB = "Neptune-ConsolidateSignatures"

    ID = 0

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

        self.outputDirectoryLocation = outputDirectoryLocation
        self.logDirectoryLocation = logDirectoryLocation

    """
    # =========================================================================

    GENERATE ID

    PURPOSE:
        Generates an ID unique to this JobManager.

    INPUT:
        [NONE]

    RETURN:
        [INT] [ID] - The generated unique ID.

    POST:
        [NONE]

    # =========================================================================
    """
    def generateID(self):

        self.ID += 1

        return int(self.ID)

    """
    # =========================================================================

    RUN JOBS

    PURPOSE:
        Runs the jobs provided to this function. The function returns after
        all jobs have completed running.

    INPUT:
        [JOB] [jobs] - The jobs to run in parallel.

    RETURN:
        [NONE]

    POST:
        Job submission message and errors are printed to standard output.

    # =========================================================================
    """
    @abc.abstractmethod
    def runJobs(self, jobs):
        return

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
    @abc.abstractmethod
    def createCountJob(
            self, inputLocation, outputLocation, k, parallelization):
        return

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
    @abc.abstractmethod
    def createAggregateJob(
            self, inclusionLocations, exclusionLocations,
            outputLocation, tag):
        return

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
    @abc.abstractmethod
    def createExtractJob(
            self, referenceLocation, referenceSize, rate, inclusion, inhits,
            exclusion, exhits, gap, size, GC, confidence, aggregateLocation,
            outputLocation):
        return

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
    @abc.abstractmethod
    def createDatabaseJob(
            self, inputLocations, aggregatedLocation, outputLocation):
        return

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
    @abc.abstractmethod
    def createFilterJob(
            self, inclusionDatabaseLocation, exclusionDatabaseLocation,
            inclusion, exclusion, inputLocation, filteredOutputLocation,
            sortedOutputLocation, filterLength, filterPercent, seedSize):
        return

    """
    # =========================================================================

    CREATE CONSOLIDATE JOB

    PURPOSE:
        Creates a ConsolidateSignatures job.

    INPUT:
        [(FILE LOCATION) LIST] [signatureLocations] - A list of Neptune
            signature file locations corresponding to files to consolidate.
        [4 <= INT] [seedSize] - The seed size used in alignments.
        [FILE DIRECTORY LOCATION] [outputDirectoryLocation] - The directory
            to write the output files.

    RETURN:
        [JOB] [job] - A ConsolidateSignatures job that may be passed to
            RunJobs(...).

    # =========================================================================
    """
    @abc.abstractmethod
    def createConsolidateJob(
            self, signatureLocations, seedSize, outputDirectoryLocation):
        return
