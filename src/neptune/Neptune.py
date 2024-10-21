#!/usr/bin/env python

"""
# =============================================================================

Copyright Government of Canada 2015-2017

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

__version__ = '1.2.5-python3'

import time

import os
import argparse
import sys
import shutil

import neptune.Execution as Execution
import neptune.Utility as Utility
import neptune.CountKMers as CountKMers
import neptune.ExtractSignatures as ExtractSignatures
import neptune.FilterSignatures as FilterSignatures

"""
# =============================================================================

GLOBALS

# =============================================================================
"""

PROGRAM_DESCRIPTION = "Neptune identifies signatures using an exact k-mer \
    matching strategy. Neptune locates sequence that is sufficiently present \
    in many inclusion targets and sufficiently absent from exclusion targets."

PROGRAM_USAGE = "%(prog)s -i INCLUSION [INCLUSION ...] -e EXCLUSION \n\t \
    [EXCLUSION ...] -o OUTPUT"

# FILE NAMES #

KMERS = "kmers"
INCLUSION = "inclusion"
EXCLUSION = "exclusion"

AGGREGATE = "aggregate.kmers"
RECEIPT = "receipt.txt"
CANDIDATES = "candidates"
FILTERED = "filtered"
SORTED = "sorted"
DATABASE = "database"
CONSOLIDATED = "consolidated"
LOG = "log"

# ARGUMENTS #

# Note: Some of the command line arguments are drawn from other scripts:

# CountKMers.py
# AggregateKMers.py
# ExtractSignatures.py
# FilterSignatures.py
# ConsolidateSignatures.py

LONG = "--"
SHORT = "-"

# REQUIRED ARGUMENTS #

OUTPUT = "output"
OUTPUT_LONG = LONG + OUTPUT
OUTPUT_SHORT = SHORT + "o"
OUTPUT_HELP = "The directory to place all output."

# OPTIONAL ARGUMENTS #

# Version number
VERSION = "version"
VERSION_LONG = LONG + VERSION
VERSION_SHORT = SHORT + "V"

# Number of parallel threads to use
PARALLELIZATION = "parallelization"
PARALLELIZATION_LONG = LONG + PARALLELIZATION
PARALLELIZATION_SHORT = SHORT + "p"
PARALLELIZATION_HELP = "The number of processes to run simultaneously."

"""
# =============================================================================

COUNT K-MERS
------------


PURPOSE
-------

This function prepares and runs a CountKMers.py job for every inclusion and
exclusion file.


INPUT
-----

[EXECUTION] [execution]
    The Execution object containing all of the current execution's parameters.


RETURN
------

[((FILE LOCATION) ITERATOR, (FILE LOCATION) ITERATOR) TUPLE]
[(inclusionKMerLocations, exclusionKMerLocations)]
    The returned tuple will contain two lists of the locations of the
    written inclusion and exclusion k-mers.

# =============================================================================
"""
def countKMers(execution):

    if not os.path.exists(execution.inclusionOutputDirectory):
        os.makedirs(execution.inclusionOutputDirectory)

    if not os.path.exists(execution.exclusionOutputDirectory):
        os.makedirs(execution.exclusionOutputDirectory)

    jobs = []
    inclusionKMerLocations = []
    exclusionKMerLocations = []

    # INCLUSION
    for inclusionLocation in execution.inclusionLocations:

        baseName = os.path.basename(inclusionLocation)

        outputLocation = os.path.abspath(
            os.path.join(
                execution.inclusionOutputDirectory, baseName + ".kmers"))

        inclusionKMerLocations.append(outputLocation)

        job = execution.jobManager.createCountJob(
            inclusionLocation, outputLocation,
            execution.k, execution.organization)
        jobs.append(job)

    # EXCLUSION
    for exclusionLocation in execution.exclusionLocations:

        baseName = os.path.basename(exclusionLocation)

        outputLocation = os.path.abspath(
            os.path.join(
                execution.exclusionOutputDirectory, baseName + ".kmers"))

        exclusionKMerLocations.append(outputLocation)

        job = execution.jobManager.createCountJob(
            exclusionLocation, outputLocation,
            execution.k, execution.organization)
        jobs.append(job)

    execution.jobManager.runJobs(jobs)

    return (inclusionKMerLocations, exclusionKMerLocations)


"""
# =============================================================================

AGGREGATE K-MERS
----------------


PURPOSE
-------

This function prepares and runs a AggregateKMers.py job which, in turn,
aggregates multuple prepared k-mer files into a single k-mers file.


INPUT
-----

[EXECUTION] [execution]
    The Execution object containing all of the current execution's parameters.

[(FILE LOCATION) ITERATOR] [inclusionKMerLocations]
    A list of inclusion k-mer file output locations.

[(FILE LOCATION) ITERATOR] [exclusionKMerLocations]
    A list of exclusion k-mer file output locations.


RETURN
------

[NONE]

# =============================================================================
"""
def aggregateKMers(execution, inclusionKMerLocations, exclusionKMerLocations):

    if execution.organization:

        aggregateMultipleFiles(
            execution,
            inclusionKMerLocations, exclusionKMerLocations)

    else:

        aggregateSingleFiles(
            execution, inclusionKMerLocations, exclusionKMerLocations)

    shutil.rmtree(execution.kmersOutputDirectory)


"""
# =============================================================================

AGGREGATE MULTIPLE FILES
------------------------


PURPOSE
-------

Aggregates k-mer count files produced from genomes, when there were multiple
k-mer files produced per genome.


INPUT
-----

[EXECUTION] [execution]
    The Execution object containing all of the current execution's parameters.

[(FILE LOCATION) ITERATOR] [inclusionKMerLocations]
    A list of inclusion k-mer file output locations.

[(FILE LOCATION) ITERATOR] [exclusionKMerLocations]
    A list of exclusion k-mer file output locations.


RETURN
------

[NONE]


POST
----

The k-mer count files are aggregated into a single k-mer count file. The
execution of the script will be halted until the job has finished.

# =============================================================================
"""
def aggregateMultipleFiles(execution, inclusionLocations, exclusionLocations):

    jobs = []
    outputLocations = []

    for tag in Utility.getAggregationTags(execution.organization):

        outputLocation = execution.aggregateLocation + "." + tag
        outputLocations.append(outputLocation)

        job = execution.jobManager.createAggregateJob(
            inclusionLocations, exclusionLocations,
            outputLocation, tag)
        jobs.append(job)

    execution.jobManager.runJobs(jobs)

    aggregateFile = open(execution.aggregateLocation, "w")

    for location in outputLocations:

        tempfile = open(location, "r")
        aggregateFile.write(tempfile.read())
        tempfile.close()

        os.remove(location)

    aggregateFile.close()


"""
# =============================================================================

AGGREGATE SINGLE FILES
----------------------


PURPOSE
-------

Aggregates k-mer count files produced from genomes, when there was only one
k-mer file produced per genome.


INPUT
-----

[EXECUTION] [execution]
    The Execution object containing all of the current execution's parameters.

[(FILE LOCATION) ITERATOR] [inclusionKMerLocations]
    A list of inclusion k-mer file output locations.

[(FILE LOCATION) ITERATOR] [exclusionKMerLocations]
    A list of exclusion k-mer file output locations.


RETURN
------

[NONE]


POST
----

The k-mer count files are aggregated into a single k-mer count file. The
execution of the script will be halted until the job is finished.

# =============================================================================
"""
def aggregateSingleFiles(
        execution, inclusionKMerLocations, exclusionKMerLocations):

    job = execution.jobManager.createAggregateJob(
        inclusionKMerLocations, exclusionKMerLocations,
        execution.aggregateLocation, None)

    execution.jobManager.runJobs([job])


"""
# =============================================================================

EXTRACT SIGNATURES
------------------


PURPOSE
-------

Prepares and runs an extract k-mers. This, in turn, extracts signatures from
a reference genome.


INPUT
-----

[EXECUTION] [execution]
    The Execution object containing all of the current execution's parameters.


RETURN
------

[(FILE LOCATION) ITERATOR] [outputLocations]
    The file locations of the extracted signatures, corresponding to either
    the references specified by the user, or, if none are specified, all the
    inclusion genomes.

POST
----

A job is submitted that extracts signatures from a reference genome
using information from aggregated k-mer information from inclusion and
exclusion genomes. The execution of the script will be halted until the job
has finished.

# =============================================================================
"""
def extractSignatures(execution):

    jobs = []
    outputLocations = []

    if execution.reference:
        references = execution.reference

    else:
        references = execution.inclusionLocations

    for reference in references:

        baseName = os.path.basename(reference)
        outputLocation = os.path.abspath(
            os.path.join(execution.candidatesDirectoryLocation, baseName))
        outputLocations.append(outputLocation)

        job = execution.jobManager.createExtractJob(
            reference, execution.referenceSize, execution.rate,
            execution.inclusionLocations, execution.inhits,
            execution.exclusionLocations, execution.exhits, execution.gap,
            execution.size, execution.gcContent, execution.confidence,
            execution.aggregateLocation, outputLocation)

        jobs.append(job)

    execution.jobManager.runJobs(jobs)

    return outputLocations


"""
# =============================================================================

FILTER SIGNATURES
-----------------


PURPOSE
-------

Filters the candidate signatures using the exclusion genomes.


INPUT
-----

[EXECUTION] [execution]
    The Execution object containing all of the current execution's parameters.

[(FILE LOCATION) ITERATOR] [candidateLocations]
    The location of candidate signatures.


RETURN
------

[NONE]


POST
----

A file of filtered candidates and sorted candidates will be produced. The
execution of the script will be halted until the job has finished.

# =============================================================================
"""
def filterSignatures(execution, candidateLocations):

    # MAKE DATABASES
    INCLUSION_NAME = "INCLUSION"
    EXCLUSION_NAME = "EXCLUSION"

    inclusionDatabaseLocation = os.path.abspath(
        os.path.join(execution.databaseDirectoryLocation, INCLUSION_NAME))

    inclusionAggregatedLocation = os.path.abspath(
        os.path.join(
            execution.databaseDirectoryLocation, "inclusionAggregated.fasta"))

    exclusionDatabaseLocation = os.path.abspath(
        os.path.join(execution.databaseDirectoryLocation, EXCLUSION_NAME))

    exclusionAggregatedLocation = os.path.abspath(
        os.path.join(
            execution.databaseDirectoryLocation, "exclusionAggregated.fasta"))

    inclusionDatabaseJob = execution.jobManager.createDatabaseJob(
        execution.inclusionLocations, inclusionAggregatedLocation,
        inclusionDatabaseLocation)

    exclusionDatabaseJob = execution.jobManager.createDatabaseJob(
        execution.exclusionLocations, exclusionAggregatedLocation,
        exclusionDatabaseLocation)

    execution.jobManager.runJobs([inclusionDatabaseJob, exclusionDatabaseJob])

    # Delete Aggregated Files:
    if os.path.exists(inclusionAggregatedLocation):
        os.remove(inclusionAggregatedLocation)

    if os.path.exists(exclusionAggregatedLocation):
        os.remove(exclusionAggregatedLocation)

    jobs = []
    sortedLocations = []

    # Filtering
    for candidateLocation in candidateLocations:

        baseName = os.path.basename(candidateLocation)
        filteredLocation = os.path.abspath(
            os.path.join(execution.filteredDirectoryLocation, baseName))
        sortedLocation = os.path.abspath(
            os.path.join(execution.sortedDirectoryLocation, baseName))
        sortedLocations.append(sortedLocation)

        job = execution.jobManager.createFilterJob(
            inclusionDatabaseLocation, exclusionDatabaseLocation,
            execution.inclusionLocations, execution.exclusionLocations,
            candidateLocation, filteredLocation, sortedLocation,
            execution.filterLength, execution.filterPercent,
            execution.seedSize)

        jobs.append(job)

    execution.jobManager.runJobs(jobs)

    shutil.rmtree(execution.databaseDirectoryLocation)

    return sortedLocations


"""
# =============================================================================

CONSOLIDATE SIGNATURES
----------------------


PURPOSE
-------

Consolidates signatures from several Neptune signature files into a single
representative Neptune signature file, determined by the signature score
and sequence similarity of all the contributing signatures.


INPUT
-----

[EXECUTION] [execution]
    The Execution object containing all of the current execution's parameters.

[(FILE LOCATION) LIST] [sortedLocations]
    The file locations of candidate signatures.


RETURN
------

[NONE]


POST
----

The produced signatures will be consolidated into a single signature file.

# =============================================================================
"""
def consolidateSignatures(execution, sortedLocations):

    job = execution.jobManager.createConsolidateJob(
        sortedLocations, execution.seedSize,
        execution.consolidatedDirectoryLocation)

    execution.jobManager.runJobs([job])


"""
# =============================================================================

EXECUTE

# =============================================================================
"""
def execute(execution):

    # --- K-MER COUNTING ---
    print("k-mer Counting...")
    start = time.perf_counter()
    inclusionKMerLocations, exclusionKMerLocations = countKMers(execution)
    end = time.perf_counter()
    print(str(end - start) + " seconds\n")

    # --- K-MER AGGREGATION ---
    print("k-mer Aggregation...")
    start = time.perf_counter()
    aggregateKMers(execution, inclusionKMerLocations, exclusionKMerLocations)
    end = time.perf_counter()
    print(str(end - start) + " seconds\n")

    # --- SIGNATURE EXTRACTION ---
    print("Signature Extraction...")
    start = time.perf_counter()
    candidateLocations = extractSignatures(execution)
    end = time.perf_counter()
    print(str(end - start) + " seconds\n")

    # --- SIGNATURE FILTERING ---
    print("Signature Filtering...")
    start = time.perf_counter()
    sortedLocations = filterSignatures(execution, candidateLocations)
    end = time.perf_counter()
    print(str(end - start) + " seconds\n")

    # Are all the signature files empty?
    if (all((os.stat(location).st_size == 0)
            for location in sortedLocations)):

        # Yes -- they are all empty.
        print("NOTICE: No signatures were identified.\n")

    else:

        # --- CONSOLIDATE SIGNATURES ---
        print("Consolidate Signatures...")
        start = time.perf_counter()
        consolidateSignatures(execution, sortedLocations)
        end = time.perf_counter()
        print(str(end - start) + " seconds\n")

    execution.produceReceipt()

    print("Complete!")


"""
# =============================================================================

EXECUTE PARALLEL

# =============================================================================
"""
def executeParallel(parameters):

    import neptune.JobManagerParallel as JobManagerParallel

    outputDirectoryLocation = os.path.abspath(parameters[OUTPUT])
    logDirectoryLocation = os.path.abspath(
        os.path.join(outputDirectoryLocation, LOG))

    parallel = parameters[PARALLELIZATION]

    jobManager = JobManagerParallel.JobManagerParallel(
        outputDirectoryLocation, logDirectoryLocation, parallel)

    execution = Execution.Execution(jobManager, parameters)
    execute(execution)


"""
# =============================================================================

PARSE

# =============================================================================
"""
def parse(parameters):

    executeParallel(parameters)


"""
# =============================================================================

MAIN

# =============================================================================
"""
def main():

    # --- PARSER --- #
    parser = argparse.ArgumentParser(
        description=PROGRAM_DESCRIPTION,
        usage=PROGRAM_USAGE)

    # --- VERSION --- #
    parser.add_argument(
        VERSION_SHORT,
        VERSION_LONG,
        action='version',
        version='%(prog)s ' + str(__version__))

    # --- REQUIRED --- #
    required = parser.add_argument_group("REQUIRED")

    required.add_argument(
        ExtractSignatures.INCLUSION_SHORT,
        ExtractSignatures.INCLUSION_LONG,
        dest=ExtractSignatures.INCLUSION,
        help=ExtractSignatures.INCLUSION_HELP,
        type=str, required=True, nargs='+')

    required.add_argument(
        ExtractSignatures.EXCLUSION_SHORT,
        ExtractSignatures.EXCLUSION_LONG,
        dest=ExtractSignatures.EXCLUSION,
        help=ExtractSignatures.EXCLUSION_HELP,
        type=str, required=True, nargs='+')

    required.add_argument(
        OUTPUT_SHORT,
        OUTPUT_LONG,
        dest=OUTPUT,
        help=OUTPUT_HELP,
        type=str, required=True)

    # --- KMERS --- #
    kmers = parser.add_argument_group("KMERS")

    kmers.add_argument(
        CountKMers.KMER_SHORT,
        CountKMers.KMER_LONG,
        dest=CountKMers.KMER,
        help=CountKMers.KMER_HELP,
        type=int, required=False)

    kmers.add_argument(
        CountKMers.ORGANIZATION_LONG,
        dest=CountKMers.ORGANIZATION,
        help=CountKMers.ORGANIZATION_HELP,
        type=int, default=3)

    # --- FILTERING --- #
    filtering = parser.add_argument_group("FILTERING")

    filtering.add_argument(
        FilterSignatures.FILTER_PERCENT_LONG,
        dest=FilterSignatures.FILTER_PERCENT,
        help=FilterSignatures.FILTER_PERCENT_HELP,
        type=float, required=False)

    filtering.add_argument(
        FilterSignatures.FILTER_LENGTH_LONG,
        dest=FilterSignatures.FILTER_LENGTH,
        help=FilterSignatures.FILTER_LENGTH_HELP,
        type=float, required=False)

    filtering.add_argument(
        FilterSignatures.SEED_SIZE_LONG,
        dest=FilterSignatures.SEED_SIZE,
        help=FilterSignatures.SEED_SIZE_HELP,
        type=int, required=False)

    # --- EXTRACTION --- #
    extraction = parser.add_argument_group("EXTRACTION")

    extraction.add_argument(
        ExtractSignatures.REFERENCE_SHORT,
        ExtractSignatures.REFERENCE_LONG,
        dest=ExtractSignatures.REFERENCE,
        help=ExtractSignatures.REFERENCE_HELP,
        type=str, required=False, nargs='+')

    extraction.add_argument(
        ExtractSignatures.REFERENCE_SIZE_LONG,
        dest=ExtractSignatures.REFERENCE_SIZE,
        help=ExtractSignatures.REFERENCE_SIZE_HELP,
        type=int, required=False)

    extraction.add_argument(
        ExtractSignatures.RATE_LONG,
        dest=ExtractSignatures.RATE,
        help=ExtractSignatures.RATE_HELP,
        type=float, required=False)

    extraction.add_argument(
        ExtractSignatures.INHITS_LONG,
        dest=ExtractSignatures.INHITS,
        help=ExtractSignatures.INHITS_HELP,
        type=int, required=False)

    extraction.add_argument(
        ExtractSignatures.EXHITS_LONG,
        dest=ExtractSignatures.EXHITS,
        help=ExtractSignatures.EXHITS_HELP,
        type=int, required=False)

    extraction.add_argument(
        ExtractSignatures.GAP_LONG,
        dest=ExtractSignatures.GAP,
        help=ExtractSignatures.GAP_HELP,
        type=int, required=False)

    extraction.add_argument(
        ExtractSignatures.SIZE_LONG,
        dest=ExtractSignatures.SIZE,
        help=ExtractSignatures.SIZE_HELP,
        type=int, required=False)

    extraction.add_argument(
        ExtractSignatures.GC_LONG,
        dest=ExtractSignatures.GC_CONTENT,
        help=ExtractSignatures.GC_HELP,
        type=float, required=False)

    extraction.add_argument(
        ExtractSignatures.CONFIDENCE_LONG,
        dest=ExtractSignatures.CONFIDENCE,
        help=ExtractSignatures.CONFIDENCE_HELP,
        type=float, required=False)

    # --- PARALLELIZATION --- #
    parallelization = parser.add_argument_group("PARALLELIZATION")

    parallelization.add_argument(
        PARALLELIZATION_SHORT,
        PARALLELIZATION_LONG,
        dest=PARALLELIZATION,
        help=PARALLELIZATION_HELP,
        type=int, default=8)

    args = parser.parse_args()
    parameters = vars(args)

    print("Neptune v" + str(__version__) + "\n")
    parse(parameters)


"""
# =============================================================================
# =============================================================================
"""
if __name__ == '__main__':

    main()
