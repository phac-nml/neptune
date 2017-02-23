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

__version__ = '1.2.3'

import time

import os
import argparse
import sys
import shutil

import Execution
import Utility
import CountKMers
import ExtractSignatures
import FilterSignatures

"""
# =============================================================================

GLOBALS

# =============================================================================
"""

# FILE NAMES
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

# ARGUMENT NAMES
OUTPUT = "output"
DRMAA = "drmaa"
VERSION = "version"

PARALLELIZATION = "parallelization"
DEFAULT_SPECIFICATION = "default-specification"
COUNT_SPECIFICATION = "count-specification"
AGGREGATE_SPECIFICATION = "aggregate-specification"
EXTRACT_SPECIFICATION = "extract-specification"
DATABASE_SPECIFICATION = "database-specification"
FILTER_SPECIFICATION = "filter-specification"
CONSOLIDATE_SPECIFICATION = "consolidate-specification"

# ARGUMENTS
LONG = "--"

OUTPUT_LONG = LONG + OUTPUT
DRMAA_LONG = LONG + DRMAA
VERSION_LONG = LONG + VERSION
PARALLELIZATION_LONG = LONG + PARALLELIZATION

DEFAULT_SPECIFICATION_LONG = LONG + DEFAULT_SPECIFICATION
COUNT_SPECIFICATION_LONG = LONG + COUNT_SPECIFICATION
AGGREGATE_SPECIFICATION_LONG = LONG + AGGREGATE_SPECIFICATION
EXTRACT_SPECIFICATION_LONG = LONG + EXTRACT_SPECIFICATION
DATABASE_SPECIFICATION_LONG = LONG + DATABASE_SPECIFICATION
FILTER_SPECIFICATION_LONG = LONG + FILTER_SPECIFICATION
CONSOLIDATE_SPECIFICATION_LONG = LONG + CONSOLIDATE_SPECIFICATION

SHORT = "-"

OUTPUT_SHORT = SHORT + "o"
PARALLELIZATION_SHORT = SHORT + "p"
VERSION_SHORT = SHORT + "V"

"""
# =============================================================================

COUNT K-MERS

PURPOSE:
    This function prepares and runs a CountKMers job for every inclusion
    and exclusion file. These jobs are run using the DRMAA framework.

INPUT:
    [EXECUTION] [execution] - The Execution object containing all of the
        current execution's parameters.

RETURN:
    [((FILE LOCATION) ITERATOR, (FILE LOCATION) ITERATOR) TUPLE]
    [(inclusionKMerLocations, exclusionKMerLocations)] -
        The returned tuple will contain two lists of the locations of the
        written inclusion and exclusion k-mers.

POST:
    A DRMAA job is submitted for each inclusion and exclusion file. The
    execution will be halted until all the jobs are completed.

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

PURPOSE:
    This function prepares and runs a AggregateKMers.py job which, in turn,
    aggregates multuple prepared k-mer files into a single k-mers file.

INPUT:
    [EXECUTION] [execution] - The Execution object containing all of the
        current execution's parameters.
    [(FILE LOCATION) ITERATOR] [inclusionKMerLocations] - A list of inclusion
        k-mer file output locations.
    [(FILE LOCATION) ITERATOR] [exclusionKMerLocations] - A list of exclusion
        k-mer file output locations.

RETURN:
    [NONE]

POST:
    A DRMAA job is submitted that aggregates the prepared inclusion and
    exclusion k-mers. The execution of the script will be halted until the
    job has finished.

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

PURPOSE:
    Aggregates k-mer count files produced from genomes, such that
    there was multiple k-mer files produced per genome.

INPUT:
    [EXECUTION] [execution] - The Execution object containing all of the
        current execution's parameters.
    [(FILE LOCATION) ITERATOR] [inclusionKMerLocations] - A list of inclusion
        k-mer file output locations.
    [(FILE LOCATION) ITERATOR] [exclusionKMerLocations] - A list of exclusion
        k-mer file output locations.

RETURN:
    [NONE]

POST:
    The k-mer count files are aggregated into a single k-mer count file.

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

PURPOSE:
    Aggregates k-mer count files produced from genomes, such that
    there was only one k-mer file produced per genome.

INPUT:
    [EXECUTION] [execution] - The Execution object containing all of the
        current execution's parameters.
    [(FILE LOCATION) ITERATOR] [inclusionKMerLocations] - A list of inclusion
        k-mer file output locations.
    [(FILE LOCATION) ITERATOR] [exclusionKMerLocations] - A list of exclusion
        k-mer file output locations.

RETURN:
    [NONE]

POST:
    The k-mer count files are aggregated into a single k-mer count file.

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

PURPOSE:
    Prepares and runs an extract k-mers job using DRMAA. This, in turn,
    extracts signatures from a reference genome.

INPUT:
    [EXECUTION] [execution] - The Execution object containing all of the
        current execution's parameters.
POST:
    A DRMAA job is submitted that extracts signatures from a reference genome
    using information from aggregated k-mer information from inclusion and
    exclusion genomes. The execution of the script will be halted until the
    job has finished.

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

PURPOSE:
    Filters the candidate signatures using the exclusion genomes.

INPUT:
    [EXECUTION] [execution] - The Execution object containing all of the
        current execution's parameters.
    [(FILE LOCATION) ITERATOR] [candidateLocations] - The location of candidate
        signatures.

RETURN:
    [NONE]

POST:
    A file of filtered candidates and sorted candidates will be produced.

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

PURPOSE:
    Consolidates signatures from several Neptune signature files into a single
    representative Neptune signature file, determined by the signature score
    and sequence similarity of all the contributing signatures.

INPUT:
    [EXECUTION] [execution] - The Execution object containing all of the
        current execution's parameters.
    [(FILE LOCATION) LIST] [sortedLocations] - The file locations of candidate
        signatures.

RETURN:
    [NONE]

POST:
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

    print("Neptune v" + str(__version__) + "\n")

    # --- K-MER COUNTING ---
    print("k-mer Counting...")
    start = time.clock()
    inclusionKMerLocations, exclusionKMerLocations = countKMers(execution)
    end = time.clock()
    print(str(end - start) + " seconds\n")

    # --- K-MER AGGREGATION ---
    print("k-mer Aggregation...")
    start = time.clock()
    aggregateKMers(execution, inclusionKMerLocations, exclusionKMerLocations)
    end = time.clock()
    print(str(end - start) + " seconds\n")

    # --- SIGNATURE EXTRACTION ---
    print("Signature Extraction...")
    start = time.clock()
    candidateLocations = extractSignatures(execution)
    end = time.clock()
    print(str(end - start) + " seconds\n")

    # --- SIGNATURE FILTERING ---
    print("Signature Filtering...")
    start = time.clock()
    sortedLocations = filterSignatures(execution, candidateLocations)
    end = time.clock()
    print(str(end - start) + " seconds\n")

    # Are all the signature files empty?
    if (all((os.stat(location).st_size == 0)
            for location in sortedLocations)):

        # Yes -- they are all empty.
        print("NOTICE: No signatures were identified.\n")

    else:

        # --- CONSOLIDATE SIGNATURES ---
        print("Consolidate Signatures...")
        start = time.clock()
        consolidateSignatures(execution, sortedLocations)
        end = time.clock()
        print(str(end - start) + " seconds\n")

    execution.produceReceipt()

    print("Complete!")


"""
# =============================================================================

EXECUTE DRMAA

# =============================================================================
"""
def executeDRMAA(parameters):

    import drmaa
    import JobManagerDRMAA

    with drmaa.Session() as session:

        outputDirectoryLocation = os.path.abspath(parameters[OUTPUT])

        logDirectoryLocation = os.path.abspath(
            os.path.join(outputDirectoryLocation, LOG))

        jobManager = JobManagerDRMAA.JobManagerDRMAA(
            outputDirectoryLocation, logDirectoryLocation,
            session, parameters[DEFAULT_SPECIFICATION])

        execution = Execution.Execution(jobManager, parameters)
        execute(execution)


"""
# =============================================================================

EXECUTE PARALLEL

# =============================================================================
"""
def executeParallel(parameters):

    import JobManagerParallel

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

    if parameters.get(DRMAA):
        executeDRMAA(parameters)

    else:
        executeParallel(parameters)


"""
# =============================================================================

MAIN

# =============================================================================
"""
def main():

    # --- PARSER ---
    parser = argparse.ArgumentParser(
        description='Neptune identifies signatures using an exact k-mer \
        matching strategy. Neptune locates sequence that is sufficiently \
        present in many inclusion targets and sufficiently absent from \
        exclusion targets.',

        usage="%(prog)s -i INCLUSION [INCLUSION ...] -e EXCLUSION \n\t" +
        "[EXCLUSION ...] -o OUTPUT")

    # --- VERSION ---
    parser.add_argument(
        VERSION_SHORT,
        VERSION_LONG,
        action='version',
        version='%(prog)s ' + str(__version__))

    # --- REQUIRED ---
    required = parser.add_argument_group("REQUIRED")

    required.add_argument(
        ExtractSignatures.INCLUSION_SHORT,
        ExtractSignatures.INCLUSION_LONG,
        dest=ExtractSignatures.INCLUSION,
        help="FASTA inclusion genome(s) to investigate for signatures",
        type=str, required=True, nargs='+')

    required.add_argument(
        ExtractSignatures.EXCLUSION_SHORT,
        ExtractSignatures.EXCLUSION_LONG,
        dest=ExtractSignatures.EXCLUSION,
        help="FASTA exclusion genome(s) that will be a background for the \
        inclusion genome(s)",
        type=str, required=True, nargs='+')

    required.add_argument(
        OUTPUT_SHORT,
        OUTPUT_LONG,
        dest=OUTPUT,
        help="output directory",
        type=str, required=True)

    # --- KMERS ---
    kmers = parser.add_argument_group("KMERS")

    kmers.add_argument(
        CountKMers.KMER_SHORT,
        CountKMers.KMER_LONG,
        dest=CountKMers.KMER,
        help="k-mer size",
        type=int, required=False)

    kmers.add_argument(
        CountKMers.ORGANIZATION_LONG,
        dest=CountKMers.ORGANIZATION,
        help="number of base positions used in k-mer organization; \
        affects the the speed of k-mer aggregation",
        type=int, default=3)

    # --- FILTERING ---
    filtering = parser.add_argument_group("FILTERING")

    filtering.add_argument(
        FilterSignatures.FILTER_PERCENT_LONG,
        dest=FilterSignatures.FILTER_PERCENT,
        help="the maximum percent identity of an exclusion hit; removes \
        candidates that have exclusion hit matches higher than this value",
        type=float, required=False)

    filtering.add_argument(
        FilterSignatures.FILTER_LENGTH_LONG,
        dest=FilterSignatures.FILTER_LENGTH,
        help="the maximum shared fractional length of an exclusion hit \
            with a candidate; remove candidates that have exclusion hit \
            matches longer than this value",
        type=float, required=False)

    filtering.add_argument(
        FilterSignatures.SEED_SIZE_LONG,
        dest=FilterSignatures.SEED_SIZE,
        help="the seed size used during alignment",
        type=int, required=False)

    # --- EXTRACTION ---
    extraction = parser.add_argument_group("EXTRACTION")

    extraction.add_argument(
        ExtractSignatures.REFERENCE_SHORT,
        ExtractSignatures.REFERENCE_LONG,
        dest=ExtractSignatures.REFERENCE,
        help="FASTA reference(s) from which to extract signatures",
        type=str, required=False, nargs='+')

    extraction.add_argument(
        ExtractSignatures.REFERENCE_SIZE_LONG,
        dest=ExtractSignatures.REFERENCE_SIZE,
        help="estimated reference size",
        type=int, required=False)

    extraction.add_argument(
        ExtractSignatures.RATE_LONG,
        dest=ExtractSignatures.RATE,
        help="probability of homologous bases not matching",
        type=float, required=False)

    extraction.add_argument(
        ExtractSignatures.INHITS_LONG,
        dest=ExtractSignatures.INHITS,
        help="minimum inclusion hits to start or continue building candidate \
        a signature",
        type=int, required=False)

    extraction.add_argument(
        ExtractSignatures.EXHITS_LONG,
        dest=ExtractSignatures.EXHITS,
        help="minimum exclusion hits to terminate a candidate signature",
        type=int, required=False)

    extraction.add_argument(
        ExtractSignatures.GAP_LONG,
        dest=ExtractSignatures.GAP,
        help="maximum number of consecutive k-mers in a candidate \
            without an inclusion hit before terminating the candidate",
        type=int, required=False)

    extraction.add_argument(
        ExtractSignatures.SIZE_LONG,
        dest=ExtractSignatures.SIZE,
        help="minimum candidate signature size",
        type=int, required=False)

    extraction.add_argument(
        ExtractSignatures.GC_LONG,
        dest=ExtractSignatures.GC_CONTENT,
        help="the GC-content of the environment",
        type=float, required=False)

    extraction.add_argument(
        ExtractSignatures.CONFIDENCE_LONG,
        dest=ExtractSignatures.CONFIDENCE,
        help="statistical confidence level used when extending k-mer gaps, \
        determining minimum inclusion hit observations",
        type=float, required=False)

    # --- PARALLELIZATION ---
    parallelization = parser.add_argument_group("PARALLELIZATION")

    parallelization.add_argument(
        PARALLELIZATION_SHORT,
        PARALLELIZATION_LONG,
        dest=PARALLELIZATION,
        help="the maximum number of parallel worker processes to create \
            (non-DRMAA mode)",
        type=int, default=8)

    # --- DRMAA ---
    drmaa = parser.add_argument_group("DRMAA")

    drmaa.add_argument(
        DRMAA_LONG,
        dest=DRMAA,
        help="runs Neptune in DRMAA mode",
        action='store_true')

    drmaa.add_argument(
        DEFAULT_SPECIFICATION_LONG,
        dest=DEFAULT_SPECIFICATION,
        help="DRM-specific parameters for all jobs",
        type=str, required=False)

    drmaa.add_argument(
        COUNT_SPECIFICATION_LONG,
        dest=COUNT_SPECIFICATION,
        help="DRM-specific parameters for k-mer counting",
        type=str, required=False)

    drmaa.add_argument(
        AGGREGATE_SPECIFICATION_LONG,
        dest=AGGREGATE_SPECIFICATION,
        help="DRM-specific parameters for k-mer aggregation",
        type=str, required=False)

    drmaa.add_argument(
        EXTRACT_SPECIFICATION_LONG,
        dest=EXTRACT_SPECIFICATION,
        help="DRM-specific parameters for signature extraction",
        type=str, required=False)

    drmaa.add_argument(
        DATABASE_SPECIFICATION_LONG,
        dest=DATABASE_SPECIFICATION,
        help="DRM-specific parameters for database construction",
        type=str, required=False)

    drmaa.add_argument(
        FILTER_SPECIFICATION_LONG,
        dest=FILTER_SPECIFICATION,
        help="DRM-specific parameters for signature filtering",
        type=str, required=False)

    drmaa.add_argument(
        CONSOLIDATE_SPECIFICATION_LONG,
        dest=CONSOLIDATE_SPECIFICATION,
        help="DRM-specific parameters for signature filtering",
        type=str, required=False)

    # --- ArgParse Work-Around ---
    for i in range(len(sys.argv)):
        if ((sys.argv[i] == DEFAULT_SPECIFICATION_LONG or
                sys.argv[i] == COUNT_SPECIFICATION_LONG or
                sys.argv[i] == AGGREGATE_SPECIFICATION_LONG or
                sys.argv[i] == DATABASE_SPECIFICATION_LONG or
                sys.argv[i] == FILTER_SPECIFICATION_LONG or
                sys.argv[i] == CONSOLIDATE_SPECIFICATION_LONG)):

            sys.argv[i + 1] = (
                str(sys.argv[i + 1])[:0] + " " +
                str(sys.argv[i + 1])[0:])

    args = parser.parse_args()
    parameters = vars(args)
    parse(parameters)


"""
# =============================================================================
# =============================================================================
"""
if __name__ == '__main__':

    main()
