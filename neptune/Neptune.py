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

import drmaa
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
PARALLELIZATION = "parallelization"
DEFAULT_SPECIFICATION = "defaultSpecification"
COUNT_SPECIFICATION = "countSpecification"
AGGREGATE_SPECIFICATION = "aggregateSpecification"
EXTRACT_SPECIFICATION = "extractSpecification"
DATABASE_SPECIFICATION = "databaseSpecification"
FILTER_SPECIFICATION = "filterSpecification"
CONSOLIDATE_SPECIFICATION = "consolidateSpecification"

# ARGUMENTS
LONG = "--"

OUTPUT_LONG = LONG + OUTPUT
DEFAULT_SPECIFICATION_LONG = LONG + DEFAULT_SPECIFICATION
COUNT_SPECIFICATION_LONG = LONG + COUNT_SPECIFICATION
AGGREGATE_SPECIFICATION_LONG = LONG + AGGREGATE_SPECIFICATION
EXTRACT_SPECIFICATION_LONG = LONG + EXTRACT_SPECIFICATION
DATABASE_SPECIFICATION_LONG = LONG + DATABASE_SPECIFICATION
FILTER_SPECIFICATION_LONG = LONG + FILTER_SPECIFICATION
CONSOLIDATE_SPECIFICATION_LONG = LONG + CONSOLIDATE_SPECIFICATION

SHORT = "-"

OUTPUT_SHORT = SHORT + "o"

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
    [(STRING ITERATOR, STRING ITERATOR) TUPLE]
    [(inclusionKMerLocations, exclusionKMerLocations)] -
        The returned tuple will contain two lists of the locations of the
        written inclusion and exclusion k-mers.

POST:
    A DRMAA job is submitted for each inclusion and exclusion file. The
    execution will be halted until all the jobs are completed.

# =============================================================================
"""
def countKMers(execution):

    print("CountKMers starting ...")

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
            execution.k, execution.parallelization)
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
            execution.k, execution.parallelization)
        jobs.append(job)

    execution.jobManager.runJobs(jobs)

    print("CountKMers finished!")

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
    [STRING ITERATOR] [inclusionKMerLocations] - A list of inclusion k-mer
        file output locations.
    [STRING ITERATOR] [exclusionKMerLocations] - A list of exclusion k-mer
        file output locations.

RETURN:
    [NONE]

POST:
    A DRMAA job is submitted that aggregates the prepared inclusion and
    exclusion k-mers. The execution of the script will be halted until the
    job has finished.

# =============================================================================
"""
def aggregateKMers(execution, inclusionKMerLocations, exclusionKMerLocations):

    print("AggregateKMers starting ...")

    if execution.parallelization:

        aggregateMultipleFiles(
            execution,
            inclusionKMerLocations, exclusionKMerLocations)

    else:

        aggregateSingleFiles(
            execution, inclusionKMerLocations, exclusionKMerLocations)

    shutil.rmtree(execution.kmersOutputDirectory)

    print("AggregateKMers finished!")

"""
# =============================================================================

AGGREGATE MULTIPLE FILES

PURPOSE:
    Aggregates k-mer count files produced from genomes, such that
    there was multiple k-mer files produced per genome.

INPUT:
    [EXECUTION] [execution] - The Execution object containing all of the
        current execution's parameters.
    [STRING ITERATOR] [inclusionKMerLocations] - A list of inclusion k-mer
        file output locations.
    [STRING ITERATOR] [exclusionKMerLocations] - A list of exclusion k-mer
        file output locations.

RETURN:
    [NONE]

POST:
    The k-mer count files are aggregated into a single k-mer count file.

# =============================================================================
"""
def aggregateMultipleFiles(execution, inclusionLocations, exclusionLocations):

    jobs = []
    outputLocations = []

    for tag in Utility.getAggregationTags(execution.parallelization):

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
    [STRING ITERATOR] [inclusionKMerLocations] - A list of inclusion k-mer
        file output locations.
    [STRING ITERATOR] [exclusionKMerLocations] - A list of exclusion k-mer
        file output locations.

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

    print("ExtractSignatures starting ...")

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

    print("ExtractSignatures finished!")

    return outputLocations

"""
# =============================================================================

FILTER SIGNATURES

PURPOSE:
    Filters the candidate signatures using the exclusion genomes.

INPUT:
    [EXECUTION] [execution] - The Execution object containing all of the
        current execution's parameters.
    [STRING ITERATOR] [candidateLocations] - The location of candidate
        signatures.

RETURN:
    [NONE]

POST:
    A file of filtered candidates will be produced.

# =============================================================================
"""
def filterSignatures(execution, candidateLocations):

    print("Filtering signatures ...")

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

    print("Filtering finished!")

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

    print("Consolidating signatures ...")

    job = execution.jobManager.createConsolidateJob(
        sortedLocations, execution.seedSize,
        execution.consolidatedDirectoryLocation)

    execution.jobManager.runJobs([job])

    print("Consolidating finished!")

"""
# =============================================================================

MAIN

# =============================================================================
"""
def main():

    # --- Parser ---
    parser = argparse.ArgumentParser(
        description='Neptune locates DNA signatures using an exact k-mer \
        matching strategy. Neptune locates sequence that is sufficiently \
        represented in many inclusion targets and sufficiently absent from \
        exclusion targets.')

    parser.add_argument(
        ExtractSignatures.REFERENCE_SHORT,
        ExtractSignatures.REFERENCE_LONG,
        dest=ExtractSignatures.REFERENCE,
        help="FASTA reference(s) from which to extract signatures",
        type=str, required=False, nargs='+')

    parser.add_argument(
        ExtractSignatures.REFERENCE_SIZE_SHORT,
        ExtractSignatures.REFERENCE_SIZE_LONG,
        dest=ExtractSignatures.REFERENCE_SIZE,
        help="estimated total reference size",
        type=int, required=False)

    parser.add_argument(
        CountKMers.KMER_SHORT,
        CountKMers.KMER_LONG,
        dest=CountKMers.KMER,
        help="k-mer size",
        type=int, required=False)

    parser.add_argument(
        ExtractSignatures.RATE_SHORT,
        ExtractSignatures.RATE_LONG,
        dest=ExtractSignatures.RATE,
        help="probability of homologous bases not matching",
        type=float, required=False)

    parser.add_argument(
        ExtractSignatures.INCLUSION_SHORT,
        ExtractSignatures.INCLUSION_LONG,
        dest=ExtractSignatures.INCLUSION,
        help="inclusion genomes",
        type=str, required=True, nargs='+')

    parser.add_argument(
        ExtractSignatures.INHITS_SHORT,
        ExtractSignatures.INHITS_LONG,
        dest=ExtractSignatures.INHITS,
        help="minimum inclusion hits to build candidate",
        type=int, required=False)

    parser.add_argument(
        ExtractSignatures.EXCLUSION_SHORT,
        ExtractSignatures.EXCLUSION_LONG,
        dest=ExtractSignatures.EXCLUSION,
        help="exclusion genome(s)",
        type=str, required=True, nargs='+')

    parser.add_argument(
        ExtractSignatures.EXHITS_SHORT,
        ExtractSignatures.EXHITS_LONG,
        dest=ExtractSignatures.EXHITS,
        help="minimum exclusion hits to remove candidate",
        type=int, required=False)

    parser.add_argument(
        ExtractSignatures.GAP_SHORT,
        ExtractSignatures.GAP_LONG,
        dest=ExtractSignatures.GAP,
        help="maximum number of consecutive k-mers in a candidate \
            without an inclusion hit",
        type=int, required=False)

    parser.add_argument(
        ExtractSignatures.SIZE_SHORT,
        ExtractSignatures.SIZE_LONG,
        dest=ExtractSignatures.SIZE,
        help="minimum candidate size",
        type=int, required=False)

    parser.add_argument(
        ExtractSignatures.GC_SHORT,
        ExtractSignatures.GC_LONG,
        dest=ExtractSignatures.GC_CONTENT,
        help="the GC-content of the environment",
        type=float, required=False)

    parser.add_argument(
        ExtractSignatures.CONFIDENCE_SHORT,
        ExtractSignatures.CONFIDENCE_LONG,
        dest=ExtractSignatures.CONFIDENCE,
        help="statistical confidence level",
        type=float, required=False)

    parser.add_argument(
        FilterSignatures.FILTER_PERCENT_SHORT,
        FilterSignatures.FILTER_PERCENT_LONG,
        dest=FilterSignatures.FILTER_PERCENT,
        help="the maximum percent identity of an exclusion hit; removes \
        candidates that have exclusion hit matches higher than this",
        type=float, required=False)

    parser.add_argument(
        FilterSignatures.FILTER_LENGTH_SHORT,
        FilterSignatures.FILTER_LENGTH_LONG,
        dest=FilterSignatures.FILTER_LENGTH,
        help="the maximum shared fractional length of an exclusion hit \
            with a candidate; remove candidates that have exclusion hit \
            matches longer than this",
        type=float, required=False)

    parser.add_argument(
        FilterSignatures.SEED_SIZE_SHORT,
        FilterSignatures.SEED_SIZE_LONG,
        dest=FilterSignatures.SEED_SIZE,
        help="the seed size used during alignment",
        type=int, required=False)

    parser.add_argument(
        OUTPUT_SHORT,
        OUTPUT_LONG,
        dest=OUTPUT,
        help="output directory",
        type=str, required=True)

    parser.add_argument(
        CountKMers.PARALLEL_SHORT,
        CountKMers.PARALLEL_LONG,
        dest=CountKMers.PARALLEL,
        help="number of base positions used in parallelization",
        type=int, default=0)

    parser.add_argument(
        DEFAULT_SPECIFICATION_LONG,
        dest=DEFAULT_SPECIFICATION,
        help="DRM-specific parameters for all jobs",
        type=str, required=False)

    parser.add_argument(
        COUNT_SPECIFICATION_LONG,
        dest=COUNT_SPECIFICATION,
        help="DRM-specific parameters for k-mer counting",
        type=str, required=False)

    parser.add_argument(
        AGGREGATE_SPECIFICATION_LONG,
        dest=AGGREGATE_SPECIFICATION,
        help="DRM-specific parameters for k-mer aggregation",
        type=str, required=False)

    parser.add_argument(
        EXTRACT_SPECIFICATION_LONG,
        dest=EXTRACT_SPECIFICATION,
        help="DRM-specific parameters for signature extraction",
        type=str, required=False)

    parser.add_argument(
        DATABASE_SPECIFICATION_LONG,
        dest=DATABASE_SPECIFICATION,
        help="DRM-specific parameters for database construction",
        type=str, required=False)

    parser.add_argument(
        FILTER_SPECIFICATION_LONG,
        dest=FILTER_SPECIFICATION,
        help="DRM-specific parameters for signature filtering",
        type=str, required=False)

    parser.add_argument(
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
                str(sys.argv[i + 1])[:0] + " "
                + str(sys.argv[i + 1])[0:])

    args = parser.parse_args()

    # --- Job Control ---
    with drmaa.Session() as session:

        execution = Execution.Execution(session, args)

        # --- K-MER COUNTING ---
        inclusionKMerLocations, exclusionKMerLocations = countKMers(execution)

        # --- K-MER AGGREGATION ---
        aggregateKMers(
            execution, inclusionKMerLocations, exclusionKMerLocations)

        # --- SIGNATURE EXTRACTION ---
        candidateLocations = extractSignatures(execution)

        # --- SIGNATURE FILTERING ---
        sortedLocations = filterSignatures(execution, candidateLocations)

        # Are all the signature files empty?
        if(all((os.stat(location).st_size == 0)
                for location in sortedLocations)):

            print "NOTICE: No signatures were identified."
            return

        # --- CONSOLIDATE SIGNATURES ---
        consolidateSignatures(execution, sortedLocations)

        execution.produceReceipt()

if __name__ == '__main__':

    main()
