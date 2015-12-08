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

Author: Eric Marinier
Date: 22 April 2015

# =============================================================================
"""

import drmaa
import os
import argparse
import sys
import math
from scipy.misc import comb

import Utility
import JobManager
import CountKMers
import ExtractSignatures
import FilterSignatures

"""
# =============================================================================

GLOBALS

# =============================================================================
"""

KMERS = "kmers"
INCLUSION = "inclusion"
EXCLUSION = "exclusion"

AGGREGATE = "aggregate.kmers"
RECEIPT = "receipt.txt"
CANDIDATES = "candidates"
FILTERED = "filtered"
SORTED = "sorted"
DATABASE = "database"
LOG = "log"

# NAMES
OUTPUT = "output"
PARALLELIZATION = "parallelization"
DEFAULT_SPECIFICATION = "defaultSpecification"
COUNT_SPECIFICATION = "countSpecification"
AGGREGATE_SPECIFICATION = "aggregateSpecification"
EXTRACT_SPECIFICATION = "extractSpecification"
DATABASE_SPECIFICATION = "databaseSpecification"
FILTER_SPECIFICATION = "filterSpecification"

# ARGUMENTS
LONG = "--"

OUTPUT_LONG = LONG + OUTPUT
DEFAULT_SPECIFICATION_LONG = LONG + DEFAULT_SPECIFICATION
COUNT_SPECIFICATION_LONG = LONG + COUNT_SPECIFICATION
AGGREGATE_SPECIFICATION_LONG = LONG + AGGREGATE_SPECIFICATION
EXTRACT_SPECIFICATION_LONG = LONG + EXTRACT_SPECIFICATION
DATABASE_SPECIFICATION_LONG = LONG + DATABASE_SPECIFICATION
FILTER_SPECIFICATION_LONG = LONG + FILTER_SPECIFICATION

SHORT = "-"

OUTPUT_SHORT = SHORT + "o"

EXPECTED_HITS_THRESHOLD = 0.05

"""
# =============================================================================

EXECUTION

# =============================================================================
"""
class Execution():

    def __init__(self, session, args):

        # -- rate --
        # 0.0 <= q <= 1.0
        if (args.rate is not None and
                (float(args.rate) < 0.0 or
                    float(args.rate) > 1.0)):
            raise RuntimeError("The rate is out of range.")

        self.rate = args.rate

        # -- minimum inclusion hits --
        # 1 <= inhits
        if (args.inhits is not None and
                (int(args.inhits) < 1)):
            raise RuntimeError("The inclusion hits is out of range.")

        self.inhits = args.inhits

        # -- minimum exclusion hits --
        # 1 <= exhits
        if (args.exhits is not None and
                (int(args.exhits) < 1)):
            raise RuntimeError("The exclusion hits is out of range.")

        self.exhits = args.exhits

        # -- maximum gap size --
        # 1 <= gap
        if (args.gap is not None and
                (int(args.gap) < 1)):
            raise RuntimeError("The gap size is out of range.")

        self.gap = args.gap

        # -- minimum signature size --
        # 1 <= size
        if (args.size is not None and
                (int(args.size) < 1)):
            raise RuntimeError("The signature size is out of range.")

        self.size = args.size

        # -- GC-content --
        # 0.0 <= gc <= 1.0
        if (args.gcContent is not None and
                (float(args.gcContent) < 0.0 or
                    float(args.gcContent) > 1.0)):
            raise RuntimeError("The GC-content is out of range.")

        self.gcContent = args.gcContent

        # -- statistical confidence --
        # 0.0 < confidence < 1.0
        if (args.confidence is not None and
                (float(args.confidence) <= 0.0 or
                    float(args.confidence) >= 1.0)):
            raise RuntimeError("The statistical confidence is out of range.")

        self.confidence = args.confidence

        # -- filter length --
        # 0.0 <= filterLength <= 1.0
        if (args.filterLength is not None and
                (float(args.filterLength) < 0.0 or
                    float(args.filterLength) > 1.0)):
            raise RuntimeError("The filter length is out of range.")

        self.filterLength = args.filterLength

        # -- filter percent --
        # 0.0 <= filterPercent <= 1.0
        if (args.filterPercent is not None and
                (float(args.filterPercent) < 0.0 or
                    float(args.filterPercent) > 1.0)):
            raise RuntimeError("The filter percent is out of range.")

        self.filterPercent = args.filterPercent

        # -- seed size --
        # 4 <= seedSize
        if (args.seedSize is not None and
                (int(args.seedSize) < 4)):
            raise RuntimeError("The seed size is out of range.")

        self.seedSize = args.seedSize

        # -- parallelization --
        # 1 <= parallelization
        if (args.parallelization is not None and
                (int(args.parallelization) < 0)):
            raise RuntimeError("The parallelization is out of range.")

        self.parallelization = args.parallelization

        # -- inclusion locations --
        # inclusion exists
        if args.inclusion is None:
            raise RuntimeError("Inclusion sequence(s) are missing.")

        self.inclusionLocations = []
        self.expandInput(args.inclusion, self.inclusionLocations)

        if len(args.inclusion) is 0:
            raise RuntimeError("Inclusion sequence(s) are missing.")

        # -- exclusion locations --
        # exclusion exists
        if args.exclusion is None:
            raise RuntimeError("Exclusion sequence(s) are missing.")

        self.exclusionLocations = []
        self.expandInput(args.exclusion, self.exclusionLocations)

        if len(args.exclusion) is 0:
            raise RuntimeError("exclusion sequence(s) are missing.")

        # -- reference locations --
        self.reference = args.reference

        # -- reference size --
        self.referenceSize = args.referenceSize

        # -- output locations --
        # output exists
        if args.output is None:
            raise RuntimeError("The output directory is missing.")

        # -- k-mer --
        # 1 <= k
        if (args.kmer is not None and
                (int(args.kmer) < 1)):
            raise RuntimeError("The k-mer size is out of range.")

        elif args.kmer is not None:
            self.k = int(args.kmer)

        else:
            self.estimateKMerSize()

        self.outputDirectoryLocation = os.path.abspath(args.output)

        if not os.path.exists(self.outputDirectoryLocation):
            os.makedirs(self.outputDirectoryLocation)

        self.candidatesDirectoryLocation = os.path.abspath(
            os.path.join(self.outputDirectoryLocation, CANDIDATES))
        if not os.path.exists(self.candidatesDirectoryLocation):
            os.makedirs(self.candidatesDirectoryLocation)

        self.filteredDirectoryLocation = os.path.abspath(
            os.path.join(self.outputDirectoryLocation, FILTERED))
        if not os.path.exists(self.filteredDirectoryLocation):
            os.makedirs(self.filteredDirectoryLocation)

        self.sortedDirectoryLocation = os.path.abspath(
            os.path.join(self.outputDirectoryLocation, SORTED))
        if not os.path.exists(self.sortedDirectoryLocation):
            os.makedirs(self.sortedDirectoryLocation)

        self.databaseDirectoryLocation = os.path.abspath(
            os.path.join(self.outputDirectoryLocation, DATABASE))
        if not os.path.exists(self.databaseDirectoryLocation):
            os.makedirs(self.databaseDirectoryLocation)

        self.kmersOutputDirectory = os.path.abspath(
            os.path.join(self.outputDirectoryLocation, KMERS))
        self.inclusionOutputDirectory = os.path.abspath(
            os.path.join(self.kmersOutputDirectory, INCLUSION))
        self.exclusionOutputDirectory = os.path.abspath(
            os.path.join(self.kmersOutputDirectory, EXCLUSION))

        self.aggregateLocation = os.path.abspath(
            os.path.join(self.outputDirectoryLocation, AGGREGATE))

        self.logDirectoryLocation = os.path.abspath(
            os.path.join(self.outputDirectoryLocation, LOG))
        if not os.path.exists(self.logDirectoryLocation):
            os.makedirs(self.logDirectoryLocation)

        # -- job manager --
        self.jobManager = JobManager.JobManager(
            session, self.outputDirectoryLocation,
            self.logDirectoryLocation, args.defaultSpecification)

        # -- job specifications --
        if args.countSpecification:
            self.jobManager.setCountSpecification(
                args.countSpecification)

        if args.aggregateSpecification:
            self.jobManager.setAggregateSpecification(
                args.aggregateSpecification)

        if args.extractSpecification:
            self.jobManager.setExtractSpecification(
                args.extractSpecification)

        if args.databaseSpecification:
            self.jobManager.setDatabaseSpecification(
                args.databaseSpecification)

        if args.filterSpecification:
            self.jobManager.setFilterSpecification(
                args.filterSpecification)

    """
    # =========================================================================

    EXPAND INPUT

    PURPOSE:
        Expands the file input locations to include all files within all
        directories.

        The directories are not included in the final list. However, all
            non-directories located within the directories are included as
            individual files.

    INPUT:
        [STRING ITERATOR] [inputLocations] - The file input locations.
        [STRING LIST] [result] - The list to fill with the expanded input
            entries.

    RETURN:
        [NONE]

    POST:
        The passed [result] parameter will contain the expanded input.

    # =========================================================================
    """
    @staticmethod
    def expandInput(inputLocations, result):

        # Expand directories into files:
        for location in inputLocations:

            if os.path.isdir(location):

                onlyfiles = [
                    os.path.join(location, f)
                    for f in os.listdir(location)
                    if os.path.isfile(os.path.join(location, f))]

                result += onlyfiles

            else:

                result += [location]

        # Convert to absolute path:
        for i in range(len(result)):

            result[i] = os.path.abspath(result[i])

    """
    # =========================================================================

    CALCULATE EXPECTED K-MER HITS

    PURPOSE:
        Calculates the expected number of arbitrary k-mer matches with a
        genome.

        This is calculating P(k_x = k_y) * ((gs - k + 1) C (2)).

    INPUT:
        [0 <= FLOAT <= 1] [GC] - The GC content.
        [0 <= FLOAT <= 1] [GS] - The genome size.
        [1 <= INT] [K] - The k-mer size.

    RETURN:
        [FLOAT >= 0] [expected] - The expected number of arbitrary k-mer
            matches.

    # =========================================================================
    """
    @staticmethod
    def calculateExpectedKMerHits(gc, gs, k):

        # P(k_x = k_y) * ((gs - k + 1) C (2)) -- from manuscript
        a = 2.0 * math.pow((1.0 - gc) / 2.0, 2.0)
        b = 2.0 * math.pow(gc / 2.0, 2.0)
        c = math.pow(a + b, k)
        d = comb((gs - k + 1), (2), False)
        expected = c * d

        return expected

    """
    # =========================================================================

    ESTIMATE K-MER SIZE

    PURPOSE:
        Estimates the appropriate k-mer size.

    INPUT:
        [NONE]

    RETURN:
        [NONE]

    POST:
        The member variable [self.k] is assigned an appriorate value of k,
        determined by the most extreme GC-content of all references and the
        largest size of any reference.
    # =========================================================================
    """
    def estimateKMerSize(self):

        print "Estimating k-mer size ..."

        maxGenomeSize = 1
        maxGCContent = 0.5  # least extreme GC-content

        for inclusionLocation in self.inclusionLocations:

            inclusionFile = open(inclusionLocation, 'r')

            size = 0
            sumGC = 0
            sumAT = 0
            gcContent = 0

            for line in inclusionFile:

                if line[0] != ">":
                    line = line.strip()

                    size += len(line)
                    sumGC += line.count('G') + line.count('C')
                    sumAT += line.count('A') + line.count('T')

            gcContent = float(sumGC) / float(sumGC + sumAT)

            # SIZE
            if size > maxGenomeSize:
                maxGenomeSize = size

            # GC-CONTENT
            if gcContent > maxGCContent:
                maxGCContent = gcContent

            elif (1.0 - gcContent) > maxGCContent:
                maxGCContent = (1.0 - gcContent)

            inclusionFile.close()

        """
        NOTE:

        When GC-content is 1.0 (worst),
        and the length of genome is 10^80 (atoms in observable universe),
        and the threshold is 0.05,
        then the k required is 535.
        Therefore, k = 535 is the maximum we use.

        """
        for k in range(3, 535, 2):

            expected = self.calculateExpectedKMerHits(
                maxGCContent, maxGenomeSize, k)

            if expected < EXPECTED_HITS_THRESHOLD:

                self.k = k
                print "k = " + str(self.k)
                return

        # No suitable k estimated.
        raise RuntimeError("ERROR: No suitable value for k determined. \n")

    """
    # =========================================================================

    PRODUCE RECEIPT

    PURPOSE:
        Produces an execution receipt. This receipt provides information about
        the execution parameters.

    INPUT:
        [NONE]

    RETURN:
        [NONE]

    POST:
        An execution receipt is printed to standard output.

    # =========================================================================
    """
    def produceReceipt(self):

        receiptLocation = os.path.abspath(
            os.path.join(self.outputDirectoryLocation, RECEIPT))
        receiptFile = open(receiptLocation, "w")

        receiptFile.write(
            "==============================================================\n")
        receiptFile.write("RUN RECEIPT\n")
        receiptFile.write(
            "==============================================================\n")
        receiptFile.write("\n")
        receiptFile.write("-- Command Line -- \n")
        receiptFile.write("\n")

        for arg in sys.argv:
            receiptFile.write(str(arg) + " ")
        receiptFile.write("\n")

        receiptFile.write("\n")
        receiptFile.write("-- General -- \n")
        receiptFile.write("\n")
        receiptFile.write(
            "k = "
            + str(self.k) + "\n")
        receiptFile.write(
            "rate = "
            + str(self.rate) + "\n")
        receiptFile.write(
            "min inclusion = "
            + str(self.inhits) + "\n")
        receiptFile.write(
            "min exclusion = "
            + str(self.exhits) + "\n")
        receiptFile.write(
            "max gap size = "
            + str(self.gap) + "\n")
        receiptFile.write(
            "min signature size = "
            + str(self.size) + "\n")
        receiptFile.write(
            "GC-content = "
            + str(self.gcContent) + "\n")
        receiptFile.write(
            "filter length = "
            + str(self.filterLength) + "\n")
        receiptFile.write(
            "filter percent = "
            + str(self.filterPercent) + "\n")
        receiptFile.write(
            "parallelization = "
            + str(self.parallelization) + "\n")
        receiptFile.write(
            "referenceSize = "
            + str(self.referenceSize) + "\n")
        receiptFile.write("\n")

        receiptFile.write("-- Files -- \n")
        receiptFile.write("\n")
        receiptFile.write("inclusion targets = \n")

        for location in self.inclusionLocations:
            receiptFile.write("\t" + str(location) + "\n")

        receiptFile.write("exclusion targets = \n")

        for location in self.exclusionLocations:
            receiptFile.write("\t" + str(location) + "\n")

        receiptFile.write("reference = " + str(self.reference) + "\n")
        receiptFile.write(
            "output = "
            + str(self.outputDirectoryLocation) + "\n")
        receiptFile.write("\n")

        receiptFile.write("-- DRMAA -- \n")
        receiptFile.write("\n")
        receiptFile.write(
            "countSpecification = "
            + str(self.jobManager.countSpecification) + "\n")
        receiptFile.write(
            "aggregateSpecification = "
            + str(self.jobManager.aggregateSpecification) + "\n")
        receiptFile.write(
            "extractSpecification = "
            + str(self.jobManager.extractSpecification) + "\n")
        receiptFile.write(
            "databaseSpecification = "
            + str(self.jobManager.databaseSpecification) + "\n")
        receiptFile.write(
            "filterSpecification = "
            + str(self.jobManager.filterSpecification) + "\n")
        receiptFile.write("\n")

        receiptFile.close()

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

    # Filtering
    for candidateLocation in candidateLocations:

        baseName = os.path.basename(candidateLocation)
        filteredLocation = os.path.abspath(
            os.path.join(execution.filteredDirectoryLocation, baseName))
        sortedLocation = os.path.abspath(
            os.path.join(execution.sortedDirectoryLocation, baseName))

        job = execution.jobManager.createFilterJob(
            inclusionDatabaseLocation, exclusionDatabaseLocation,
            execution.inclusionLocations, execution.exclusionLocations,
            candidateLocation, filteredLocation, sortedLocation,
            execution.filterLength, execution.filterPercent,
            execution.seedSize)

        jobs.append(job)

    execution.jobManager.runJobs(jobs)

    print("Filtering finished!")

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

    os.rmdir(execution.inclusionOutputDirectory)
    os.rmdir(execution.exclusionOutputDirectory)
    os.rmdir(execution.kmersOutputDirectory)

    print("AggregateKMers finished!")

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

    # --- ArgParse Work-Around ---
    for i in range(len(sys.argv)):
        if ((sys.argv[i] == DEFAULT_SPECIFICATION_LONG or
                sys.argv[i] == COUNT_SPECIFICATION_LONG or
                sys.argv[i] == AGGREGATE_SPECIFICATION_LONG or
                sys.argv[i] == DATABASE_SPECIFICATION_LONG or
                sys.argv[i] == FILTER_SPECIFICATION_LONG)):

            sys.argv[i + 1] = (
                str(sys.argv[i + 1])[:0] + " "
                + str(sys.argv[i + 1])[0:])

    args = parser.parse_args()

    # --- Job Control ---
    with drmaa.Session() as session:

        execution = Execution(session, args)

        # --- K-MER COUNTING ---
        inclusionKMerLocations, exclusionKMerLocations = countKMers(execution)

        # --- K-MER AGGREGATION ---
        aggregateKMers(
            execution, inclusionKMerLocations, exclusionKMerLocations)

        # --- SIGNATURE EXTRACTION ---
        candidateLocations = extractSignatures(execution)

        # --- SIGNATURE FILTERING ---
        filterSignatures(execution, candidateLocations)

        execution.produceReceipt()

if __name__ == '__main__':

    main()
