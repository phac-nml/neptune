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

This script counts k-mers in a FASTA or multi FASTA file and writes them to
either one or multiple k-mer files. Only the lexicographically smaller of a
k-mer and it's reverse complement is reported. The k-mers are reported in
sorted order.

Multiple output files are written only when the degree of parallelization is
specified. The number of output files is: (4^[parallelization]). Increasing the
degree of parallelization will not speed up k-mer counting. However, it may
improve the speed of tools which use these k-mers, such as AggregateKMers.

INPUT:

reference.fasta

OUTPUT (single file):

reference.kmers

OUTPUT (multiple files):

reference.kmers.A
reference.kmers.C
reference.kmers.G
reference.kmers.T

Such that all k-mers in reference.kmers.A begin with "A".

script.py -h
script.py -k K -i INPUT -o OUTPUT [-p PARALLELIZATION]

EXAMPLES:

script.py -k 21 -i path/to/reference.FASTA -o path/to/output.kmers
script.py -k 21 -i reference.FASTA -o output.kmers -p 3

# =============================================================================
"""

import argparse
import os
import operator

import Utility

"""
# =============================================================================

GLOBALS

# =============================================================================
"""

# NAMES
INPUT = "input"
OUTPUT = "output"
KMER = "kmer"
PARALLEL = "parallelization"

# ARGUMENTS
LONG = "--"

INPUT_LONG = LONG + INPUT
OUTPUT_LONG = LONG + OUTPUT
KMER_LONG = LONG + KMER
PARALLEL_LONG = LONG + PARALLEL

SHORT = "-"

INPUT_SHORT = SHORT + "i"
OUTPUT_SHORT = SHORT + "o"
KMER_SHORT = SHORT + "k"
PARALLEL_SHORT = SHORT + "p"

# DEFAULTS

PARALLELIZATION_DEFAULT = 0

"""
# =============================================================================

WRITE MULTIPLE FILES

PURPOSE:
    Writes the k-mers to several output files. The number of output files is
    determined by the degree of organization.

    The output file names and locations are determined according to the
    provided base name and the degree of organization.

INPUT:
    [(STRING, INT) ITERABLE] [kmers] - The k-mers to write to several files.
    [STRING] [outputLocation] - The base file path to write the output files.
    [INT >= 0] [organization] - The degree of organization.
        This will produce 4^[organization] output files.

RETURN:
    [NONE]

POST:
    There will be several output files written. The number of output files
    will be 4^[organization] and they will be named according to the
    [outputLocation] parameter.

    These files will be appended with a sequence tag, determined automatically
    based on the degree of [organization].

    There will additionally be a file for all sequences which begin with
    special characters.

# =============================================================================
"""
def writeMultipleFiles(kmers, outputLocation, organization):

    outputFiles = {}
    tags = Utility.getAggregationTags(organization)

    # initialize output files
    for tag in tags:

        outputName = outputLocation + "." + tag
        outputFiles[tag] = open(outputName, 'w')

    # write k-mers to output
    for item in kmers:

        kmer = str(item[0])
        count = str(item[1])

        tag = kmer[:organization]

        if tag in outputFiles:
            outputFiles[tag].write(kmer + " " + count + "\n")

        # handle special characters
        else:
            outputFiles[Utility.AGGREGATE_OTHER].write(
                kmer + " " + count + "\n")

    # close files
    for item in outputFiles:

        outputFiles[item].close()

"""
# =============================================================================

WRITE SINGLE FILE

PURPOSE:
    Writes the passed k-mers to the output location.

INPUT:
    [(STRING, INT) ITERABLE] [kmers] - The k-mers to write to files.
    [STRING] [outputFile] - A writable file-like object to write output.

RETURN:
    [NONE]

POST:
    The k-mers will written to the [outputFile] in the order they are
    iterated in their data structure.

# =============================================================================
"""
def writeSingleFile(kmers, outputFile):

    for kmer in kmers:

        outputFile.write(str(kmer[0]) + " " + str(kmer[1]) + "\n")

"""
# =============================================================================

COUNT

PURPOSE:
    Counts the k-mers in a file and outputs the k-mer counts.

INPUT:
    [FILE LOCATION] [inputLocation] - The location of the input.
    [FILE LOCATION] [outputLocation] - The output location.
    [INT >= 0] [k] - The k-mer size.
    [INT >= 0] [organization] - The degree of organization.
        This is responsible for the number of output files.

RETURN:
    [NONE]

POST:
    The k-mers located in the input will be output to the [outputLocation].

# =============================================================================
"""
def count(inputLocation, outputLocation, k, organization):

    # check input file
    if not os.path.isfile(inputLocation):
        raise RuntimeError(
            "ERROR: Could not open input file: " + inputLocation + "\n")

    inputFile = open(inputLocation, 'r')

    references = Utility.buildReferences(inputFile)

    kmers = {}

    # iterate all references
    for ref in references:

        # next reference
        reference = references[ref]

        # every kmer in reference
        for i in range(len(reference.strip()) - k + 1):

            # k-mer and reverse complement
            kmer = reference[i:i + k]
            reverse = Utility.reverseComplement(kmer)

            kmer = min(kmer, reverse)

            if kmer in kmers:
                kmers[kmer] += 1

            else:
                kmers[kmer] = 1

    # sort k-mers
    sortedKMers = sorted(kmers.items(), key=operator.itemgetter(0))

    # write k-mers out
    if organization == 0:
        outputFile = open(outputLocation, 'w')
        writeSingleFile(sortedKMers, outputFile)
        outputFile.close()

    else:
        writeMultipleFiles(sortedKMers, outputLocation, organization)

    # close input file
    inputFile.close()

"""
# =============================================================================

PARSE

# =============================================================================
"""
def parse(parameters):

    inputLocation = parameters[INPUT]
    outputLocation = parameters[OUTPUT]
    k = parameters[KMER]

    parallelization = parameters[PARALLEL] \
        if parameters[PARALLEL] else PARALLELIZATION_DEFAULT

    count(inputLocation, outputLocation, k, parallelization)

"""
# =============================================================================

MAIN

# =============================================================================
"""
def main():

    # description
    parser = argparse.ArgumentParser(
        description='Counts k-mers in a FASTA or multi FASTA file and \
        writes them to either one or multiple files. Only the \
        lexicographically smaller of a k-mer and its reverse is \
        reported. The k-mers are reported in sorted order.')

    # input location
    parser.add_argument(
        INPUT_SHORT,
        INPUT_LONG,
        dest=INPUT,
        help="input FASTA file location",
        type=str, required=True)

    # output location
    parser.add_argument(
        OUTPUT_SHORT,
        OUTPUT_LONG,
        dest=OUTPUT,
        help="output file location",
        type=str, required=True)

    # k-mer size
    parser.add_argument(
        KMER_SHORT,
        KMER_LONG,
        dest=KMER,
        help="k-mer size",
        type=int, required=True)

    # parallelization
    parser.add_argument(
        PARALLEL_SHORT,
        PARALLEL_LONG,
        dest=PARALLEL,
        help="degree of parallelization; produces 4^[parallelization] files",
        type=int)

    args = parser.parse_args()
    parameters = vars(args)
    parse(parameters)

"""
# =============================================================================
# =============================================================================
"""
if __name__ == '__main__':

    main()
