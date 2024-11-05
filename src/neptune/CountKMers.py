#!/usr/bin/env python

"""
# =============================================================================

Copyright Government of Canada 2015-2024

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

Multiple output files are written only when the degree of organization is
specified. The number of output files is: (4^[organization]). Increasing the
degree of organization will not speed up k-mer counting. However, it may
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
script.py -k K -i INPUT -o OUTPUT [-p ORGANIZATION]

EXAMPLES:

script.py -k 21 -i path/to/reference.FASTA -o path/to/output.kmers
script.py -k 21 -i reference.FASTA -o output.kmers -p 3

# =============================================================================
"""

import argparse
import os
import operator

import neptune.Utility as Utility

"""
# =============================================================================

GLOBALS

# =============================================================================
"""

PROGRAM_DESCRIPTION = "Counts k-mers in a FASTA or multi FASTA file and \
    writes them to either one or multiple files. Only the lexicographically \
    smaller of a k-mer and its reverse is reported. The k-mers are reported \
    in sorted order."

# DEFAULTS

ORGANIZATION_DEFAULT = 0

# ARGUMENTS

LONG = "--"
SHORT = "-"

# REQUIRED ARGUMENTS #

# Input
INPUT = "input"
INPUT_LONG = LONG + INPUT
INPUT_SHORT = SHORT + "i"
INPUT_HELP = "The file location of a FASTA file from which to count k-mers."

# Output
OUTPUT = "output"
OUTPUT_LONG = LONG + OUTPUT
OUTPUT_SHORT = SHORT + "o"
OUTPUT_HELP = "The file location to write the k-mers."

# k-mer
KMER = "kmer"
KMER_LONG = LONG + KMER
KMER_SHORT = SHORT + "k"
KMER_HELP = "The size of the k-mers."

# OPTIONAL ARGUMENTS #

# Organization
ORGANIZATION = "organization"
ORGANIZATION_LONG = LONG + ORGANIZATION
ORGANIZATION_HELP = "The degree of k-mer organization in the output files.\
    This exploits the four-character alphabet of nucleotides to produce \
    several k-mer output files, with all k-mers in a file beginning with the \
    same short sequence of nucleotides. This parameter determines the number \
    of nucleotides to use and will produce 4^X output files, where X is the \
    number of nucleotides specified by this parameter. The number of output \
    files directly corresponds to the amount of parallelization in the k-mer \
    aggregation process."

"""
# =============================================================================

WRITE MULTIPLE FILES
--------------------


PURPOSE
-------

Writes the k-mers to several output files. The number of output files is
determined by the degree of organization.

The output file names and locations are determined according to the provided
base name and the degree of organization.


INPUT
-----

[(STRING, INT) ITERABLE] [kmers]
    The k-mers to write to several files.

[STRING] [outputLocation]
    The base file path to write the output files.

[INT >= 0] [organization]
    The degree of organization. This will produce 4^[organization] output
    files.


RETURN
------

[NONE]

POST
----

There will be several output files written. The number of output files will be
4^[organization] and they will be named according to the [outputLocation]
parameter.

These files will be appended with a sequence tag, determined automatically
based on the degree of [organization]. This sequence tag reveals which
nucleotides must match in the initial sequence of all k-mers in the file.

There will additionally be a file for all sequences which begin with special
characters (e.g. N).

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
-----------------


PURPOSE
-------

Writes the passed k-mers to the output location.


INPUT
-----

[(STRING, INT) ITERABLE] [kmers]
    The k-mers to write to files.

[STRING] [outputFile]
    A writable file-like object to write output.


RETURN
------

[NONE]


POST
----

The k-mers will written to the [outputFile] in the order they are iterated in
their [kmers] data structure.

# =============================================================================
"""
def writeSingleFile(kmers, outputFile):

    for kmer in kmers:

        outputFile.write(str(kmer[0]) + " " + str(kmer[1]) + "\n")


"""
# =============================================================================

COUNT
-----


PURPOSE
-------

Counts the k-mers in a file and outputs the k-mer counts.


INPUT
-----

[FILE LOCATION] [inputLocation]
    The location of the input in FASTA format.

[FILE LOCATION] [outputLocation]
    The output location of the k-mer counts.

[INT >= 1] [k]
    The k-mer size.

[INT >= 0] [organization]
    The degree of organization. This is responsible for the number of output
    files.


RETURN
------

[NONE]


POST
----

The k-mers located in the input will be output to the [outputLocation] in
sorted order.

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

    organization = parameters[ORGANIZATION] \
        if parameters[ORGANIZATION] else ORGANIZATION_DEFAULT

    count(inputLocation, outputLocation, k, organization)


"""
# =============================================================================

MAIN

# =============================================================================
"""
def main():

    # description
    parser = argparse.ArgumentParser(description=PROGRAM_DESCRIPTION)

    # input location
    parser.add_argument(
        INPUT_SHORT,
        INPUT_LONG,
        dest=INPUT,
        help=INPUT_HELP,
        type=str, required=True)

    # output location
    parser.add_argument(
        OUTPUT_SHORT,
        OUTPUT_LONG,
        dest=OUTPUT,
        help=OUTPUT_HELP,
        type=str, required=True)

    # k-mer size
    parser.add_argument(
        KMER_SHORT,
        KMER_LONG,
        dest=KMER,
        help="k-mer size",
        type=int, required=True)

    # organization
    parser.add_argument(
        ORGANIZATION_LONG,
        dest=ORGANIZATION,
        help=ORGANIZATION_HELP,
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
