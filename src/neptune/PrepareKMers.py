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

"""
# =============================================================================

This script prepares k-mers for further analysis by producing a list of
lexicographically sorted and distinct k-mers.

INPUT:

GCT
ATA
TTT
GGG

OUTPUT:

AAA
AGC
ATA
CCC

USAGE:

script.py -h
script.py -i INPUT -o OUTPUT

EXAMPLE:

script.py -i dirty.kmers -o clean.kmers

# =============================================================================
"""

import os
import sys
import argparse

from Bio.Seq import Seq

"""
# =============================================================================

GLOBALS

# =============================================================================
"""

PROGRAM_DESCRIPTION = "Produces a sorted list of distinct k-mers. Only the \
    lexicographically first of a k-mer and its reverse complement is reported."

# ARGUMENTS #
LONG = "--"
SHORT = "-"

# REQUIRED #

# Input
INPUT = "input"
INPUT_LONG = LONG + INPUT
INPUT_SHORT = SHORT + "i"
INPUT_HELP = "The file location of the k-mers to prepare for processing."

# Output
OUTPUT = "output"
OUTPUT_LONG = LONG + OUTPUT
OUTPUT_SHORT = SHORT + "o"
OUTPUT_HELP = "The file location to write prepared k-mers."

"""
# =============================================================================

CONVERT
-------


PURPOSE
-------

Converts a k-mer or k-mer line to the lexicographically smaller of [kmer] and
its [reverseComplement(kmer)].


INPUT
-----

[STRING] [kmer]
    The k-mer to convert in nucleotide characters.


RETURN
------

[STRING] [kmer]
    The lexicographically smaller of [kmer] and [reverseComplement(kmer)].

# =============================================================================
"""
def convert(kmer):

    reverse = str(Seq(kmer).reverse_complement())
    reverse = reverse.upper()

    kmer = kmer if kmer < reverse else reverse

    return kmer


"""
# =============================================================================

MAIN

# =============================================================================
"""
def main(kmersInputFile, kmersOutputFile):

    kmers = []

    for line in kmersInputFile:

        kmer = line.split()[0].upper()
        kmers.append(convert(kmer))

    kmers = sorted(set(kmers))

    for kmer in kmers:
        output = str(kmer) + "\n"
        kmersOutputFile.write(output)


"""
# =============================================================================

PARSING

# =============================================================================
"""
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=PROGRAM_DESCRIPTION)

    parser.add_argument(
        INPUT_SHORT,
        INPUT_LONG,
        dest=INPUT,
        help=INPUT_HELP,
        type=str, required=True)

    parser.add_argument(
        OUTPUT_SHORT,
        OUTPUT_LONG,
        dest=OUTPUT,
        help=OUTPUT_HELP,
        type=str, required=True)

    args = parser.parse_args()

    kmersInputLocation = args.input
    kmersOutputLocation = args.output

    if not os.path.isfile(kmersInputLocation):
        sys.stderr.write("ERROR: Could not open input file.\n")
        exit(1)

    kmersInputFile = open(kmersInputLocation, 'r')
    kmersOutputFile = open(kmersOutputLocation, 'w')

    main(kmersInputFile, kmersOutputFile)

    kmersInputFile.close()
    kmersOutputFile.close()

    exit(0)
