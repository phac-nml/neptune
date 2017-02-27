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

This script aggregates the k-mers from one or more inclusion files with one or
more exclusion files. The output provides a count of the number of distinct
k-mer observaions per each inclusion and exclusion files.

The input files must have one k-mer per line and each k-mer must be preceded
by no spaces or characters. Any characters, including k-mer counts, following
the k-mers will be ignored.

If the delete flag is used, then all input files will be deleted after they
aggregated.

INPUT (one file):

AAAAA
AAAAC
AAAAG
AAAAT

The output will be in the following format:

[k-mer] [inclusion counts] [exclusion counts]

OUTPUT:

AAAAA 3 1
AAAAC 1 2
AAAAG 3 3
AAAAT 3 0

USAGE:

script.py -h
script.py -i [INCLUSION] [...] -e [EXCLUSION] [...] -o [OUTPUT] [--delete]

EXAMPLE:

script.py -i inclusion1.kmers inclusion2.kmers -e exclusion1.kmers -o out.kmers
script.py -i inclusion1.kmers -e exclusion1.kmers -o out.kmers --delete

# =============================================================================
"""

import os
import argparse

"""
# =============================================================================

GLOBALS

# =============================================================================
"""

# NAMES
INCLUSION = "inclusion"
EXCLUSION = "exclusion"
OUTPUT = "output"
DELETE = "delete"

# ARGUMENTS
LONG = "--"

INCLUSION_LONG = LONG + INCLUSION
EXCLUSION_LONG = LONG + EXCLUSION
OUTPUT_LONG = LONG + OUTPUT
DELETE_LONG = LONG + DELETE

SHORT = "-"

INCLUSION_SHORT = SHORT + "i"
EXCLUSION_SHORT = SHORT + "e"
OUTPUT_SHORT = SHORT + "o"

# DEFAULTS

DELETE_DEFAULT = False

"""
# =============================================================================

FIND SMALLEST

PURPOSE:
    Locates the lexicographically smallest element in a list of strings.

    The function ignores empty strings and uses a sentinel value to determine
    if all the strings are empty. This sentinel value must always be
    lexicographically larger than all string values.

INPUT:
    [STRING ITERABLE] [strings] - An iterable strings object.
    [STRING] [SENTINEL] - A sentinel value, which must be lexicographically
        larger than all [strings].

RETURN:
    [STRING] [smallest] - The lexicographically smallest value of all the
        strings or SENTINEL if all strings are empty.

# =============================================================================
"""
def findSmallest(strings, SENTINEL):

    smallest = SENTINEL     # the lexicographically smallest observed string

    for string in strings:
        if string != "" and string < smallest:
            smallest = string

    return smallest


"""
# =============================================================================

AGGREGATE KMER

PURPOSE:
    Computes the counts of the k-mer, determined by the number of observations
    across all files.

    This function assumes: len([kmers]) == len([files])

    The function observes all the k-mers and compares them with the k-mer
    parameter. When there is a match, the function increases the count and
    advances the file associated with the k-mer.

    The [kmers] parameter corresponds to the heads of all the k-mer files. The
    files MUST necessarily be read and advanced when there is a k-mer match
    found in a corresponding [kmers] array. This is because each k-mer in each
    file is only ever observed once.

INPUT:
    [STRING] [kmer] - The k-mer to compare against all other k-mers and count
        when observed.
    [STRING LIST] [kmers] - A list of strings, which should understood as the
        head k-mer of files.
    [FILE LIST] [files] - A list of open files associated with the [kmers]
        list. It is assumed: len([kmers]) == len([files])

RETURN:
    [INT >= 0] [count] - The number of exact [kmer] matches found in the list
        of [kmers].

# =============================================================================
"""
def aggregateKMer(kmer, kmers, files):

    count = 0

    # iterate over all k-mers
    for i in range(len(kmers)):

        # check for empty string
        if kmers[i] != "" and kmers[i] == kmer:

            count += 1

            # advance file
            line = files[i].readline()

            # check for end of file and assign next k-mer
            kmers[i] = line.split()[0] if line.split() else line

    return count


"""
# =============================================================================

AGGREGATE

PURPOSE:
    Aggregates the k-mers in the inclusion and exclusion files and produces a
    file containing the k-mers and their inclusion and exclusion counts. The
    input inclusion and exclusion files must contain only distinct and
    lexicographically sorted k-mers.

INPUT:
    [FILE LIST] [inclusionLocations] - The list of openable inclusion k-mer
        files locations.
    [FILE LIST] [exclusionLocations] - The list of openable exclusion k-mer
        files locations.
    [FILE] [outputLocation] - The file location to write the aggregated k-mers.
    [BOOL] [delete] - Whether or not to delete the [inclusionLocations] and
        [exclusionLocations] after aggregation is complete.

    NOTE: The input files must contain only distinct and lexicographically
    sorted k-mers. These k-mers must appear first on every line and be
    preceded by no spaces or special characters. Any characters following
    the k-mers, including k-mer counts, will be ignored.

POST:
    The k-mers and their aggregate counts value will be written to a file at
    the [outputLocation].

# =============================================================================
"""
def aggregate(inclusionLocations, exclusionLocations, outputLocation, delete):

    SENTINEL = "~"              # sentinel value

    inclusionKMers = []         # current k-mer of inclusion files
    exclusionKMers = []         # current k-mer of exclusion files

    # open files
    inclusionFiles = []
    exclusionFiles = []

    # open inclusion files
    for location in inclusionLocations:

        if not os.path.isfile(location):
            raise RuntimeError(
                "ERROR: Could not open inclusion file: " +
                str(location) + "\n")

        inclusionFiles.append(open(location, 'r'))

    # open exclusion files
    for location in exclusionLocations:

        if not os.path.isfile(location):
            raise RuntimeError(
                "ERROR: Could not open exclusion file: " +
                str(location) + "\n")

        exclusionFiles.append(open(location, 'r'))

    outputFile = open(outputLocation, 'w')

    # initialize k-mers:
    for inclusionFile in inclusionFiles:

        line = inclusionFile.readline()

        # check for end of file and assign next k-mer
        kmer = line.split()[0] if line.split() else line
        inclusionKMers.append(kmer)

    for exclusionFile in exclusionFiles:

        line = exclusionFile.readline()

        # check for end of file and assign next k-mer
        kmer = line.split()[0] if line.split() else line
        exclusionKMers.append(kmer)

    # aggregate values:
    while(1):

        kmer = findSmallest((inclusionKMers + exclusionKMers), SENTINEL)

        # inclusion aggregation:
        incounts = aggregateKMer(kmer, inclusionKMers, inclusionFiles)

        # exclusion aggregation:
        excounts = aggregateKMer(kmer, exclusionKMers, exclusionFiles)

        # are all files at end of file?
        if kmer == SENTINEL:
            break

        # write aggregated k-mer to output
        outputString = str(kmer) + " " + str(incounts) + " " + str(excounts)
        outputFile.write(outputString + "\n")

    # close files
    for inclusion in inclusionFiles:

        inclusion.close()

    for exclusion in exclusionFiles:

        exclusion.close()

    outputFile.close()

    # delete input files
    if delete:

        for filename in inclusionLocations + exclusionLocations:

            if os.path.exists(filename):
                    os.remove(filename)


"""
# =============================================================================

PARSE

# =============================================================================
"""
def parse(parameters):

    inclusionLocations = parameters[INCLUSION]
    exclusionLocations = parameters[EXCLUSION]
    outputLocation = parameters[OUTPUT]

    delete = parameters[DELETE] \
        if parameters[DELETE] else DELETE_DEFAULT

    # aggregate
    aggregate(inclusionLocations, exclusionLocations, outputLocation, delete)


"""
# =============================================================================

MAIN

# =============================================================================
"""
def main():

    # description
    parser = argparse.ArgumentParser(
        description='Aggregates one or more inclusion k-mer files \
        with one or more exclusion files. The number of distinct \
        inclusion and exclusion k-mer observations per file will be \
        reported immediately following each k-mer in the output.')

    # inclusion k-mers
    parser.add_argument(
        INCLUSION_SHORT,
        INCLUSION_LONG,
        dest=INCLUSION,
        help="inclusion k-mer file(s)",
        type=str, required=True, nargs='+')

    # exclusion k-mers
    parser.add_argument(
        EXCLUSION_SHORT,
        EXCLUSION_LONG,
        dest=EXCLUSION,
        help="exclusion k-mer file(s)",
        type=str, required=True, nargs='+')

    # output k-mers
    parser.add_argument(
        OUTPUT_SHORT,
        OUTPUT_LONG,
        dest=OUTPUT,
        help="output file path", type=str,
        required=True)

    # delete input files
    parser.add_argument(
        DELETE_LONG,
        dest=DELETE,
        help="delete input flag",
        action='store_true')

    args = parser.parse_args()
    parameters = vars(args)
    parse(parameters)


"""
# =============================================================================
# =============================================================================
"""
if __name__ == '__main__':

    main()
