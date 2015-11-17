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
Date: 30 April 2015

This script provides shared utility to other scripts.

# =============================================================================
"""

import math

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

AGGREGATE_OTHER = "__OTHER__"

"""
# =============================================================================

GET AGGREGATION TAGS

PURPOSE:
    Produces a list of aggregation tags associate with the specified degree
    of parallelization.

INPUT:
    [INT >= 0] [parallelization] - The degree of parallelization.

RETURN:
    [STRING ITERABLE] [tags] - An iterable object of sequence tags.

        There will be (4^[parallelization]) tags produced and an
        "AGGREGATE_OTHER" tag.

        These tags will contain A, C, G, and T, and arranged in lexicographic
        order with the "AGGREGATE_OTHER" tag at the very end.

# =============================================================================
"""
def getAggregationTags(parallelization):

    if parallelization < 0:
        raise RuntimeError("The degree of parallelization is out of range.")

    instances = int(math.pow(4.0, float(parallelization)))
    tags = []

    # initialize
    for i in range(instances):

        tag = generateSequence(i, parallelization)
        tags.append(tag)

    tags.append(AGGREGATE_OTHER)

    return tags

"""
# =============================================================================

GENERATE SEQUENCE

PURPOSE:
    Generates and returns a nucleotide sequence of a specified
    length from a passed integer value.

INPUT:
    [INT >= 0] [integer] - The integer to convert to sequence.
    [INT >= 0] [length] - The length of the sequence to return.

RETURN:
    [STRING] [sequence] - A sequence of length [length] generated from the
        integer [integer].

        Sequences are genrated such that, for all 0 <= x < 4^[length]:

        generateSequence(x, [length]) < generateSequence(x+1, [length])

        when considered lexicographically.

# =============================================================================
"""
def generateSequence(integer, length):

    if integer < 0:
        raise RuntimeError("The sequence integer is out of range.")

    if length < 0:
        raise RuntimeError("The sequence length is out of range.")

    current = integer
    sequence = ""

    for i in range(length):

        base = current % 4

        if base == 0:
            sequence = "A" + sequence

        elif base == 1:
            sequence = "C" + sequence

        elif base == 2:
            sequence = "G" + sequence

        elif base == 3:
            sequence = "T" + sequence

        current = current / 4

    return sequence

"""
# =============================================================================

REVERSE COMPLEMENT

PURPOSE:
    Produces and returns the reverse complement of the sequence.

INPUT:
    [STRING] [sequence] - The sequence to reverse complement.

RETURN:
    [STRING] [reverse] - The reverse complement of the passed sequence.

# =============================================================================
"""
def reverseComplement(sequence):

    reverse = str(Seq(sequence, generic_dna).reverse_complement())
    return reverse

"""
# =============================================================================

BUILD REFERENCES

PURPOSE:
    Builds string references from the reference file.

INPUT:
    [FILE] [referenceFile] - The file from which to build the string reference.

RETURN:
    [STRING ITERABLE] [references] - A list of string references.

# =============================================================================
"""
def buildReferences(referenceFile):

    references = {}

    # build references
    for line in referenceFile:

        # new reference:
        if line[0] == ">":
            tokens = (line[1:]).split()
            referenceName = tokens[0]

            references[referenceName] = ""

            print "Building: " + str(referenceName)

        # continue building reference:
        else:
            references[referenceName] += line.strip().upper()

    return references

"""
# =============================================================================

ESTIMATE REFERENCE PARAMETERS

PURPOSE:
    Estimates the size of the reference and GC-content from one or more
    reference fragments.

INPUT:
    [STRING ITERABLE] [references] - A list of string references.

RETURN:
    [TUPLE: INT, FLOAT] [size, gcContent] - An estimate for the reference size
        and GC-content.

# =============================================================================
"""
def estimateReferenceParameters(references):

    if references is None or len(references) < 1:
        raise RuntimeError("There are no references.")

    sumGC = 0
    sumAT = 0
    size = 0

    for reference in references:

        sumGC += references[reference].count('G') \
            + references[reference].count('C')

        sumAT += references[reference].count('A') \
            + references[reference].count('T')

        size += len(references[reference])

    gcContent = float(sumGC) / float(sumGC + sumAT)

    return size, gcContent
