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
Date: 17 April 2015

This script extracts candidate signatures from one or more inclusion files
within the context of one or more exclusion files. The signatures are extracted
from a reference using k-mers generated from all inputs.

INPUT (aggregated k-mers):

[k-mer] [inclusion counts] [exclusion counts]

AAAAA 3 1
AAAAC 1 2
AAAAG 3 3
AAAAT 3 0

USAGE:

script.py -h
script.py -r REFERENCE -i INCLUSION [INCLUSION ...]
        -e EXCLUSION [EXCLUSION ...]
        -k KMERS -o OUTPUT

EXAMPLE:

script.py -r inclusion1.fasta -i inclusion/* -e exclusion/*
    -k aggregated.kmers -o candidates.out

# =============================================================================
"""

import math
import argparse
import os

from Utility import reverseComplement
from Utility import buildReferences
from Utility import estimateReferenceParameters

from scipy.stats import norm

"""
# =============================================================================

GLOBALS

# =============================================================================
"""

# NAMES
REFERENCE = "reference"
REFERENCE_SIZE = "referenceSize"
RATE = "rate"
INCLUSION = "inclusion"
INHITS = "inhits"
EXCLUSION = "exclusion"
EXHITS = "exhits"
GAP = "gap"
SIZE = "size"
GC_CONTENT = "gcContent"
CONFIDENCE = "confidence"
KMERS = "kmers"
OUTPUT = "output"

# ARGUMENTS
LONG = "--"

REFERENCE_LONG = LONG + REFERENCE
REFERENCE_SIZE_LONG = LONG + REFERENCE_SIZE
RATE_LONG = LONG + RATE
INCLUSION_LONG = LONG + INCLUSION
INHITS_LONG = LONG + INHITS
EXCLUSION_LONG = LONG + EXCLUSION
EXHITS_LONG = LONG + EXHITS
GAP_LONG = LONG + GAP
SIZE_LONG = LONG + SIZE
GC_LONG = LONG + GC_CONTENT
CONFIDENCE_LONG = LONG + CONFIDENCE
KMERS_LONG = LONG + KMERS
OUTPUT_LONG = LONG + OUTPUT

SHORT = "-"

REFERENCE_SHORT = SHORT + "r"
REFERENCE_SIZE_SHORT = SHORT + "rs"
RATE_SHORT = SHORT + "q"
INCLUSION_SHORT = SHORT + "i"
INHITS_SHORT = SHORT + "ih"
EXCLUSION_SHORT = SHORT + "e"
EXHITS_SHORT = SHORT + "eh"
GAP_SHORT = SHORT + "g"
SIZE_SHORT = SHORT + "s"
GC_SHORT = SHORT + "gc"
CONFIDENCE_SHORT = SHORT + "c"
KMERS_SHORT = SHORT + "k"
OUTPUT_SHORT = SHORT + "o"

"""
# =============================================================================

REGION

# =============================================================================
"""
class Region():

    def __init__(self, sequence, reference, position):

        self.sequence = sequence
        self.reference = reference
        self.position = position

"""
# =============================================================================

EXTRACT

PURPOSE:
    Extracts candidate signatures from the set of references.

    The extraction process is guided by the inclusion and exclusion k-mer
    dictionaries. This function is concerned only about k-mer existence in
    those dictionaries.

INPUT:
    [STRING ITERABLE] [references] - A set of string references. This is
        intended to correspond to a single or multi-fasta file.
    [INT >= 1] [k] - The k-mer size.
    [STRING DICTIONARY] [inmers] - The inclusion k-mers dictionary.
    [STRING DICTIONARY] [exmers] - The inclusion k-mers dictionary.
    [INT >= 1] [size] - The minimum signature size.
    [INT >= 1] [gap] - The maximum allowable gap size.
    [FILE] [outputFile] = The output file to write candidate signatures.

POST:
    The candidate signatures will be written to the output file.

# =============================================================================
"""
def extract(references, k, inmers, exmers, size, gap, outputFile):

    # references
    if references is None or len(references) < 1:
        raise RuntimeError("There are no references.")

    # 1 <= kmerSize
    if k < 1:
        raise RuntimeError("The k-mer size is out of range.")

    # inclusion k-mers
    if inmers is None or len(inmers) < 1:
        raise RuntimeError("There are no inclusion k-mers.")

    # exclusion k-mers
    if exmers is None:
        raise RuntimeError("There are no exclusion k-mers.")

    # 1 <= size
    if size < 1:
        raise RuntimeError("The signature size is out of range.")

    # 1 <= gap
    if gap < 1:
        raise RuntimeError("The gap size is out of range.")

    # output
    if outputFile is None:
        raise RuntimeError("The output location is not specified.")

    regions = []

    # iterate all references
    for key in references:

        # next reference
        ref = references[key]

        # initialize positions
        start = -1
        end = -1

        # every kmer in reference
        for i in range(len(ref.strip()) - k + 1):

            # k-mer and reverse complement
            kmer = ref[i:i + k]
            reverse = reverseComplement(kmer)

            # kmer is in exclusion sufficiently -- break chain
            if kmer in exmers or reverse in exmers:

                # close the region if started:
                if (end - start) >= size:
                    region = Region(ref[start:end], key, start)
                    regions.append(region)

                # end the region regardless:
                start = -1
                end = -1

            # k-mer is in inclusion sufficiently -- build chain
            # (else -- don't both break and build)
            elif kmer in inmers or reverse in inmers:

                # new chain
                if start < 0 and end < 0:
                    start = i + k - 1
                    end = i + 1

                # hit with something else started
                if start >= 0 and end >= 0:

                    # gap within size? -- yes
                    if (i - (end + 1)) <= gap:
                        end = i + 1

                    # gap within size? -- no
                    else:
                        if (end - start) >= size:
                            region = Region(ref[start:end], key, start)
                            regions.append(region)

                        start = i + k - 1
                        end = i + 1

        if start >= 0 and end > 0 and (end - start) >= size:
            region = Region(ref[start:end], key, start)
            regions.append(region)

    for i in range(len(regions)):

        outputFile.write(
            ">" + str(i) + " " + str(len(regions[i].sequence))
            + " " + str(regions[i].reference) + " "
            + str(regions[i].position) + "\n")

        outputFile.write(regions[i].sequence + "\n")

"""
# =============================================================================

CALCULATE PROBABILITY THAT HOMOLOGOUS BASES MUTATE AND MATCH

P(X_M = Y_M)_H

PURPOSE:
    Calculates the probability that homologous bases both mutate to the same
    matching base, given a certain GC-content environment.

INPUT:
    [0 <= FLOAT <= 1] [GC] - The GC-content of the environment.

RETURN:
    [0 <= FLOAT <= 1] [P(X_M = Y_M)_H] - The proability that homologous bases
        mutate and match.

# =============================================================================
"""
def calculateProbHBMM(GC):

    # 0 <= GC <= 1
    if GC < 0 or GC > 1:
        raise RuntimeError("The GC-content is out of range.")

    GC = float(GC)

    # P(X_M = Y_M)_H
    probHBMM = (
        2 * math.pow(GC / (GC + 1), 2)
        + math.pow((1 - GC) / (GC + 1), 2)) * (1 - GC) \
        + (2 * math.pow((1 - GC) / (2 - GC), 2)
            + math.pow(GC / (2 - GC), 2)) * (GC)

    return probHBMM

"""
# =============================================================================

CALCULATE PROBABILITY THAT HOMOLOGOUS BASES MATCH

P(X = Y)_H

PURPOSE:
    Calculates the probability that homologous bases match, given a certain
    mutation rate and a GC-content environment.

INPUT:
    [0 <= FLOAT <= 1] [mutationRate] - The probability of an arbitrary base
        mutating.
    [0 <= FLOAT <= 1] [GC] - The GC-content of the environment.

RETURN:
    [0 <= FLOAT <= 1] [P(X = Y)_H] - The probability that homologous bases
        match.

# =============================================================================
"""
def calculateProbHBM(mutationRate, GC):

    # 0 <= mutationRate <= 1
    if mutationRate < 0 or mutationRate > 1:
        raise RuntimeError("The mutationRate is out of range.")

    # 0 <= GC <= 1
    if GC < 0 or GC > 1:
        raise RuntimeError("The GC-content is out of range.")

    M = float(mutationRate)

    # P(X_M = Y_M)_H
    probHBMM = calculateProbHBMM(GC)

    # P(X = Y)_H
    probHBM = (1 - M) * (1 - M) + (M) * (M) * probHBMM

    return probHBM

"""
# =============================================================================

CALCULATE PROBABILITY THAT HOMOLOGOUS K-MERS MATCH

P(k_X = k_Y)_H

PURPOSE:
    Calculates the probability that homologous k-mers match, given a certain
    mutation rate, GC-content, and k-mer size.

INPUT:
    [0 <= FLOAT <= 1] [mutationRate] - The probability of an arbitrary base
        mutating.
    [0 <= FLOAT <= 1] [GC] - The GC-content of the environment.
    [INT >= 1] [kmerSize] - The size of the k-mers.

RETURN:
    [0 <= FLOAT <= 1] [P(k_X = k_Y)_H] - The probability that homologous k-mers
        match.

# =============================================================================
"""
def calculateProbHKM(mutationRate, GC, kmerSize):

    # 0 <= mutationRate <= 1
    if mutationRate < 0 or mutationRate > 1:
        raise RuntimeError("The mutationRate is out of range.")

    # 0 <= GC <= 1
    if GC < 0 or GC > 1:
        raise RuntimeError("The GC-content is out of range.")

    # 1 <= kmerSize
    if kmerSize < 1:
        raise RuntimeError("The k-mer size is out of range.")

    # P(X = Y)_H
    probHBM = calculateProbHBM(mutationRate, GC)

    # P(k_X = k_Y)_H
    probHKM = math.pow(probHBM, kmerSize)

    return probHKM

"""
# =============================================================================

ESTIMATE SIGNATURE SIZE

PURPOSE:
    Estimates the minimum candidate size.

INPUT:
    [INT >= 1] [kmerSize] - The size of the k-mers.

RETURN:
    [INT] [estimate] - An estimate of the minimum candidate size.

# =============================================================================
"""
def estimateSignatureSize(kmerSize):

    FACTOR = 4.0

    estimate = FACTOR * kmerSize

    return int(estimate)

"""
# =============================================================================

ESTIMATE GAP SIZE

PURPOSE:
    This function estimates the maximum gap size before abandoning a candidate
    region.

INPUT:
    [0 <= FLOAT <= 1] [mutationRate] - The probability of an arbitrary base
        mutating.
    [0 <= FLOAT <= 1] [GC] - The GC-content of the environment.
    [INT >= 1] [kmerSize] - The size of the k-mers.
    [0 < FLOAT < 1] [confidence] - The statistical confidence.

RETURN:
    [INT] [estimate] - An estimate of the maximum allowable gap size. This will
        be rounded up to the next integer.

# =============================================================================
"""
def estimateGapSize(mutationRate, GC, kmerSize, confidence):

    # 0 <= mutationRate <= 1
    if mutationRate < 0 or mutationRate > 1:
        raise RuntimeError("The mutationRate is out of range.")

    # 0 <= GC <= 1
    if GC < 0 or GC > 1:
        raise RuntimeError("The GC-content is out of range.")

    # 1 <= kmerSize
    if kmerSize < 1:
        raise RuntimeError("The k-mer size is out of range.")

    # 0 < confidence < 1
    if confidence <= 0 or confidence >= 1:
        raise RuntimeError("The statistical confidence is out of range.")

    p = calculateProbHBM(mutationRate, GC)
    q = float(1 - p)
    k = float(kmerSize)
    c = float(confidence)

    # Feller (recurrence times of length k)
    mean = (1 - math.pow(p, k)) / (q * math.pow(p, k))
    variance = 1 / math.pow((q * math.pow(p, k)), 2) - (2 * k + 1) \
        / (q * math.pow(p, k)) - p / math.pow(q, 2)
    stdev = math.sqrt(variance)

    # Chebyshev's Inequality
    deviations = math.sqrt(1 / (1 - c))
    estimate = math.ceil(mean + deviations * stdev)

    return int(estimate)

"""
# =============================================================================

ESTIMATE EXCLUSION HITS

PURPOSE:
    This function estimates the minimum number of exclusion hits to prevent
    candidate building.

INPUT:
    [INT >= 0] [totalExclusion] - The total number of exclusion targets.
    [0 <= FLOAT <= 1] [mutationRate] - The probability of an arbitrary base
        mutating.
    [INT >= 1] [kmerSize] - The size of the k-mers, such that 1 <= kmerSize.

RETURN:
    [INT >= 0] [estimate] - An estimate of the minimum number of exclusion
        k-mer hits before rejection.

# =============================================================================
"""
def estimateExclusionHits(totalExclusion, rate, kmerSize):

    estimate = 1.0

    return estimate

"""
# =============================================================================

ESTIMATE INCLUSION HITS

PURPOSE:
    This function estimates the minimum number of inclusion hits required for
    confident candidate building.

INPUT:
    [INT >= 0] [totalInclusion] - The total number of inclusion targets.
    [0 <= FLOAT <= 1] [mutationRate] - The probability of an arbitrary base
        mutating.
    [0 <= FLOAT <= 1] [GC] - The GC-content of the environment.
    [INT >= 1] [kmerSize] - The size of the k-mers.
    [0 < FLOAT < 1] [confidence] - The statistical confidence.

RETURN:
    [INT >= 0] [estimate] - An estimate of the minimum number of inclusion
        k-mer hits.

# =============================================================================
"""
def estimateInclusionHits(
        totalInclusion, mutationRate, GC, kmerSize, confidence):

    # 0 <= mutationRate <= 1
    if mutationRate < 0 or mutationRate > 1:
        raise RuntimeError("The mutationRate is out of range.")

    # 0 <= GC <= 1
    if GC < 0 or GC > 1:
        raise RuntimeError("The GC-content is out of range.")

    # 1 <= kmerSize
    if kmerSize < 1:
        raise RuntimeError("The k-mer size is out of range.")

    # 0 < confidence < 1
    if confidence <= 0 or confidence >= 1:
        raise RuntimeError("The statistical confidence is out of range.")

    deviations = norm.ppf(confidence)   # (Percent Point Function - Normal)

    p = calculateProbHKM(mutationRate, GC, kmerSize)
    q = 1 - p
    N = float(totalInclusion)

    # Binomial Distribution
    mean = (N - 1) * p
    variance = (N - 1) * p * q
    stdev = math.sqrt(variance)

    estimate = math.floor(1 + mean - deviations * stdev)
    estimate = max(0.0, estimate)

    return estimate

"""
# =============================================================================

ESTIMATE K

PURPOSE:
    This function estimates the size of k from the k-mers file.

INPUT:
    [FILE] [kmerFile] - The file of aggregated k-mers.

RETURN:
    [INT >= 1] [estimate] - An estimate of the k-mer size.

# =============================================================================
"""
def estimateK(kmerFile):

    line = kmerFile.readline()
    kmer = line.split()[0]

    estimate = len(kmer)

    return estimate

"""
# =============================================================================

BUILD K-MERS

PURPOSE:
    This function fills the inclusion and exclusion k-mer dictionaries from a
    single aggregated k-mer file.

INPUT:
    [FILE] [kmerFile] - The file of aggregated k-mers.
    [STRING DICTIONARY] [inmers] - The inclusion k-mer dictionary to fill with
        k-mers.
    [STRING DICTIONARY] [exmers] - The exclusion k-mer dictionary to fill with
        k-mers.
    [INT >= 0] [inhits] - The minimum number of inclusion k-mers for a
        candidate.
    [INT >= 0] [exhits] - The maximum number of exclusion k-mers for a
        candidate.

POST:
    The inclusion and exclusion k-mer dictionaries will be filled with all
    k-mers found in the k-mers file with at least inhits and at most exhits
    counts for the inclusion and exclusion k-mers, respectively.

# =============================================================================
"""
def buildKMers(kmerFile, inmers, exmers, inhits, exhits):

    for line in kmerFile:

        tokens = line.split()

        kmer = tokens[0].strip()
        incount = int(tokens[1].strip())
        excount = int(tokens[2].strip())

        if incount >= inhits:
            inmers[kmer] = incount

        if excount >= exhits:
            exmers[kmer] = excount

"""
# =============================================================================

REPORT PARAMETERS

PURPOSE:
    This function outputs the parameters to standard output.

INPUT:
    [STRING] [reportLocation] - The location of the report file.
    [FILE] [referenceFile] - The single reference to extract candidates from.
    [INT >= 0] [referenceSize] - The size of the reference.
    [0 <= FLOAT <= 1] [rate] - The rate of mutations and/or errors.
    [INT >= 0] [totalInclusion] - The number of inclusion genome files.
    [INT >= 0] [totalExclusion] - The number of exclusion genome files.
    [INT >= 0] [inhits] - The minimum number of inclusion k-mer hits.
    [INT >= 0] [exhits] - The maximum number of exclusion k-mer hits.
    [INT >= 1] [k] - The size of the k-mer.
    [FILE] kmerFile - The file containing aggregated k-mers.
    [INT >= 1] gap - The maximum inclusion k-mer gap size.
    [INT >= 1] size - The minimum size of any candidate.
    [0 <= FLOAT <= 1] GC - The GC-content of the environment.

POST:
    The parameters will be written to the report file location.

# =============================================================================
"""
def reportParameters(
        reportFile, referenceLocation, referenceSize, rate,
        totalInclusion, totalExclusion, inhits, exhits,
        k, kmerLocation, gap, size, GC):

    reportFile.write("==== Parameterization Report ====\n")
    reportFile.write("\n")
    reportFile.write("Reference File = " + str(referenceLocation) + "\n")
    reportFile.write("Reference Size = " + str(referenceSize) + "\n")
    reportFile.write("GC-Content = %.2f" % (GC) + "\n")
    reportFile.write("\n")
    reportFile.write("Mutation / Error Rate = " + str(rate) + "\n")
    reportFile.write("\n")
    reportFile.write("Inclusion Genomes = " + str(totalInclusion) + "\n")
    reportFile.write("Minimum Inclusion Hits = " + str(inhits) + "\n")
    reportFile.write("\n")
    reportFile.write("Exclusion Genomes = " + str(totalExclusion) + "\n")
    reportFile.write("Maximum Exclusion Hits = " + str(exhits) + "\n")
    reportFile.write("\n")
    reportFile.write("k-mer Size = " + str(k) + "\n")
    reportFile.write("k-mer File = " + str(kmerLocation) + "\n")
    reportFile.write("\n")
    reportFile.write("Maximum k-mer Gap Size = " + str(gap) + "\n")
    reportFile.write("Minimum Signature Size = " + str(size) + "\n")

"""
# =============================================================================

MAIN

# =============================================================================
"""
def main():

    # --- Parser ---
    parser = argparse.ArgumentParser(
        description='Extracts signatures from targets \
        using k-mer information.')

    parser.add_argument(
        REFERENCE_SHORT,
        REFERENCE_LONG,
        dest=REFERENCE,
        help="FASTA reference from which to extract signatures",
        type=str, required=True)

    parser.add_argument(
        REFERENCE_SIZE_SHORT,
        REFERENCE_SIZE_LONG,
        dest=REFERENCE_SIZE,
        help="estimated total reference size",
        type=int, required=False)

    parser.add_argument(
        RATE_SHORT,
        RATE_LONG,
        dest=RATE,
        help="probability of a mutation or error at an arbitrary position",
        type=float, default=0.01)

    parser.add_argument(
        INCLUSION_SHORT,
        INCLUSION_LONG,
        dest=INCLUSION,
        help="inclusion genome(s)",
        type=str, required=True, nargs='+')

    parser.add_argument(
        INHITS_SHORT,
        INHITS_LONG,
        dest=INHITS,
        help="minimum inclusion hits to build candidate",
        type=int, required=False)

    parser.add_argument(
        EXCLUSION_SHORT,
        EXCLUSION_LONG,
        dest=EXCLUSION,
        help="exclusion genome(s)",
        type=str, required=True, nargs='+')

    parser.add_argument(
        EXHITS_SHORT,
        EXHITS_LONG,
        dest=EXHITS,
        help="minimum exclusion hits to remove candidate",
        type=int, required=False)

    parser.add_argument(
        KMERS_SHORT,
        KMERS_LONG,
        dest=KMERS,
        help="k-mer file",
        type=str, required=True)

    parser.add_argument(
        GAP_SHORT,
        GAP_LONG,
        dest=GAP,
        help="maximum number of consecutive k-mers in a candidate \
            without an inclusion hit",
        type=int, required=False)

    parser.add_argument(
        SIZE_SHORT,
        SIZE_LONG,
        dest=SIZE,
        help="minimum candidate size",
        type=int, required=False)

    parser.add_argument(
        GC_SHORT,
        GC_LONG,
        dest=GC_CONTENT,
        help="the GC-content of the environment",
        type=float, required=False)

    parser.add_argument(
        CONFIDENCE_SHORT,
        CONFIDENCE_LONG,
        dest=CONFIDENCE,
        help="statistical confidence level",
        type=float, required=False, default=0.95)

    parser.add_argument(
        OUTPUT_SHORT,
        OUTPUT_LONG,
        dest=OUTPUT,
        help="output file",
        type=str, required=True)

    args = parser.parse_args()

    # --- Reference ---
    print "\n==== Building References ====\n"

    if not os.path.isfile(args.reference):
            raise RuntimeError("ERROR: Could not open the reference file.\n")

    referenceLocation = args.reference
    referenceFile = open(referenceLocation, 'r')
    references = buildReferences(referenceFile)
    referenceFile.close()

    print "\nDone!"

    # --- Reference Size & GC-Content ---
    if not args.referenceSize or not args.gcContent:
        referenceSize, GC = estimateReferenceParameters(references)

    if args.referenceSize:
        referenceSize = args.referenceSize

    if args.gcContent:
        GC = args.gcContent

    # --- Rate ---
    rate = args.rate

    # --- Statistical Confidence ---
    confidence = args.confidence

    # --- k-mer Size ---
    if not os.path.isfile(args.kmers):
            raise RuntimeError("ERROR: Could not open k-mer file.\n")

    kmerLocation = args.kmers
    kmerFile = open(kmerLocation, 'r')
    k = estimateK(kmerFile)
    kmerFile.close()

    # --- Minimum Inclusion Hits ---
    totalInclusion = len(args.inclusion)

    if args.inhits:
        inhits = args.inhits

    else:
        inhits = estimateInclusionHits(totalInclusion, rate, GC, k, confidence)

    # --- Maximum Exclusion Hits ---
    totalExclusion = len(args.exclusion)

    if args.exhits:
        exhits = args.exhits

    else:
        exhits = estimateExclusionHits(totalExclusion, rate, k)

    # --- k-mer Tables ---
    print "\n==== Building k-mer Tables ====\n"

    kmerFile = open(args.kmers, 'r')
    inmers = {}
    exmers = {}
    buildKMers(kmerFile, inmers, exmers, inhits, exhits)
    kmerFile.close()

    print "Done!"

    # --- Gap Size ---
    if args.gap:
        gap = args.gap

    else:
        gap = estimateGapSize(rate, GC, k, confidence)

    # --- Minimum Signature Size ---
    if args.size:
        size = args.size

    else:
        size = estimateSignatureSize(k)

    # --- Report ---
    reportLocation = str(args.output) + ".report"
    reportFile = open(reportLocation, 'w')
    reportParameters(
        reportFile, referenceLocation, referenceSize, rate,
        totalInclusion, totalExclusion, inhits, exhits,
        k, kmerLocation, gap, size, GC)
    reportFile.close()

    # --- Output File ---
    print "\n==== Extracing Signatures ====\n"
    outputFile = open(args.output, 'w')
    extract(references, k, inmers, exmers, size, gap, outputFile)
    outputFile.close()

    print "Done!"

    print "\n==== Exiting ====\n"

if __name__ == '__main__':

    main()
