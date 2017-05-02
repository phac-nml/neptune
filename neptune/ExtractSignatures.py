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

EXTRACT SIGNATURES
------------------


PURPOSE
-------

This script extracts candidate signatures from one or more inclusion files
within the context of one or more exclusion files. The signatures are extracted
from a reference using k-mers generated from all inclusion and exclusion inputs.


INPUT
-----

(aggregated k-mers)

[k-mer] [inclusion counts] [exclusion counts]

AAAAA 3 1
AAAAC 1 2
AAAAG 3 3
AAAAT 3 0


USAGE
-----

ExtractSignatures.py [-h] -r REFERENCE -i INCLUSION [INCLUSION ...] -e
                            EXCLUSION [EXCLUSION ...] -k KMERS -o OUTPUT
                            [-rs REFERENCE-SIZE] [-q RATE] [-ih INHITS]
                            [-eh EXHITS] [-g GAP] [-s SIZE] [-gc GC-CONTENT]
                            [-c CONFIDENCE]


EXAMPLE
-------

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

import Signature

from scipy.stats import norm

"""
# =============================================================================

GLOBALS

# =============================================================================
"""

PROGRAM_DESCRIPTION = 'This script extracts signatures from reference targets \
    using k-mer information.'

# DEFAULTS #

RATE_DEFAULT = 0.01
CONFIDENCE_DEFAULT = 0.95

# ARGUMENTS #

LONG = "--"
SHORT = "-"

# REQUIRED ARGUMENTS #

# Reference
REFERENCE = "reference"
REFERENCE_LONG = LONG + REFERENCE
REFERENCE_SHORT = SHORT + "r"
REFERENCE_HELP = "The FASTA reference from which to extract signatures."

# Inclusion Targets
INCLUSION = "inclusion"
INCLUSION_LONG = LONG + INCLUSION
INCLUSION_SHORT = SHORT + "i"
INCLUSION_HELP = "The inclusion targets in FASTA format."

# Exclusion Targets
EXCLUSION = "exclusion"
EXCLUSION_LONG = LONG + EXCLUSION
EXCLUSION_SHORT = SHORT + "e"
EXCLUSION_HELP = "The exclusion targets in FASTA format."

# Aggregated k-mers File
KMERS = "kmers"
KMERS_LONG = LONG + KMERS
KMERS_SHORT = SHORT + "k"
KMERS_HELP = "The aggregated k-mer file produced by AggregateKMers.py."

# Output File
OUTPUT = "output"
OUTPUT_LONG = LONG + OUTPUT
OUTPUT_SHORT = SHORT + "o"
OUTPUT_HELP = "The location to output candidate signatures in FASTA format."

# OPTIONAL ARGUMENTS #

# Reference Size
REFERENCE_SIZE = "reference-size"
REFERENCE_SIZE_LONG = LONG + REFERENCE_SIZE
REFERENCE_SIZE_SHORT = SHORT + "rs"
REFERENCE_SIZE_HELP = "The estimated total size in nucleotides of the \
    reference. This will be calculated if not specified."

# Rate (of Errors, Mutations)
RATE = "rate"
RATE_LONG = LONG + RATE
RATE_SHORT = SHORT + "q"
RATE_HELP = "The probability of a mutation or error at an arbitrary position. \
    The default value is " + str(RATE_DEFAULT) + "."

# Minimum Inclusion Hits
INHITS = "inhits"
INHITS_LONG = LONG + INHITS
INHITS_SHORT = SHORT + "ih"
INHITS_HELP = "The minimum number of inclusion targets that must contain a \
    k-mer observed in the reference to begin or continue building candidate \
    signatures. This will be calculated if not specified."

# Maximum Exclusion Hits
EXHITS = "exhits"
EXHITS_LONG = LONG + EXHITS
EXHITS_SHORT = SHORT + "eh"
EXHITS_HELP = "The maximum allowable number of exclusion targets that may \
    contain a k-mer observed in the reference before terminating the \
    construction of a candidate signature. This will be calculated if \
    not specified."

# k-mer Gap Size
GAP = "gap"
GAP_LONG = LONG + GAP
GAP_SHORT = SHORT + "g"
GAP_HELP = "The maximum number of consecutive k-mers observed in the \
    reference during signature candidate construction that fail to have \
    enough inclusion hits before terminating the construction of a candidate \
    signature. This will be calculated if not specified and is determined \
    from the size of k and the rate."

# Minimum Signature Size
SIZE = "size"
SIZE_LONG = LONG + SIZE
SIZE_SHORT = SHORT + "s"
SIZE_HELP = "The minimum size of all reported candidate signatures. \
    Identified candidate signatures shorter than this value will be discard."

# GC Content
GC_CONTENT = "gc-content"
GC_LONG = LONG + GC_CONTENT
GC_SHORT = SHORT + "gc"
GC_HELP = "The average GC-content of all inclusion and exclusion targets. \
    This will be calculated from inclusion and exclusion targets if not \
    specified."

# Statistical Confidence
CONFIDENCE = "confidence"
CONFIDENCE_LONG = LONG + CONFIDENCE
CONFIDENCE_SHORT = SHORT + "c"
CONFIDENCE_HELP = "The statistical confidence level in decision making \
    involving probabilities when producing candidate signatures."

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
-------


PURPOSE
-------

Extracts candidate signatures from the set of references.

The extraction process is guided by the inclusion and exclusion k-mer
dictionaries and the reference. This function is concerned only about k-mer
existence in those dictionaries.


INPUT
-----

[STRING ITERABLE] [references]
    An iterable object of string references. This is intended to correspond to
    a list of single or multi-fasta files.

[INT >= 1] [k]
    The k-mer size.

[KMER DICTIONARY] [inmers]
    The inclusion k-mers dictionary.

[KMER DICTIONARY] [exmers]
    The exclusion k-mers dictionary.

[INT >= 1] [size]
    The minimum signature size in characters.

[INT >= 1] [gap]
    The maximum allowable gap size in k-mers. Note that this works in k-mer
    space, not in sequence character spaces.

[FILE] [outputFile]
    The output file to write candidate signatures.


POST
----

The candidate signatures will be written to the [outputFile].

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

        signature = Signature.Signature(
            i, 0.0, 0.0, 0.0, regions[i].sequence,
            regions[i].reference, regions[i].position)

        Signature.writeSignature(signature, outputFile)


"""
# =============================================================================

CALCULATE PROBABILITY THAT HOMOLOGOUS BASES MUTATE AND MATCH
------------------------------------------------------------

P(X_M = Y_M)_H


PURPOSE
-------

Calculates the probability that homologous bases both mutate to the same
matching base, given a certain GC content.


INPUT
-----

[0 <= FLOAT <= 1] [GC]
    The GC-content of the environment. A high GC means that it is more likely
    that homologuous bases both mutate to the same base, because most mutations
    will be to either G's or C's.


RETURN
------

[0 <= FLOAT <= 1] [P(X_M = Y_M)_H]
    The proability that homologous bases both mutate and match each other.

# =============================================================================
"""
def calculateProbHBMM(GC):

    # 0 <= GC <= 1
    if GC < 0 or GC > 1:
        raise RuntimeError("The GC-content is out of range.")

    GC = float(GC)

    # P(X_M = Y_M)_H
    probHBMM = (
        2 * math.pow(GC / (GC + 1), 2) +
        math.pow((1 - GC) / (GC + 1), 2)) * (1 - GC) +\
        (2 * math.pow((1 - GC) / (2 - GC), 2) +
            math.pow(GC / (2 - GC), 2)) * (GC)

    return probHBMM


"""
# =============================================================================

CALCULATE PROBABILITY THAT HOMOLOGOUS BASES MATCH
-------------------------------------------------

P(X = Y)_H


PURPOSE
-------

Calculates the probability that homologous bases match, given a certain
mutation rate and a GC-content environment.


INPUT
-----

[0 <= FLOAT <= 1] [mutationRate]
    The probability of an arbitrary base mutating.

[0 <= FLOAT <= 1] [GC]
    The GC content of the environment (or genomes).


RETURN
------

[0 <= FLOAT <= 1] [P(X = Y)_H]
    The probability that homologous bases match.

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
--------------------------------------------------

P(k_X = k_Y)_H


PURPOSE
-------

Calculates the probability that homologous k-mers match, given a certain
mutation rate, GC-content, and k-mer size.


INPUT
-----

[0 <= FLOAT <= 1] [mutationRate]
    The probability of an arbitrary base mutating.

[0 <= FLOAT <= 1] [GC]
    The GC content of the environment.

[INT >= 1] [kmerSize]
    The size of the k-mers.


RETURN
------

[0 <= FLOAT <= 1] [P(k_X = k_Y)_H]
    The probability that homologous k-mers match.

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
-----------------------


PURPOSE
-------

Estimates the minimum candidate size. This function uses a constant factor and
multiplies it by the k-mer size to determine a minimum default signature size.

A minimum signature size is often useful because longer signatures tend to be
more specific that shorter signatures.

INPUT
-----

[INT >= 1] [kmerSize]
    The size of the k-mers.


RETURN
------

[INT] [estimate]
    An estimate of the minimum candidate signature size.

# =============================================================================
"""
def estimateSignatureSize(kmerSize):

    FACTOR = 4.0

    estimate = FACTOR * kmerSize

    return int(estimate)


"""
# =============================================================================

ESTIMATE GAP SIZE
-----------------


PURPOSE
-------

Estimates the maximum k-mer gap size before abandoning a candidate signature
region. When k-mer gap sizes become too long, it becomes more likely that we
are observing multiple signatures, rather than one signature.

INPUT
-----

[0 <= FLOAT <= 1] [mutationRate]
    The probability of an arbitrary base mutating.

[0 <= FLOAT <= 1] [GC]
    The GC-content of the environment.

[INT >= 1] [kmerSize]
    The size of the k-mers.

[0 < FLOAT < 1] [confidence]
    The statistical confidence level in decision making involving
    probabilities when producing candidate signatures.


RETURN
------

[INT] [estimate]
    An estimate of the maximum allowable k-mer gap size. This will be rounded
    up to the next integer.

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
-----------------------


PURPOSE
-------

Estimates the minimum number of exclusion hits to prevent candidate building.
This function currently defaults to returning 1, meaning any exclusion hit
will immediately end signature construction. This approach achieves maximum
k-mer specificity. However, the sequence specificity is usually not 100%.

This function is left open to accept several parameters which might be used
in the future to determine a more nuanced estimate.


INPUT
-----

[INT >= 0] [totalExclusion]
    The total number of exclusion targets.

[0 <= FLOAT <= 1] [rate]
    The probability of an arbitrary base mismatching (ex: SNV).

[INT >= 1] [kmerSize]
    The size of the k-mers.


RETURN
------

[INT >= 1] [estimate]
    An estimate of the minimum number of exclusion k-mer hits before rejection.

# =============================================================================
"""
def estimateExclusionHits(totalExclusion, rate, kmerSize):

    ESTIMATE = int(1)

    return ESTIMATE


"""
# =============================================================================

ESTIMATE INCLUSION HITS
-----------------------


PURPOSE
-------

Estimates the minimum number of inclusion hits required for confident
candidate signature building.


INPUT
-----

[INT >= 0] [totalInclusion]
    The total number of inclusion targets.

[0 <= FLOAT <= 1] [mutationRate]
    The probability of an arbitrary base mutating.

[0 <= FLOAT <= 1] [GC]
    The GC content of the environment.

[INT >= 1] [kmerSize]
    The size of the k-mers.

[0 < FLOAT < 1] [confidence]
    The statistical confidence level in decision making involving
    probabilities when producing candidate signatures.


RETURN
------

[INT >= 0] [estimate]
    An estimate of the minimum number of inclusion k-mer hits. This is the
    number of inclusion targets that must share a k-mer observed in a reference
    for it to be considered in candidate signature construction.

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
----------


PURPOSE
-------

Estimates the size of k from the k-mers file. This function will read the first
k-mer in the file and assume the k-mers are all the same length.


INPUT
-----

[FILE] [kmerFile]
    The file of aggregated k-mers. This file must be open and ready to be
    read from the start of the file.


RETURN
------

[INT >= 1] [estimate]
    An estimate of the k-mer size.


POST
----

One line will have been read from the file and the file will not be closed.

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
------------


PURPOSE
-------

Builds the inclusion and exclusion k-mer dictionaries from a single aggregated
k-mer file. This will likely cause a significant amount of memory to be
allocated.


INPUT
-----

[FILE] [kmerFile]
    A readable file-like object of aggregated k-mers.

[(STRING KMER) -> (INT) DICTIONARY] [inmers]
    The inclusion k-mer dictionary to fill with k-mers.

[(STRING KMER) -> (INT) DICTIONARY] [exmers]
    The exclusion k-mer dictionary to fill with k-mers.

[INT >= 0] [inhits]
    The minimum number of inclusion targets that must contain a k-mer observed
    in the reference to begin or continue building candidate signatures.

[INT >= 0] [exhits]
    The maximum allowable number of exclusion targets that may contain a k-mer
    observed in the reference before terminating the construction of a
    candidate signature.

POST:
    The inclusion and exclusion k-mer dictionaries will be filled with all
    k-mers found in the k-mers file with at least [inhits] and at least
    [exhits] counts for the inclusion and exclusion k-mers, respectively.

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
-----------------


PURPOSE
-------

This function outputs the parameters to standard output. This can be useful
for debugging.


INPUT
-----

[FILE] [reportFile]
    The writable file-like object to write the report.

[FILE LOCATION] [referenceLocation]
    The single reference to extract candidates from.

[INT >= 0] [referenceSize]
    The size of the reference.

[0 <= FLOAT <= 1] [rate]
    The rate of mutations and/or errors.

[INT >= 0] [totalInclusion]
    The number of inclusion genome files.

[INT >= 0] [totalExclusion]
    The number of exclusion genome files.

[INT >= 0] [inhits]
    The minimum number of inclusion k-mer hits.

[INT >= 0] [exhits]
    The maximum number of exclusion k-mer hits.

[INT >= 1] [k]
    The size of the k-mer.

[FILE LOCATION] [kmerLocation]
    The file containing aggregated k-mers.

[INT >= 1] [gap]
    The maximum inclusion k-mer gap size.

[INT >= 1] [size]
    The minimum size of any candidate.

[0 <= FLOAT <= 1] [GC]
    The average GC content of all the targets..


POST
----

The parameters will be written to [reportFile].

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
    reportFile.write("SNV Rate = " + str(rate) + "\n")
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

PARSE

# =============================================================================
"""
def parse(parameters):

    # --- Reference ---
    if not os.path.isfile(parameters[REFERENCE]):
        raise RuntimeError("ERROR: Could not open the reference file.\n")

    referenceLocation = parameters[REFERENCE]
    referenceFile = open(referenceLocation, 'r')
    references = buildReferences(referenceFile)
    referenceFile.close()

    # --- Reference Size & GC-Content ---
    if not parameters[REFERENCE_SIZE] or not parameters[GC_CONTENT]:
        referenceSize, GC = estimateReferenceParameters(references)

    if parameters[REFERENCE_SIZE]:
        referenceSize = parameters[REFERENCE_SIZE]

    if parameters[GC_CONTENT]:
        GC = parameters[GC_CONTENT]

    # --- Rate ---
    rate = parameters.get(RATE) \
        if parameters.get(RATE) else RATE_DEFAULT

    # --- Statistical Confidence ---
    confidence = parameters.get(CONFIDENCE) \
        if parameters.get(CONFIDENCE) else CONFIDENCE_DEFAULT

    # --- k-mer Size ---
    if not os.path.isfile(parameters[KMERS]):
            raise RuntimeError("ERROR: Could not open k-mer file.\n")

    kmerLocation = parameters[KMERS]
    kmerFile = open(kmerLocation, 'r')
    k = estimateK(kmerFile)
    kmerFile.close()

    # --- Minimum Inclusion Hits ---
    totalInclusion = len(parameters[INCLUSION])

    if parameters[INHITS]:
        inhits = parameters[INHITS]

    else:
        inhits = estimateInclusionHits(totalInclusion, rate, GC, k, confidence)

    # --- Maximum Exclusion Hits ---
    totalExclusion = len(parameters[EXCLUSION])

    if parameters[EXHITS]:
        exhits = parameters[EXHITS]

    else:
        exhits = estimateExclusionHits(totalExclusion, rate, k)

    # --- k-mer Tables ---
    kmerFile = open(parameters[KMERS], 'r')
    inmers = {}
    exmers = {}
    buildKMers(kmerFile, inmers, exmers, inhits, exhits)
    kmerFile.close()

    # --- Gap Size ---
    if parameters[GAP]:
        gap = parameters[GAP]

    else:
        gap = estimateGapSize(rate, GC, k, confidence)

    # --- Minimum Signature Size ---
    if parameters[SIZE]:
        size = parameters[SIZE]

    else:
        size = estimateSignatureSize(k)

    # --- Report ---
    reportLocation = str(parameters[OUTPUT]) + ".report"
    reportFile = open(reportLocation, 'w')
    reportParameters(
        reportFile, referenceLocation, referenceSize, rate,
        totalInclusion, totalExclusion, inhits, exhits,
        k, kmerLocation, gap, size, GC)
    reportFile.close()

    # --- Extraction ---
    outputFile = open(parameters[OUTPUT], 'w')
    extract(references, k, inmers, exmers, size, gap, outputFile)
    outputFile.close()


"""
# =============================================================================

MAIN

# =============================================================================
"""
def main():

    parser = argparse.ArgumentParser(description=PROGRAM_DESCRIPTION)

    # REQUIRED #

    parser.add_argument(
        REFERENCE_SHORT,
        REFERENCE_LONG,
        dest=REFERENCE,
        help=REFERENCE_HELP,
        type=str, required=True)

    parser.add_argument(
        INCLUSION_SHORT,
        INCLUSION_LONG,
        dest=INCLUSION,
        help=INCLUSION_HELP,
        type=str, required=True, nargs='+')

    parser.add_argument(
        EXCLUSION_SHORT,
        EXCLUSION_LONG,
        dest=EXCLUSION,
        help=EXCLUSION_HELP,
        type=str, required=True, nargs='+')

    parser.add_argument(
        KMERS_SHORT,
        KMERS_LONG,
        dest=KMERS,
        help=KMERS_HELP,
        type=str, required=True)

    parser.add_argument(
        OUTPUT_SHORT,
        OUTPUT_LONG,
        dest=OUTPUT,
        help=OUTPUT_HELP,
        type=str, required=True)

    # OPTIONAL #

    parser.add_argument(
        REFERENCE_SIZE_SHORT,
        REFERENCE_SIZE_LONG,
        dest=REFERENCE_SIZE,
        help=REFERENCE_SIZE_HELP,
        type=int, required=False)

    parser.add_argument(
        RATE_SHORT,
        RATE_LONG,
        dest=RATE,
        help=RATE_HELP,
        type=float, required=False)

    parser.add_argument(
        INHITS_SHORT,
        INHITS_LONG,
        dest=INHITS,
        help=INHITS_HELP,
        type=int, required=False)

    parser.add_argument(
        EXHITS_SHORT,
        EXHITS_LONG,
        dest=EXHITS,
        help=EXHITS_HELP,
        type=int, required=False)

    parser.add_argument(
        GAP_SHORT,
        GAP_LONG,
        dest=GAP,
        help=GAP_HELP,
        type=int, required=False)

    parser.add_argument(
        SIZE_SHORT,
        SIZE_LONG,
        dest=SIZE,
        help=SIZE_HELP,
        type=int, required=False)

    parser.add_argument(
        GC_SHORT,
        GC_LONG,
        dest=GC_CONTENT,
        help=GC_HELP,
        type=float, required=False)

    parser.add_argument(
        CONFIDENCE_SHORT,
        CONFIDENCE_LONG,
        dest=CONFIDENCE,
        help=CONFIDENCE_HELP,
        type=float, required=False)

    args = parser.parse_args()
    parameters = vars(args)
    parse(parameters)


"""
# =============================================================================
# =============================================================================
"""
if __name__ == '__main__':

    main()
