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
Date: 13 May 2015

INPUT:

> 1
AACCGGAT
> 2
GGAGAATTCT
> 3
TTTAGAGC

OUTPUT:

> 2
GGAGAATTCT

USAGE:

FilterSignatures.py -h
FilterSignatures.py -d DATABASE -i INPUT -o OUTPUT
    [-fp FILTER_PERCENT] [-fl FILTER_LENGTH]

EXAMPLE:

script.py -d exclusion.db -i inclusion1.candidates -o inclusion1.filtered
script.py -d ex.db -i in1.can -o in1.fil -fp 0.70 -fl 0.50

# =============================================================================
"""

import argparse
import subprocess
import operator

"""
# =============================================================================

GLOBALS

# =============================================================================
"""

# NAMES
INCLUSION_DATABASE = "inclusionDatabase"
EXCLUSION_DATABASE = "exclusionDatabase"
INCLUSION = "inclusion"
EXCLUSION = "exclusion"
INPUT = "input"
FILTERED_OUTPUT = "filteredOutput"
SORTED_OUTPUT = "sortedOutput"
FILTER_PERCENT = "filterPercent"
FILTER_LENGTH = "filterLength"

# ARGUMENTS
LONG = "--"

INCLUSION_DATABASE_LONG = LONG + INCLUSION_DATABASE
EXCLUSION_DATABASE_LONG = LONG + EXCLUSION_DATABASE
INCLUSION_LONG = LONG + INCLUSION
EXCLUSION_LONG = LONG + EXCLUSION
INPUT_LONG = LONG + INPUT
FILTERED_OUTPUT_LONG = LONG + FILTERED_OUTPUT
SORTED_OUTPUT_LONG = LONG + SORTED_OUTPUT
FILTER_PERCENT_LONG = LONG + FILTER_PERCENT
FILTER_LENGTH_LONG = LONG + FILTER_LENGTH

SHORT = "-"

INCLUSION_DATABASE_SHORT = SHORT + "dbin"
EXCLUSION_DATABASE_SHORT = SHORT + "dbex"
INCLUSION_SHORT = SHORT + "i"
EXCLUSION_SHORT = SHORT + "e"
INPUT_SHORT = SHORT + "r"
FILTERED_OUTPUT_SHORT = SHORT + "fo"
SORTED_OUTPUT_SHORT = SHORT + "so"
FILTER_PERCENT_SHORT = LONG + "fp"
FILTER_LENGTH_SHORT = SHORT + "fl"

overallScore = {}
inclusionScore = {}
exclusionScore = {}

"""
# =============================================================================

HIT

# =============================================================================
"""
class Hit():

    def __init__(self, line):

        tokens = line.split()

        ID = tokens[0]
        length = tokens[1]
        reference = tokens[2]
        alignmentLength = tokens[3]
        percentIdentity = tokens[4]
        alignmentScore = tokens[5]

        self.ID = str(ID).strip()
        self.length = int(length)
        self.reference = str(reference).strip()
        self.alignmentLength = int(alignmentLength)
        self.percentIdentity = float(percentIdentity)
        self.alignmentScore = float(alignmentScore)

        self.neptuneScore = float(0.0)

"""
# =============================================================================

SIGNATURE

# =============================================================================
"""
class Signature():

    def __init__(self, ID, sequence, reference, position):

        self.ID = ID
        self.sequence = sequence.strip()
        self.length = len(sequence.strip())
        self.reference = reference
        self.position = position

    """
    # =========================================================================

    READ SIGNATURES

    PURPOSE:
        Reads a signature file into a signature dictionary.

    INPUT:
        [STRING] [fileLocation] - The file location of the signatures.

    RETURN:
        [(STRING) -> (SIGNATURE) DICTIONARY] - A dictionary mapping string IDs
            to signature objects.

    # =========================================================================
    """
    @staticmethod
    def readSignatures(fileLocation):

        signaturesFile = open(fileLocation, 'r')
        signatures = {}

        while True:

            line1 = signaturesFile.readline()
            line2 = signaturesFile.readline()

            if not line2:
                break

            tokens = (line1[1:]).split()
            ID = tokens[0]
            reference = tokens[2]
            position = tokens[3]

            sequence = line2

            signatures[ID] = Signature(ID, sequence, reference, position)

        signaturesFile.close()

        return signatures

    """
    # =========================================================================

    WRITE SIGNATURES

    PURPOSE:
        Writes the signature to the destination. This function is designed to
        be symmetric with the read signatures function.

    INPUT:
        [SIGNATURE] [signature] - The signature to write.
        [WRITABLE] [destination] - An open and writable object.

    RETURN:
        [NONE]

    POST:
        The passed signature will be written to the destination.

    # =========================================================================
    """
    @staticmethod
    def writeSignature(signature, destination):

        destination.write(
            ">" + str(signature.ID) + " " + str(signature.length)
            + " " + str(signature.reference) + " "
            + str(signature.position) + "\n")

        destination.write(str(signature.sequence) + "\n")

"""
# =============================================================================

QUERY DATABASE

PURPOSE:
    Queries the database with a specified query.

INPUT:
    [STRING] [databaseLocation] - The file location of the database.
    [STRING] [queryLocation] - The file location of the query.
    [STRING] [outputLocation] - The file location of the output.
    [0 <= FLOAT <= 1] [filterPercent] - The maximum percent identity of an
        exclusion hit with a candidate.

POST:
    A query file will be created in the output directory.

RETURN:
    [STRING] [queryOutputLocation] - The file location of the query output.

# =============================================================================
"""
def queryDatabase(
        databaseLocation, queryLocation, outputLocation, filterPercent):

    # COMMAND LINE
    COMMAND = "blastn"

    DATABASE = "-db"
    QUERY = "-query"
    OUTPUT = "-out"
    OUTPUT_FORMAT = "-outfmt"
    OUTPUT_FORMAT_STRING = "6 qseqid qlen sseqid length pident score"
    PERCENT_IDENTITY = "-perc_identity"
    WORD_SIZE = "-word_size"
    WORD_SIZE_VALUE = 11
    DUST = "-dust"
    DUST_VALUE = "no"

    queryOutputLocation = outputLocation + ".query"

    args = [
        COMMAND,
        DATABASE, databaseLocation,
        QUERY, queryLocation,
        OUTPUT, queryOutputLocation,
        OUTPUT_FORMAT, OUTPUT_FORMAT_STRING,
        PERCENT_IDENTITY, str(filterPercent),
        WORD_SIZE, str(WORD_SIZE_VALUE),
        DUST, DUST_VALUE]

    subprocess.check_call(args)

    return queryOutputLocation

"""
# =============================================================================

UPDATE HIT OVERALL DICTIONARY

PURPOSE:
    Updates the passed dictionary with the passed hit. This function is looking
    for the best hit for a query, regardless of the subject (reference).

INPUT:
    [HIT] [hit] - The hit object associated with the hit.
    [(STRING) -> (HIT) DICTIONARY] - The best overall query hit dictionary.
        This dictionary maps signature (query) IDs to hit objects.

RETURN:
    [NONE]

POST:
    The passed dictionary object will be updated.

# =============================================================================
"""
def updateHitOverallDictionary(hit, dictionary):

    key = hit.ID

    if key not in dictionary:

        dictionary[key] = hit

    elif hit.alignmentScore > dictionary[key].alignmentScore:

        dictionary[key] = hit

"""
# =============================================================================

UPDATE HIT PAIR DICTIONARY

PURPOSE:
    Updates the passed dictionary with the passed hit. This function is looking
    for the best hit for a (query, reference) pair.

INPUT:
    [HIT] [hit] - The hit object associated with the hit.
    [(STRING, STRING) -> (HIT) DICTIONARY] - The best overall query hit
        dictionary. This dictionary maps (query, reference) tuples to hit
        objects.

RETURN:
    [NONE]

POST:
    The passed dictionary object will be updated.

# =============================================================================
"""
def updateHitPairDictionary(hit, dictionary):

    key = (hit.ID, hit.reference)

    if key not in dictionary:

        hit.neptuneScore = (
            (float(hit.alignmentLength) / float(hit.length))
            * (float(hit.percentIdentity) / float(100)))

        dictionary[key] = hit

    elif hit.alignmentScore > dictionary[key].alignmentScore:

        hit.neptuneScore = (
            (float(hit.alignmentLength) / float(hit.length))
            * (float(hit.percentIdentity) / float(100)))

        dictionary[key] = hit

"""
# =============================================================================

UPDATE EXCLUSION SCORES

PURPOSE:
    Updates the exclusion score and negative component of the overall score.

INPUT:
    [(STRING, STRING) -> (HIT) DICTIONARY] - The best overall query hit
        dictionary. This dictionary maps (query, reference) tuples to hit
        objects.
    [INT >= 1] [totalExclusion] - The total number of exclusion files.

RETURN:
    [NONE]

POST:
    The exclusion scores and negative component of the overall scores will be
    updated.

# =============================================================================
"""
def updateExclusionScores(dictionary, totalExclusion):

    for key in dictionary:

        hit = dictionary[key]
        ID = hit.ID
        neptuneScore = hit.neptuneScore

        if ID not in overallScore:
            overallScore[ID] = (
                -float(neptuneScore) / float(totalExclusion))

        else:
            overallScore[ID] -= (
                float(neptuneScore) / float(totalExclusion))

        if ID not in exclusionScore:
            exclusionScore[ID] = (
                -float(neptuneScore) / float(totalExclusion))

        else:
            exclusionScore[ID] -= (
                float(neptuneScore) / float(totalExclusion))

"""
# =============================================================================

UPDATE INCLUSION SCORES

PURPOSE:
    Updates the inclusion score and positive component of the overall score.

INPUT:
    [(STRING, STRING) -> (HIT) DICTIONARY] - The best overall query hit
        dictionary. This dictionary maps (query, reference) tuples to hit
        objects.
    [INT >= 1] [totalInclusion] - The total number of inclusion files.

RETURN:
    [NONE]

POST:
    The inclusion scores and positive component of the overall scores will be
    updated.

# =============================================================================
"""
def updateInclusionScores(dictionary, totalInclusion):

    for key in dictionary:

        hit = dictionary[key]
        ID = hit.ID
        neptuneScore = hit.neptuneScore

        if ID not in overallScore:
            overallScore[ID] = (
                float(neptuneScore) / float(totalInclusion))

        else:
            overallScore[ID] += (
                float(neptuneScore) / float(totalInclusion))

        if ID not in inclusionScore:
            inclusionScore[ID] = (
                float(neptuneScore) / float(totalInclusion))

        else:
            inclusionScore[ID] += (
                float(neptuneScore) / float(totalInclusion))

"""
# =============================================================================

# =============================================================================
"""
def reportCandidates(
        candidatesLocation, outputLocation,
        hitOverallDictionary, filterLength):

    outputFile = open(outputLocation, 'w')
    candidateSignatures = Signature.readSignatures(candidatesLocation)

    for ID in candidateSignatures:

        signature = candidateSignatures[ID]

        if ID in hitOverallDictionary:

            hit = hitOverallDictionary[ID]

            if (float(hit.alignmentLength) / float(signature.length)
                    < float(filterLength)):
                Signature.writeSignature(signature, outputFile)

        else:
            Signature.writeSignature(signature, outputFile)

    outputFile.close()

"""
# =============================================================================

# =============================================================================
"""
def reportSorted(filteredLocation, outputLocation, sortedSignatureIDs):

    filteredSignatures = Signature.readSignatures(filteredLocation)

    # REPORT SORTED SIGNATURES
    outputFile = open(outputLocation, 'w')

    for ID in sortedSignatureIDs:

        if ID in filteredSignatures:

            signature = filteredSignatures[ID]

            outputFile.write(">" + str(signature.ID))

            if ID in overallScore:
                outputFile.write(
                    " score=" + "{0:.2f}".format(
                        round(overallScore[ID], 2)))
            else:
                outputFile.write(" score=0.00")

            if ID in inclusionScore:
                outputFile.write(
                    " in=" + "{0:.2f}".format(
                        round(inclusionScore[ID], 2)))
            else:
                outputFile.write(" in=0.00")

            if ID in exclusionScore:
                outputFile.write(
                    " ex=" + "{0:.2f}".format(
                        round(exclusionScore[ID], 2)))
            else:
                outputFile.write(" ex=0.00")

            outputFile.write(" len=" + str(signature.length))
            outputFile.write(" ref=" + str(signature.reference))
            outputFile.write(" pos=" + str(signature.position))
            outputFile.write("\n")

            outputFile.write(signature.sequence + "\n")

    outputFile.close()

"""
# =============================================================================

REPORT SIGNATURES

PURPOSE:
    Reports the candidate signatures which are not found in the list of
    filterable signatures.

INPUT:
    [STRING] [candidatesLocation] - The file location of the candidate
        signatures.
    [STRING] [exclusionQueryLocation] - The file location of the signaturess
        to filter.
    [STRING] [outputLocation] - The location for the output of the filtered
        signatures.
    [0 <= FLOAT 0 <= 1] [filterLength] - The maximum percent length of an
        exclusion hit with a candidate.
    [1 <= INT] [totalExclusion] - The total number of exclusion targets.

POST:
    A file of filterted signatures will be produced.

RETURN:
    [STRING] [outputLocation] - The file location of the output.

# =============================================================================
"""
def reportSignatures(
        candidateLocation, exclusionQueryLocation, outputLocation,
        filterLength, totalExclusion):

    exclusionQueryFile = open(exclusionQueryLocation, 'r')

    hitOverallDictionary = {}
    hitPairDictionary = {}

    # LOAD FILTER
    for line in exclusionQueryFile:

        hit = Hit(line)
        updateHitOverallDictionary(hit, hitOverallDictionary)
        updateHitPairDictionary(hit, hitPairDictionary)

    exclusionQueryFile.close()

    updateExclusionScores(hitPairDictionary, totalExclusion)

    reportCandidates(
        candidateLocation, outputLocation, hitOverallDictionary, filterLength)

    return outputLocation

"""
# =============================================================================

SORT

PURPOSE:
    Sorts the filtered signatures according to their signature score.

INPUT:
    [STRING] [filteredLocation] - The file location of the filtered signatures.
    [STRING] [inclusionQueryLocation] - The file location of the filtered
        signatures against the inclusion targets.
    [STRING] [outputLocation] - The file output location of the sorted
        signatures.
    [0 <= FLOAT 0 <= 1] [filterLength] - The minimum query alignment length
        for the signature to be considered a hit and used in scoring.
    [0 <= INT] [totalInclusion] - The total number of inclusion targets.

POST:
    The signatures will be written to the [outputLocation] in score-descending
    order.

RETURN:
    [STRING] [outputLocation] - The location of the output file.

# =============================================================================
"""
def sortSignatures(
        filteredLocation, inclusionQueryLocation,
        outputLocation, totalInclusion):

    inclusionQueryFile = open(inclusionQueryLocation, 'r')

    hitPairDictionary = {}

    # LOAD FILTER
    for line in inclusionQueryFile:

        hit = Hit(line)
        updateHitPairDictionary(hit, hitPairDictionary)

    inclusionQueryFile.close()

    updateInclusionScores(hitPairDictionary, totalInclusion)

    sortedSignatureIDs = [ID for (ID, score) in sorted(
        overallScore.items(), key=operator.itemgetter(1), reverse=True)]

    reportSorted(filteredLocation, outputLocation, sortedSignatureIDs)

    return outputLocation

"""
# =============================================================================

FILTER SIGNATURES

PURPOSE:
    Filters and prepares the candidate signatures for output to the user.

INPUT:
    [STRING] [inclusionDatabaseLocation] - The file location of the
        inclusion database.
    [STRING] [exclusionDatabaseLocation] - The file location of the
        exclusion database.
    [0 <= INT] [totalInclusion] - The total number of inclusion targets.
    [0 <= INT] [totalExclusion] - The total number of exclusion targets.
    [STRING] [candidateLocation] - The location of the candidate signatures.
    [STRING] [filteredOutputLocation] - The file location to write the
        filtered signatures.
    [STRING] [sortedOutputLocation] - The file location to write the sorted
        signatures.
    [0 <= FLOAT 0 <= 1] [filterLength] - The minimum query alignment length
        for the signature to be considered a hit and used in scoring.
    [0 <= FLOAT 0 <= 1] [filterLength] - The maximum percent length of an
        exclusion hit with a candidate.
POST:
    Filtered signatures will be written to [filteredOutputLocation] and
    sorted signatures will be written to [sortedOutputLocation].

RETURN:
    NONE

# =============================================================================
"""
def filterSignatures(
        inclusionDatabaseLocation, exclusionDatabaseLocation,
        totalInclusion, totalExclusion, candidateLocation,
        filteredOutputLocation, sortedOutputLocation, filterLength,
        filterPercent):

    # BLASTN - EXCLUSION
    exclusionQueryLocation = queryDatabase(
        exclusionDatabaseLocation, candidateLocation, filteredOutputLocation,
        filterPercent)

    # FILTER
    filteredLocation = reportSignatures(
        candidateLocation, exclusionQueryLocation, filteredOutputLocation,
        filterLength, totalExclusion)

    # BLASTN - INCLUSION
    inclusionQueryLocation = queryDatabase(
        inclusionDatabaseLocation, filteredLocation,
        sortedOutputLocation, filterPercent)

    # SORT
    sortSignatures(
        filteredLocation, inclusionQueryLocation,
        sortedOutputLocation, totalInclusion)

"""
# =============================================================================

MAIN

# =============================================================================
"""
def main():

    # --- Parser ---
    parser = argparse.ArgumentParser(
        description='Filters signature candidates which match sufficiently to \
        any target in an exclusion database. Any signatures which to do not \
        align sufficiently to any exclusion target are written to file and \
        are considered filtered signatures.')

    parser.add_argument(
        INCLUSION_DATABASE_SHORT,
        INCLUSION_DATABASE_LONG,
        dest=INCLUSION_DATABASE,
        help="inclusion database location",
        type=str, required=True)

    parser.add_argument(
        EXCLUSION_DATABASE_SHORT,
        EXCLUSION_DATABASE_LONG,
        dest=EXCLUSION_DATABASE,
        help="exclusion database location",
        type=str, required=True)

    parser.add_argument(
        INCLUSION_SHORT,
        INCLUSION_LONG,
        dest=INCLUSION,
        help="inclusion genome(s)",
        type=str, required=True, nargs='+')

    parser.add_argument(
        EXCLUSION_SHORT,
        EXCLUSION_LONG,
        dest=EXCLUSION,
        help="exclusion genome(s)",
        type=str, required=True, nargs='+')

    parser.add_argument(
        INPUT_SHORT,
        INPUT_LONG,
        dest=INPUT,
        help="candidates location",
        type=str, required=True)

    parser.add_argument(
        FILTERED_OUTPUT_SHORT,
        FILTERED_OUTPUT_LONG,
        dest=FILTERED_OUTPUT,
        help="filtered output location",
        type=str, required=True)

    parser.add_argument(
        SORTED_OUTPUT_SHORT,
        SORTED_OUTPUT_LONG,
        dest=SORTED_OUTPUT,
        help="sorted output location",
        type=str, required=True)

    parser.add_argument(
        FILTER_PERCENT_SHORT,
        FILTER_PERCENT_LONG,
        dest=FILTER_PERCENT,
        help="the maximum percent identity of an exclusion hit",
        type=float, required=False, default=0.50)

    parser.add_argument(
        FILTER_LENGTH_SHORT,
        FILTER_LENGTH_LONG,
        dest=FILTER_LENGTH,
        help="the maximum shared fractional length of an exclusion hit \
            with a candidate",
        type=float, required=False, default=0.50)

    args = parser.parse_args()

    inclusionDatabaseLocation = args.inclusionDatabase
    exclusionDatabaseLocation = args.exclusionDatabase
    totalInclusion = len(args.inclusion)
    totalExclusion = len(args.exclusion)
    inputLocation = args.input
    filteredOutputLocation = args.filteredOutput
    sortedOutputLocation = args.sortedOutput
    filterLength = args.filterLength
    filterPercent = args.filterPercent

    filterSignatures(
        inclusionDatabaseLocation, exclusionDatabaseLocation,
        totalInclusion, totalExclusion, inputLocation,
        filteredOutputLocation, sortedOutputLocation, filterLength,
        filterPercent)

if __name__ == '__main__':

    main()
