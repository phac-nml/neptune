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
import operator

import Database
import Signature

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
SEED_SIZE = "seedSize"

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
SEED_SIZE_LONG = LONG + SEED_SIZE

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
SEED_SIZE_SHORT = SHORT + "ss"

# DEFAULTS

FILTER_PERCENT_DEFAULT = 0.50
FILTER_LENGTH_DEFAULT = 0.50
SEED_SIZE_DEFAULT = 11

"""
# =============================================================================

FILTER SIGNATURES

# =============================================================================
"""
class FilterSignatures():

    """
    # =========================================================================

    INIT

    PURPOSE:
        Constructs the FilterSignatures object.

    INPUT:
        [FILE LOCATION] [candidatesLocation] - The location of the candidate
            signatures.
        [FILE LOCATION] [filteredLocation] - The file location to write the
            filtered signatures.
        [FILE LOCATION] [sortedLocation] - The file location to write the
            sorted signatures.
        [1 <= INT] [totalInclusion] - The total number of inclusion targets.
        [1 <= INT] [totalExclusion] - The total number of exclusion targets.
        [0 <= FLOAT 0 <= 1] [filterLength] - The minimum query alignment length
            for the signature to be considered a hit and used in scoring.

    # =========================================================================
    """
    def __init__(
            self, candidatesLocation, filteredLocation, sortedLocation,
            totalInclusion, totalExclusion, filterLength):

        self.candidatesLocation = candidatesLocation    # Candidate signatures.
        self.filteredLocation = filteredLocation        # Filtered signatures.
        self.sortedLocation = sortedLocation            # Sorted signatures.

        self.totalInclusion = totalInclusion        # Number of inclusion.
        self.totalExclusion = totalExclusion        # Number of exclusion.

        self.filterLength = filterLength            # The min filtering size.

        self.hitOverallDictionary = {}  # The best (hit.ID, X).
        self.hitPairDictionary = {}     # The best (hit.ID, hit.reference).

        self.overallScore = {}          # The overall score of a signature.
        self.inclusionScore = {}        # The inclusion score of a signature.
        self.exclusionScore = {}        # The exclusion score of a signature.

    """
    # =========================================================================

    UPDATE HIT OVERALL DICTIONARY

    PURPOSE:
        The overall score dictionary will be updated with the passed hit. This
        function is looking to keep the best hit for a query, regardless of the
        subject (reference).

    INPUT:
        [HIT] [hit] - The database hit object associated with the hit.

    RETURN:
        [NONE]

    POST:
        The [self.overallScore] dictionary will be updated.

    # =========================================================================
    """
    def updateHitOverallDictionary(self, hit):

        key = hit.ID                                # The key: ID of the hit.
        dictionary = self.hitOverallDictionary      # The dictionary to update.

        if key not in dictionary:

            dictionary[key] = hit

        elif hit.alignmentScore > dictionary[key].alignmentScore:

            dictionary[key] = hit

    """
    # =========================================================================

    UPDATE HIT PAIR DICTIONARY

    PURPOSE:
        The pair score dictionary will be updated with the passed hit. This
        function is looking for the best hit for a (hit, reference) pair.

    INPUT:
        [HIT] [hit] - The hit object associated with the hit.

    RETURN:
        [NONE]

    POST:
        The [self.hitPairDictionary] object will be updated.

    # =========================================================================
    """
    def updateHitPairDictionary(self, hit):

        key = (hit.ID, hit.reference)           # The key: ID and reference.
        dictionary = self.hitPairDictionary     # The dictionary to update.

        if key not in dictionary:

            hit.neptuneScore = (
                (float(hit.alignmentLength) / float(hit.length)) *
                (float(hit.percentIdentity) / float(100)))

            dictionary[key] = hit

        elif hit.alignmentScore > dictionary[key].alignmentScore:

            hit.neptuneScore = (
                (float(hit.alignmentLength) / float(hit.length)) *
                (float(hit.percentIdentity) / float(100)))

            dictionary[key] = hit

    """
    # =========================================================================

    UPDATE EXCLUSION SCORES

    PURPOSE:
        Updates the exclusion score and negative component of the overall
        score.

    RETURN:
        [NONE]

    POST:
        The exclusion scores and negative component of the overall scores of
        [self.overallScore] and [self.exclusionScore] will be updated using the
        [self.hitPairDictionary].

    # =========================================================================
    """
    def updateExclusionScores(self):

        dictionary = self.hitPairDictionary     # The pair dictionary.

        for key in dictionary:

            hit = dictionary[key]
            ID = hit.ID
            neptuneScore = hit.neptuneScore

            if ID not in self.overallScore:
                self.overallScore[ID] = (
                    -float(neptuneScore) / float(self.totalExclusion))

            else:
                self.overallScore[ID] -= (
                    float(neptuneScore) / float(self.totalExclusion))

            if ID not in self.exclusionScore:
                self.exclusionScore[ID] = (
                    -float(neptuneScore) / float(self.totalExclusion))

            else:
                self.exclusionScore[ID] -= (
                    float(neptuneScore) / float(self.totalExclusion))

    """
    # =========================================================================

    UPDATE INCLUSION SCORES

    PURPOSE:
        Updates the inclusion score and positive component of the overall
        score.

    RETURN:
        [NONE]

    POST:
        The exclusion scores and negative component of the overall scores of
        [self.overallScore] and [self.inclusionScore] will be updated using the
        [self.hitPairDictionary].

    # =========================================================================
    """
    def updateInclusionScores(self):

        dictionary = self.hitPairDictionary     # The pair dictionary.

        for key in dictionary:

            hit = dictionary[key]
            ID = hit.ID
            neptuneScore = hit.neptuneScore

            if ID not in self.overallScore:
                self.overallScore[ID] = (
                    float(neptuneScore) / float(self.totalInclusion))

            else:
                self.overallScore[ID] += (
                    float(neptuneScore) / float(self.totalInclusion))

            if ID not in self.inclusionScore:
                self.inclusionScore[ID] = (
                    float(neptuneScore) / float(self.totalInclusion))

            else:
                self.inclusionScore[ID] += (
                    float(neptuneScore) / float(self.totalInclusion))

    """
    # =========================================================================

    REPORT FILTERED CANDIDATES

    PURPOSE:
        Reports all of the filtered candidate signatures to output. These are
        the candidate signatures that were are not immediately removed by
        filtering against the exclusion database.

    RETURN:
        [NONE]

    POST:
        The filtered candidates will be written to the [self.filtedLocation].

    # =========================================================================
    """
    def reportFilteredCandidates(self):

        outputFile = open(self.filteredLocation, 'w')   # The output.
        candidateSignatures = Signature.readSignatures(
            self.candidatesLocation)                    # The input.
        dictionary = self.hitOverallDictionary          # Overall dictionary.

        for ID in candidateSignatures:

            signature = candidateSignatures[ID]

            if ID in dictionary:

                hit = dictionary[ID]

                if (float(hit.alignmentLength) / float(signature.length) <
                        float(self.filterLength)):
                    Signature.writeSignature(signature, outputFile)

            else:

                Signature.writeSignature(signature, outputFile)

        outputFile.close()

    """
    # =========================================================================

    REPORT SORTED

    PURPOSE:
        Reports all of the filtered signatures in sorted order to output. This
        function does not sort the signatures. The sorted order is informed by
        the [sortedSignatureIDs] parameter.

    INPUT:
        [(SIGNATURE ID) LIST] [sortedSignatureIDs] - An iterable list, in
            sorted order, of signature IDs.

    RETURN:
        [NONE]

    POST:
        The sorted signatures will be written to the [self.filteredLocation].

    # =========================================================================
    """
    def reportSorted(self, sortedSignatureIDs):

        filteredSignatures = (
            Signature.readSignatures(self.filteredLocation))    # The output.

        # REPORT SORTED SIGNATURES
        outputFile = open(self.sortedLocation, 'w')

        for ID in sortedSignatureIDs:

            if ID in filteredSignatures:

                signature = filteredSignatures[ID]

                # -- Score Signature -- #
                if ID in self.overallScore:
                    signature.score = self.overallScore[ID]
                else:
                    signature.score = 0.0

                if ID in self.inclusionScore:
                    signature.inscore = self.inclusionScore[ID]
                else:
                    signature.inscore = 0.0

                if ID in self.exclusionScore:
                    signature.exscore = self.exclusionScore[ID]
                else:
                    signature.exscore = 0.0

                Signature.writeSignature(signature, outputFile)

        outputFile.close()

    """
    # =========================================================================

    REPORT SIGNATURES

    PURPOSE:
        Reports the candidate signatures which are not found in the list of
        filterable signatures.

    INPUT:
        [FILE LOCATION] [exclusionQueryLocation] - The file location of the
            signatures to filter.

    RETURN:
        [NONE]

    POST:
        A file of filterted signatures will be produced at
        [self.filteredLocation].

    # =========================================================================
    """
    def reportSignatures(self, exclusionQueryLocation):

        exclusionQueryFile = open(exclusionQueryLocation, 'r')

        # LOAD FILTER
        for line in exclusionQueryFile:

            hit = Database.Hit(line)
            self.updateHitOverallDictionary(hit)
            self.updateHitPairDictionary(hit)

        exclusionQueryFile.close()

        self.updateExclusionScores()

        self.reportFilteredCandidates()

    """
    # =========================================================================

    SORT

    PURPOSE:
        Sorts the filtered signatures according to their signature score.

    INPUT:
        [FILE LOCATION] [inclusionQueryLocation] - The file location of the
            filtered signatures against the inclusion targets.

    RETURN:
        [NONE]

    POST:
        The [self.filteredLocation] signatures will be written to the
        [self.sortedLocation] in score-descending order.

    # =========================================================================
    """
    def sortSignatures(self, inclusionQueryLocation):

        inclusionQueryFile = open(inclusionQueryLocation, 'r')

        # LOAD FILTER
        for line in inclusionQueryFile:

            hit = Database.Hit(line)
            self.updateHitPairDictionary(hit)

        inclusionQueryFile.close()

        self.updateInclusionScores()

        sortedSignatureIDs = [ID for (ID, score) in sorted(
            self.overallScore.items(), key=operator.itemgetter(1),
            reverse=True)]

        self.reportSorted(sortedSignatureIDs)

"""
# =============================================================================

FILTER SIGNATURES

PURPOSE:
    Reports the filtered signatures, both as a list of candidates and in sorted
    order.

INPUT:
    [FILE LOCATION] [inclusionDatabaseLocation] - The file location of the
        inclusion database.
    [FILE LOCATION] [exclusionDatabaseLocation] - The file location of the
        exclusion database.
    [1 <= INT] [totalInclusion] - The total number of inclusion targets.
    [1 <= INT] [totalExclusion] - The total number of exclusion targets.
    [FILE LOCATION] [candidatesLocation] - The location of the candidate
        signatures.
    [FILE LOCATION] [filteredOutputLocation] - The file location to write the
        filtered signatures.
    [FILE LOCATION] [sortedOutputLocation] - The file location to write the
        sorted signatures.
    [0 <= FLOAT 0 <= 1] [filterLength] - The minimum query alignment length
        for the signature to be considered a hit and used in scoring.
    [0 <= FLOAT 0 <= 1] [filterPercent] - The maximum percent identity of an
        exclusion hit with a candidate.
    [4 <= INT] [seedSize] - The seed size used in alignments.

RETURN:
    [NONE]

POST:
    Filtered signatures will be written to [filteredOutputLocation] and
    sorted signatures will be written to [sortedOutputLocation].

# =============================================================================
"""
def filterSignatures(
        inclusionDatabaseLocation, exclusionDatabaseLocation,
        totalInclusion, totalExclusion, candidatesLocation,
        filteredOutputLocation, sortedOutputLocation, filterLength,
        filterPercent, seedSize):

    filterSignatures = FilterSignatures(
        candidatesLocation, filteredOutputLocation, sortedOutputLocation,
        totalInclusion, totalExclusion, filterLength)

    # QUERY DB - EXCLUSION
    exclusionQueryLocation = Database.queryDatabase(
        exclusionDatabaseLocation, candidatesLocation,
        filteredOutputLocation, filterPercent, seedSize)

    # FILTER
    filterSignatures.reportSignatures(exclusionQueryLocation)

    # QUERY DB - INCLUSION
    inclusionQueryLocation = Database.queryDatabase(
        inclusionDatabaseLocation, filteredOutputLocation,
        sortedOutputLocation, filterPercent, seedSize)

    # SORT
    filterSignatures.sortSignatures(inclusionQueryLocation)

"""
# =============================================================================

PARSE

# =============================================================================
"""
def parse(parameters):

    inclusionDatabaseLocation = parameters[INCLUSION_DATABASE]
    exclusionDatabaseLocation = parameters[EXCLUSION_DATABASE]
    totalInclusion = len(parameters[INCLUSION])
    totalExclusion = len(parameters[EXCLUSION])
    inputLocation = parameters[INPUT]
    filteredOutputLocation = parameters[FILTERED_OUTPUT]
    sortedOutputLocation = parameters[SORTED_OUTPUT]

    filterLength = parameters[FILTER_LENGTH] \
        if parameters[FILTER_LENGTH] else FILTER_LENGTH_DEFAULT

    filterPercent = parameters[FILTER_PERCENT] \
        if parameters[FILTER_PERCENT] else FILTER_PERCENT_DEFAULT

    seedSize = parameters[SEED_SIZE] \
        if parameters[SEED_SIZE] else SEED_SIZE_DEFAULT

    filterSignatures(
        inclusionDatabaseLocation, exclusionDatabaseLocation,
        totalInclusion, totalExclusion, inputLocation,
        filteredOutputLocation, sortedOutputLocation, filterLength,
        filterPercent, seedSize)

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
        type=float, required=False)

    parser.add_argument(
        FILTER_LENGTH_SHORT,
        FILTER_LENGTH_LONG,
        dest=FILTER_LENGTH,
        help="the maximum shared fractional length of an exclusion hit \
            with a candidate",
        type=float, required=False)

    parser.add_argument(
        SEED_SIZE_SHORT,
        SEED_SIZE_LONG,
        dest=SEED_SIZE,
        help="the seed size used during alignment",
        type=int, required=False)

    args = parser.parse_args()
    parameters = vars(args)
    parse(parameters)

"""
# =============================================================================
# =============================================================================
"""
if __name__ == '__main__':

    main()
