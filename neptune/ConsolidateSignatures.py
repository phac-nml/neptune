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

CONSOLIDATE SIGNATURES

# =============================================================================
"""

import math
import argparse
import os

import sys
import operator
import subprocess

import Signature
import Database
import Utility

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

"""
# =============================================================================

GLOBALS

# =============================================================================
"""

# NAMES
SIGNATURES = "signatures"
OUTPUT = "output"
CONSOLIDATED_DATABASE = "database"

# ARGUMENTS
LONG = "--"

SIGNATURES_LONG = LONG + SIGNATURES
CONSOLIDATED_DATABASE_LONG = LONG + CONSOLIDATED_DATABASE
OUTPUT_LONG = LONG + OUTPUT

SHORT = "-"

SIGNATURES_SHORT = SHORT + "s"
CONSOLIDATED_DATABASE_SHORT = SHORT + CONSOLIDATED_DATABASE
OUTPUT_SHORT = SHORT + "o"

# OTHER
COMPILED_SIGNATURES = "compiled.fasta"
COMPILED_DATABASE = "compiled.db"
COMPILED_DATABASE_QUERY = COMPILED_DATABASE + ".query"
CONSOLIDATED_SIGNATURES = "consolidated.fasta"

"""
# =============================================================================

COMPILE SIGNATURES

PURPOSE:

INPUT:
    [[STRING ID] -> [SIGNATURE] DICTIONARY] compiledSignatures
    [signatureLocations]
    [outputLocation]

RETURN:
    [sortedIDs]

POST:

# =============================================================================
"""
def compileSignatures(compiledSignatures, signatureLocations, outputLocation):    

    fileID = 0

    for location in signatureLocations:

        signatures = Signature.readSignatures(location)

        for signatureID in signatures:

            compileID = str(fileID) + "." + signatureID
            compiledSignatures[compileID] = signatures[signatureID]
            compiledSignatures[compileID].ID = compileID

        fileID += 1

    sortedSignatures = sorted(
        compiledSignatures.iteritems(),
        key=lambda (k,v): v.score, reverse=True)

    sortedIDs = [item[0] for item in sortedSignatures]

    outputFile = open(outputLocation, 'w')

    for item in sortedSignatures:

        ID = item[0]
        Signature.writeSignature(compiledSignatures[ID], outputFile)

    outputFile.close()

    return sortedIDs

"""
# =============================================================================

PRODUCE SIGNATURES

PURPOSE:

INPUT:
    [[STRING ID] -> [SIGNATURE] DICTIONARY] compiledSignatures
    [STRING LIST] [sortedIDs]
    [FILE LOCATION] [queryLocation]
    [outputFile]

RETURN:
    [NONE]

POST:
    The consolidated signatures will be written to the [outputFile].

# =============================================================================
"""
def produceSignatures(
        compiledSignatures, sortedIDs, queryLocation, outputFile):

    queryFile = open(queryLocation)

    hits = {}
    outputSignatures = {}

    # Build a list of all query hits.
    for line in queryFile:

        hit = Database.Hit(line)

        # only care if longer
        if (float(hit.alignmentLength) / float(hit.length) < float(0.50)):
            continue

        if hit.ID in hits:
            hits[hit.ID].append(hit.reference)

        else:
            hits[hit.ID] = [hit.reference]

    # Build a list of output signatures.
    for signatureID in sortedIDs:

        # Is the signature close to anything already output?
        if(all((ID not in outputSignatures) for ID in hits[signatureID])):

            outputSignatures[signatureID] = compiledSignatures[signatureID]
            Signature.writeSignature(compiledSignatures[signatureID], outputFile)

"""
# =============================================================================

CONSOLIDATE SIGNATURES

PURPOSE:
    Consolidates signatures from several FASTA signature files into a single
    representative signature file, determined by signature score and sequence
    similarity.

INPUT:
    [FILE LOCATION LIST] [signatureLocations] - A list of signature file
        locations to consolidate.
    [FILE DIRECTORY LOCATION] [outputDirectoryLocation] - The directory to
        write output files.        

RETURN:
    [NONE]

POST:
    The signatures will be consolidated and written to the [outputFile].

# =============================================================================
"""
def consolidateSignatures(signatureLocations, outputDirectoryLocation):

    # --- Compile ---
    compiledSignatures = {}
    compiledSignatureLocation = COMPILED_SIGNATURES
    sortedIDs = compileSignatures(
        compiledSignatures, signatureLocations, compiledSignatureLocation)

    # --- Build Database ---
    databaseLocation = os.path.join(
        outputDirectoryLocation, COMPILED_DATABASE)
    queryLocation = os.path.join(
        outputDirectoryLocation, COMPILED_DATABASE_QUERY)

    Database.createDatabaseJob(compiledSignatureLocation, databaseLocation)
    Database.queryDatabase(
        databaseLocation, compiledSignatureLocation, queryLocation, 0.50)

    # --- Produce Consolidated Signatures ---
    outputLocation = os.path.join(
        outputDirectoryLocation, CONSOLIDATED_SIGNATURES)
    outputFile = open(outputLocation, 'w')

    produceSignatures(
        compiledSignatures, sortedIDs, queryLocation, outputFile)

    outputFile.close()

    print "\n==== Exiting ====\n"    

"""
# =============================================================================

MAIN

# =============================================================================
"""
def main():

    # --- Parser ---
    parser = argparse.ArgumentParser(
        description='Consolidates signatures from several FASTA signature \
        files into a single representative signature file, determined by \
        signature score and sequence similarity.')

    parser.add_argument(
        SIGNATURES_SHORT,
        SIGNATURES_LONG,
        dest=SIGNATURES,
        help="file locations of all signatures to consolidate",
        type=str, required=True, nargs='+')

    parser.add_argument(
        OUTPUT_SHORT,
        OUTPUT_LONG,
        dest=OUTPUT,
        help="output directory",
        type=str, required=True)

    args = parser.parse_args()

    signatureLocations = []
    Utility.expandInput(args.signatures, signatureLocations)
    outputDirectoryLocation = args.output

    consolidateSignatures(signatureLocations, outputDirectoryLocation)

if __name__ == '__main__':

    main()
