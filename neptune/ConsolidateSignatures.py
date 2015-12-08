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
Date: 7 December 2015

This script consolidates multiple signatures files, produced by Neptune, into
a single file containing the best signatures from all files. This script
attempts to avoid overlapping signatures. However, this is not guaranteed.

# =============================================================================
"""

import argparse
import os

import Signature
import Database
import Utility

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
    Compiles the signatures from several Neptune signature files into a single
    dictionary containing all signatures.

INPUT:
    [[STRING ID] -> [SIGNATURE] DICTIONARY] [compiledSignatures] - An initially
        empty dictionary that will be filled with signatures located in within
        the [signatureLocations] files.
    [FILE LOCATION LIST] [signatureLocations] - A list of signature file
        locations.

RETURN:
    [[STRING ID] -> [SIGNATURE] DICTIONARY] [compiledSignatures] - A dictionary
        containing all compiled signatures. This dictionary is the same object
        as the initially passed [compiledSignatures] dictionary.

POST:
    The [compiledSignatures] dictionary will be filled the signatures.

# =============================================================================
"""
def compileSignatures(compiledSignatures, signatureLocations):

    fileID = 0

    # -- Read Files -- #
    for location in signatureLocations:

        signatures = Signature.readSignatures(location)

        for signatureID in signatures:

            compileID = str(fileID) + "." + signatureID
            compiledSignatures[compileID] = signatures[signatureID]
            compiledSignatures[compileID].ID = compileID

        fileID += 1

    return compiledSignatures

"""
# =============================================================================

PRODUCE SIGNATURES

PURPOSE:
    Produces a list of consolidated signatures by outputting to file while
    avoiding outputting duplicate signatures.

INPUT:
    [SIGNATURE LIST] [sortedSignatures] - A list of signatures, sorted by their
        corresponding Neptune signature scores.
    [FILE] [queryFile] - A readable BLASTN query file.
    [FILE] [destination] - A writable file-like object.

RETURN:
    [NONE]

POST:
    The list of consolidated signatures will be written to the [destination].

# =============================================================================
"""
def produceSignatures(sortedSignatures, queryFile, destination):

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
    for signature in sortedSignatures:

        # Is the signature close to anything already output?
        if(all((ID not in outputSignatures) for ID in hits[signature.ID])):

            outputSignatures[signature.ID] = signature
            Signature.writeSignature(signature, destination)

"""
# =============================================================================

CONSOLIDATE SIGNATURES

PURPOSE:
    Consolidates signatures from several Neptune signature files into a single
    representative Neptune signature file, determined by signature score and
    sequence similarity.

INPUT:
    [FILE LOCATION LIST] [signatureLocations] - A list of Neptune signature
        file locations corresponding to files to consolidate.
    [FILE DIRECTORY LOCATION] [outputDirectoryLocation] - The directory to
        write the output files.

RETURN:
    [NONE]

POST:
    The signatures files will be written to several locations within the
    [outputDirectoryLocation].

# =============================================================================
"""
def consolidateSignatures(signatureLocations, outputDirectoryLocation):

    # --- Compile Signatures --- #
    compiledSignatures = {}
    compileSignatures(compiledSignatures, signatureLocations)

    # -- Sort Signatures -- #
    sortedSignatures = Signature.sortSignatures(compiledSignatures)

    # -- Write Signatures -- #
    compiledSignatureLocation = os.path.join(
        outputDirectoryLocation, COMPILED_SIGNATURES)

    compiledSignatureFile = open(compiledSignatureLocation, 'w')
    Signature.writeSignatures(sortedSignatures, compiledSignatureFile)
    compiledSignatureFile.close()

    # --- Build and Query Database --- #
    databaseLocation = os.path.join(
        outputDirectoryLocation, COMPILED_DATABASE)
    queryLocation = os.path.join(
        outputDirectoryLocation, COMPILED_DATABASE_QUERY)

    Database.createDatabaseJob(compiledSignatureLocation, databaseLocation)
    Database.queryDatabase(
        databaseLocation, compiledSignatureLocation, queryLocation, 0.50, 11)

    # --- Produce Signatures --- #
    outputLocation = os.path.join(
        outputDirectoryLocation, CONSOLIDATED_SIGNATURES)
    outputFile = open(outputLocation, 'w')
    queryFile = open(queryLocation, 'r')

    produceSignatures(sortedSignatures, queryFile, outputFile)

    outputFile.close()
    queryFile.close()

    print "\n==== Exiting ====\n"

"""
# =============================================================================

MAIN

# =============================================================================
"""
def main():

    # --- Parser ---
    parser = argparse.ArgumentParser(
        description='Consolidates signatures from several Neptune signature \
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
