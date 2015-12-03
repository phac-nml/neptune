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
CONSOLODATED_DATABASE = "database"

# ARGUMENTS
LONG = "--"

SIGNATURES_LONG = LONG + SIGNATURES
CONSOLODATED_DATABASE_LONG = LONG + CONSOLODATED_DATABASE
OUTPUT_LONG = LONG + OUTPUT

SHORT = "-"

SIGNATURES_SHORT = SHORT + "s"
CONSOLODATED_DATABASE_SHORT = SHORT + CONSOLODATED_DATABASE
OUTPUT_SHORT = SHORT + "o"

"""
# =============================================================================

COMPILE SIGNATURES

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

    sortedSignatures = sorted(compiledSignatures.iteritems(), key=lambda (k,v): v.score, reverse=True)

    outputFile = open(outputLocation, 'w')

    for item in sortedSignatures:

        ID = item[0]
        Signature.writeSignature(compiledSignatures[ID], outputFile)

    outputFile.close()

"""
# =============================================================================

PRODUCE CONSOLODATED SIGNATURES

# =============================================================================
"""
def produceConsolodatedSignatures(compiledSignatures, queryLocation, outputLocation):

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

    outputFile = open(outputLocation, 'w')

    # Build a list of output signatures.
    for signatureID in compiledSignatures:

        # Is the signature close to anything already output?
        if(all((ID not in outputSignatures) for ID in hits[signatureID])):

            outputSignatures[signatureID] = compiledSignatures[signatureID]
            Signature.writeSignature(compiledSignatures[signatureID], outputFile)

    outputFile.close()

"""
# =============================================================================

MAIN

# =============================================================================
"""
def main():

    # --- Parser ---
    parser = argparse.ArgumentParser(
        description='.')

    parser.add_argument(
        SIGNATURES_SHORT,
        SIGNATURES_LONG,
        dest=SIGNATURES,
        help="file locations of all signatures to consolodate",
        type=str, required=True, nargs='+')

    parser.add_argument(
        CONSOLODATED_DATABASE_SHORT,
        CONSOLODATED_DATABASE_LONG,
        dest=CONSOLODATED_DATABASE,
        help="the location to create the consolodated database location",
        type=str, required=True)

    parser.add_argument(
        OUTPUT_SHORT,
        OUTPUT_LONG,
        dest=OUTPUT,
        help="output file",
        type=str, required=True)

    args = parser.parse_args()

    signatureLocations = []
    Utility.expandInput(args.signatures, signatureLocations)

    databaseLocation = args.database
    outputLocation = args.output

    compiledSignatures = {}
    compiledSignatureLocation = "consolodatedSignatures.fasta"
    compileSignatures(compiledSignatures, signatureLocations, compiledSignatureLocation)

    Database.createDatabaseJob(compiledSignatureLocation, databaseLocation)
    Database.queryDatabase(databaseLocation, compiledSignatureLocation, "db.out", 0.50)

    produceConsolodatedSignatures(compiledSignatures, "db.out.query", outputLocation)

    print "\n==== Exiting ====\n"

if __name__ == '__main__':

    main()
