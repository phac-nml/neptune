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

MAIN

# =============================================================================
"""
def main():

    # --- Parser ---
    parser = argparse.ArgumentParser(
        description='.')

    parser.add_argument(
        "-s",
        "--signatures",
        dest="signatures",
        help="signatures",
        type=str, required=True, nargs='+')

    args = parser.parse_args()

    signatureFileLocations = []
    Utility.expandInput(args.signatures, signatureFileLocations)

    masterSignatures = {}

    fileNumber = 0

    for fileLocation in signatureFileLocations:

        signatures = Signature.readSignatures(fileLocation)

        for signatureID in signatures:
            masterID = str(fileNumber) + "." + signatureID
            masterSignatures[masterID] = signatures[signatureID]
            masterSignatures[masterID].ID = masterID

        fileNumber += 1

    sortedSignatures = sorted(masterSignatures.iteritems(), key=lambda (k,v): v.score, reverse=True)

    outputFile = open("consolodatedSignatures.fasta", 'w')

    for item in sortedSignatures:

        ID = item[0]
        Signature.writeSignature(masterSignatures[ID], outputFile)

    outputFile.close()

    # make db !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!?

    # COMMAND LINE
    COMMAND = "makeblastdb"

    TYPE = "-dbtype"
    NUCLEOTIDE = "nucl"

    INPUT = "-in"
    INPUT_LOCATIONS = "consolodatedSignatures.fasta"

    TITLE = "-title"
    NAME = "DATABASE"

    OUTPUT = "-out"
    OUTPUT_LOCATION = "consolodatedSignatures.db"

    args = []

    args.append(COMMAND)

    args.append(TYPE)
    args.append(NUCLEOTIDE)

    args.append(INPUT)
    args.append(INPUT_LOCATIONS)

    args.append(TITLE)
    args.append(NAME)

    args.append(OUTPUT)
    args.append(OUTPUT_LOCATION)

    print args

    subprocess.check_call(args)

    Database.queryDatabase("consolodatedSignatures.db", "consolodatedSignatures.fasta", "db.out", 0.50)

    blastOutput = open("db.out.query")
    hits = {}
    outputSignatures = {}

    for line in blastOutput:

        hit = Database.Hit(line)

        # only care if longer
        if (float(hit.alignmentLength) / float(hit.length) < float(0.50)):
            continue

        if hit.ID in hits:
            hits[hit.ID].append(hit.reference)

        else:
            hits[hit.ID] = [hit.reference]

    for signatureID in masterSignatures:

        if(all((ID not in outputSignatures) for ID in hits[signatureID])):

            outputSignatures[signatureID] = masterSignatures[signatureID]

    sortedOutputSignatures = sorted(outputSignatures.iteritems(), key=lambda (k,v): v.score, reverse=True)

    outputFile = open("finalSignatures.fasta", 'w')

    for item in sortedOutputSignatures:

        ID = item[0]
        Signature.writeSignature(masterSignatures[ID], outputFile)

    outputFile.close()

    print "\n==== Exiting ====\n"

if __name__ == '__main__':

    main()
