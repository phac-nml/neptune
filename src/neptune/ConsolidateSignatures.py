#!/usr/bin/env python

"""
# =============================================================================

Copyright Government of Canada 2015-2017

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

This script consolidates multiple signatures files, produced by Neptune, into
a single file containing the best signatures from all files. This script
attempts to avoid overlapping signatures. However, this is not guaranteed.

# =============================================================================
"""

import argparse
import os

import neptune.Signature as Signature
import neptune.Database as Database
import neptune.Utility as Utility

"""
# =============================================================================

GLOBALS

# =============================================================================
"""

PROGRAM_DESCRIPTION = 'Consolidates signatures from several Neptune signature \
    files (FASTA format) into a single representative signature file, \
    determined by signature score and sequence similarity.'

# DEFAULTS #

SEED_SIZE_DEFAULT = 11

# ARGUMENTS #

LONG = "--"
SHORT = "-"

# REQUIRED ARGUMENTS #

# Signatures
SIGNATURES = "signatures"
SIGNATURES_LONG = LONG + SIGNATURES
SIGNATURES_SHORT = SHORT + "s"
SIGNATURES_HELP = "The file locations of all signatures to consolidate."

# Output
OUTPUT = "output"
OUTPUT_LONG = LONG + OUTPUT
OUTPUT_SHORT = SHORT + "o"
OUTPUT_HELP = "The output directory to place the consolidate signatures and \
    any additional files."

# OPTIONAL ARGUMENTS #

# Seed Size
SEED_SIZE = "seed-size"
SEED_SIZE_LONG = LONG + SEED_SIZE
SEED_SIZE_SHORT = SHORT + "ss"
SEED_SIZE_HELP = "The seed size used during sequence alignment."

# OTHER #

COMPILED_SIGNATURES = "compiled.fasta"
COMPILED_DATABASE = "compiled.db"
COMPILED_DATABASE_QUERY = COMPILED_DATABASE + ".query"
CONSOLIDATED_SIGNATURES = "consolidated.fasta"

"""
# =============================================================================

COMPILE SIGNATURES
------------------


PURPOSE
-------

Compiles the signatures from several Neptune signature files (FASTA format)
into a single dictionary containing all signatures. This may result in repeated
signatures.


INPUT
-----

[(STRING ID) -> (SIGNATURE) DICTIONARY] [compiledSignatures]
    An initially empty dictionary that will be filled with signatures located
    in within the [signatureLocations] files.

[(FILE LOCATION) LIST] [signatureLocations]
    A list of signature file locations from which to compile signatures from.


RETURN
------

[(STRING ID) -> (SIGNATURE) DICTIONARY] [compiledSignatures]
    A dictionary containing all compiled signatures. This dictionary is the
    same object as the initially passed [compiledSignatures] dictionary.


POST
----

The [compiledSignatures] dictionary will be filled with the signatures.

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
------------------


PURPOSE
-------

Produces a list of consolidated signatures by outputting signatures to a file,
while attempting to avoid outputting duplicate signatures.


INPUT
-----

[SIGNATURE LIST] [sortedSignatures]
    A list of signatures, sorted by their corresponding Neptune signature
    scores. This list of signatures may contain apparently-duplicate
    signatures.

[FILE] [blastOutputFile]
    A readable BLASTN output file. This query is the output of aligning all
    [sortedSignatures] against themselves.

[FILE] [destination]
    A writable file-like object to write the consolidated signatures.


RETURN
------

[NONE]


POST
----

The list of consolidated signatures will be written to the [destination].

# =============================================================================
"""
def produceSignatures(sortedSignatures, blastOutputFile, destination):

    hits = {}  # [SIGNATURE ID] -> [(SIGNATURE ID) LIST] // (alignments)
    outputSignatures = {}  # Collection of already-output signatures.

    # Pre-populate hits map with empty lists:
    for signature in sortedSignatures:
        hits[signature.ID] = []

    # Build a list of all query hits.
    # This creates a dictionary mapping signatures that align to each other.
    # [SIGNATURE ID] -> [(SIGNATURE ID) LIST]
    for line in blastOutputFile:

        hit = Database.Hit(line)

        # We only keep the hit if the ratio of the signature-to-alignment
        # length is sufficiently long.
        if (float(hit.alignmentLength) / float(hit.length) >= float(0.50)):
            # Append the signature ID to the existing list of IDs.
            hits[hit.ID].append(hit.reference)

    # Write the signatures to output, while maintaining a dictionary of
    # signatures that were previously written to output. This attempts to
    # avoid writing signatures appear to be duplicates or appear to overlap
    # significantly.
    for signature in sortedSignatures:

        # Is the signature close to anything already written to output?
        if (all((ID not in outputSignatures) for ID in hits[signature.ID])):

            # The signature appears to be sufficiently unique.
            # Write the signature to output and update outputed signatures.
            outputSignatures[signature.ID] = signature
            Signature.writeSignature(signature, destination)


"""
# =============================================================================

CONSOLIDATE SIGNATURES
----------------------


PURPOSE
-------

Consolidates signatures from several Neptune signature files into a single
representative Neptune signature file, determined by the signature score and
sequence similarity of all the contributing signatures.

The function compiles all the signatures into a single dictionary file, sorts
these signatures according to their Neptune signature score, and writes these
sorted signatures to a file. It then uses BLAST to query the signatures against
themselves and uses this information to report signatures in a greedy manner.
The signatures are reported in an order according to their signature score and
only if there has been no other similar signature (determined by BLAST) that
has already been reported.


INPUT
-----

[(FILE LOCATION) LIST] [signatureLocations]
    A list of Neptune signature file locations corresponding to files to
    consolidate.

[4 <= INT] [seedSize]
    The seed size used in alignments to determine similarity.

[(FILE DIRECTORY) LOCATION] [outputDirectoryLocation]
    The directory to write the output files.


RETURN
------

[NONE]

POST
----

The signatures and associated files will be written to several locations within
the [outputDirectoryLocation].

# =============================================================================
"""
def consolidateSignatures(
        signatureLocations, seedSize, outputDirectoryLocation):

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
        databaseLocation, compiledSignatureLocation,
        queryLocation, 0.50, seedSize)

    # --- Produce Signatures --- #
    outputLocation = os.path.join(
        outputDirectoryLocation, CONSOLIDATED_SIGNATURES)
    outputFile = open(outputLocation, 'w')
    queryFile = open(queryLocation, 'r')

    produceSignatures(sortedSignatures, queryFile, outputFile)

    outputFile.close()
    queryFile.close()

    # --- Clean Output --- #
    filelist = [f for f in os.listdir(outputDirectoryLocation)
                if f.startswith(COMPILED_DATABASE)]

    for f in filelist:
        os.remove(os.path.join(outputDirectoryLocation, f))

    os.remove(os.path.join(outputDirectoryLocation, COMPILED_SIGNATURES))


"""
# =============================================================================

PARSE

# =============================================================================
"""
def parse(parameters):

    signatureLocations = []
    Utility.expandInput(parameters[SIGNATURES], signatureLocations)

    outputDirectoryLocation = parameters[OUTPUT]

    seedSize = parameters[SEED_SIZE] \
        if parameters[SEED_SIZE] else SEED_SIZE_DEFAULT

    consolidateSignatures(
        signatureLocations, seedSize, outputDirectoryLocation)


"""
# =============================================================================

MAIN

# =============================================================================
"""
def main():

    # --- Parser ---
    parser = argparse.ArgumentParser(description=PROGRAM_DESCRIPTION)

    parser.add_argument(
        SIGNATURES_SHORT,
        SIGNATURES_LONG,
        dest=SIGNATURES,
        help=SIGNATURES_HELP,
        type=str, required=True, nargs='+')

    parser.add_argument(
        OUTPUT_SHORT,
        OUTPUT_LONG,
        dest=OUTPUT,
        help=OUTPUT_HELP,
        type=str, required=True)

    parser.add_argument(
        SEED_SIZE_SHORT,
        SEED_SIZE_LONG,
        dest=SEED_SIZE,
        help=SEED_SIZE_HELP,
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
