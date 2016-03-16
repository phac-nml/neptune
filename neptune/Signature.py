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

SIGNATURE

# =============================================================================
"""
class Signature():

    def __init__(
            self, ID, score, inscore, exscore, sequence, reference, position):

        self.ID = str(ID)
        self.score = float(score)
        self.inscore = float(inscore)
        self.exscore = float(exscore)
        self.sequence = str(sequence).strip()
        self.length = len(str(sequence).strip())
        self.reference = str(reference)
        self.position = int(position)

"""
# =========================================================================

READ SIGNATURES

PURPOSE:
    Reads a signature file and places the signatures into a new signature
    dictionary. This function is designed to be symmetric with the write
    signatures function.

INPUT:
    [FILE LOCATION] [fileLocation] - The file location of the signatures.

RETURN:
    [(STRING ID) -> (SIGNATURE) DICTIONARY] - A dictionary mapping string
        IDs to signature objects.

# =========================================================================
"""
def readSignatures(fileLocation):

    signaturesFile = open(fileLocation, 'r')
    signatures = {}

    while True:

        # read lines
        line1 = signaturesFile.readline()
        line2 = signaturesFile.readline()

        # reached the end of file
        if not line2:
            break

        tokens = (line1[1:]).split()

        ID = tokens[0]
        score = tokens[1].split("=")[1]
        inscore = tokens[2].split("=")[1]
        exscore = tokens[3].split("=")[1]
        reference = tokens[5].split("=")[1]
        position = tokens[6].split("=")[1]

        sequence = line2

        signature = Signature(
            ID, score, inscore, exscore, sequence, reference, position)
        signatures[ID] = signature

    signaturesFile.close()

    return signatures

"""
# =============================================================================

WRITE SIGNATURES

PURPOSE:
    Writes several signatures to the destination. This function is designed to
    be symmetric with the read signatures function.

INPUT:
    [SIGNATURE LIST] [signatures] - A list of Signature objects.
    [FILE] [destination] - An open and writable file-like object.

RETURN:
    [NONE]

POST:
    The [signatures] will be written to the [destination]

# =============================================================================
"""
def writeSignatures(signatures, destination):

    for signature in signatures:

        writeSignature(signature, destination)

"""
# =========================================================================

WRITE SIGNATURE

PURPOSE:
    Writes the signature to the destination. This function is designed to
    be symmetric with the read signatures function.

INPUT:
    [SIGNATURE] [signature] - The signature to write.
    [FILE] [destination] - An open and writable file-like object.

RETURN:
    [NONE]

POST:
    The passed signature will be written to the [destination].

# =========================================================================
"""
def writeSignature(signature, destination):

    destination.write(
        ">" + str(signature.ID) + " " +
        str("score=") + "{0:.4f}".format(signature.score) + " " +
        str("in=") + "{0:.4f}".format(abs(signature.inscore)) + " " +
        str("ex=") + "{0:.4f}".format(abs(signature.exscore)) + " " +
        str("len=") + str(signature.length) + " " +
        str("ref=") + str(signature.reference) + " " +
        str("pos=") + str(signature.position) + "\n")

    destination.write(str(signature.sequence) + "\n")

"""
# =============================================================================

SORT SIGNATURES

PURPOSE:
    Sorts the signatures by their corresponding score in descending order.

INPUT:
    [(STRING ID) -> (SIGNATURE) DICTIONARY] [signatures] - A dictionary mapping
        signature IDs to their corresponding signature.

RETURN:
    [SIGNATURE LIST] [sortedSignatures] - A list of sorted Signature objects.

POST:
    [NONE]

# =============================================================================
"""
def sortSignatures(signatures):

    sortedSignatures = sorted(
        signatures.iteritems(),
        key=lambda (k, v): v.score, reverse=True)

    return [item[1] for item in sortedSignatures]
