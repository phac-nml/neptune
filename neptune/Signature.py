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
    Reads a signature file into a signature dictionary.

INPUT:
    [STRING] [fileLocation] - The file location of the signatures.

RETURN:
    [(STRING) -> (SIGNATURE) DICTIONARY] - A dictionary mapping string IDs
        to signature objects.

# =========================================================================
"""
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


INPUT:
    [SIGNATURE LIST] [signatures] -


RETURN:


POST:

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
    [WRITABLE] [destination] - An open and writable object.

RETURN:
    [NONE]

POST:
    The passed signature will be written to the destination.

# =========================================================================
"""
def writeSignature(signature, destination):

    destination.write(
        ">" + str(signature.ID) + " "
        + str("score=") + str(signature.score) + " "
        + str("in=") + str(signature.inscore) + " "
        + str("ex=") + str(signature.exscore) + " "
        + str("len=") + str(signature.length) + " "
        + str("ref=") + str(signature.reference) + " "
        + str("pos=") + str(signature.position) + "\n")

    destination.write(str(signature.sequence) + "\n")

"""
# =============================================================================

SORT SIGNATURES

PURPOSE:


INPUT:


RETURN:


POST:

# =============================================================================
"""
def sortSignatures(signatures):

    sortedSignatures = sorted(
        signatures.iteritems(),
        key=lambda (k, v): v.score, reverse=True)

    return [item[1] for item in sortedSignatures]
