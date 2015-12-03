import subprocess

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
