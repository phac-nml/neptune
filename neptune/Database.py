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

CREATE DATABASE JOB

PURPOSE:
    Creates a build database job.

INPUT:
    [STRING ITERATOR] [inputLocations] - The location of a single FASTA file
        from which to build the database.
    [STRING] [outputLocation] - The output location of the database.

RETURN:
    [DRMAA JOB TEMPLATE] [job] - A create database job.

# =============================================================================
"""
def createDatabaseJob(inputLocation, outputLocation):

    # COMMAND LINE
    COMMAND = "makeblastdb"

    TYPE = "-dbtype"
    NUCLEOTIDE = "nucl"

    INPUT = "-in"
    INPUT_LOCATIONS = inputLocation

    TITLE = "-title"
    NAME = "DATABASE"

    OUTPUT = "-out"
    OUTPUT_LOCATION = outputLocation

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

"""
# =============================================================================

QUERY DATABASE

PURPOSE:
    Queries the database with a specified query.

INPUT:
    [FILE LOCATION] [databaseLocation] - The file location of the database.
    [FILE LOCATION] [queryLocation] - The file location of the query.
    [FILE LOCATION] [outputLocation] - The file location of the output.
    [0 <= FLOAT <= 1] [percentIdentity] - The minimum percent identity of an
        alignment.

RETURN:
    [FILE LOCATION] [outputLocation] - The file location of the query output.

POST:
    A query file will be created in the output directory.

# =============================================================================
"""
def queryDatabase(
        databaseLocation, queryLocation, outputLocation, percentIdentity):

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

    args = [
        COMMAND,
        DATABASE, databaseLocation,
        QUERY, queryLocation,
        OUTPUT, outputLocation,
        OUTPUT_FORMAT, OUTPUT_FORMAT_STRING,
        PERCENT_IDENTITY, str(percentIdentity),
        WORD_SIZE, str(WORD_SIZE_VALUE),
        DUST, DUST_VALUE]

    subprocess.check_call(args)

    return outputLocation
