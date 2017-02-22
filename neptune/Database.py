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
    Creates a build database job that will execute a Python subprocess and
    construct the database using BLAST.

INPUT:
    [FILE LOCATION] [inputLocation] - The location of a single FASTA file
        from which to build the database.
    [FILE LOCATION] [outputLocation] - The output location of the database.

POST:
    A Python subprocess is executed for the create database job. Control will
    return to the calling function after this subprocess is complete.

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

    subprocess.check_call(args)


"""
# =============================================================================

QUERY DATABASE

PURPOSE:
    Queries the database with a specified query by executing a Python
    subprocess.

INPUT:
    [FILE LOCATION] [databaseLocation] - The file location of the database.
    [FILE LOCATION] [queryLocation] - The file location of the query (FASTA).
    [FILE LOCATION] [outputLocation] - The file location to write the output.
    [0 <= FLOAT <= 1] [percentIdentity] - The minimum percent identity of an
        alignment for it to be reported.
    [4 <= INT] [seedSize] - The seed size used in query alignments.

RETURN:
    [FILE LOCATION] [outputLocation] - The file location of the query output.
        This is the same location as the passed [outputLocation]

POST:
    A query file will be created at the [outputLocation].

# =============================================================================
"""
def queryDatabase(
        databaseLocation, queryLocation, outputLocation,
        percentIdentity, seedSize):

    # COMMAND LINE
    COMMAND = "blastn"

    DATABASE = "-db"
    QUERY = "-query"
    OUTPUT = "-out"
    OUTPUT_FORMAT = "-outfmt"
    OUTPUT_FORMAT_STRING = "6 qseqid qlen sseqid length pident score"
    PERCENT_IDENTITY = "-perc_identity"
    WORD_SIZE = "-word_size"
    WORD_SIZE_VALUE = seedSize
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
