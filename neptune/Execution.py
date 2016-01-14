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

import os
import sys
import math
from scipy.misc import comb

import Neptune
import Utility
import JobManager

"""
# =============================================================================

GLOBALS

# =============================================================================
"""

EXPECTED_HITS_THRESHOLD = 0.05


"""
# =============================================================================

EXECUTION

# =============================================================================
"""
class Execution():

    def __init__(self, session, args):

        # -- rate --
        # 0.0 <= q <= 1.0
        if (args.rate is not None and
                (float(args.rate) < 0.0 or
                    float(args.rate) > 1.0)):
            raise RuntimeError("The rate is out of range.")

        self.rate = args.rate

        # -- minimum inclusion hits --
        # 1 <= inhits
        if (args.inhits is not None and
                (int(args.inhits) < 1)):
            raise RuntimeError("The inclusion hits is out of range.")

        self.inhits = args.inhits

        # -- minimum exclusion hits --
        # 1 <= exhits
        if (args.exhits is not None and
                (int(args.exhits) < 1)):
            raise RuntimeError("The exclusion hits is out of range.")

        self.exhits = args.exhits

        # -- maximum gap size --
        # 1 <= gap
        if (args.gap is not None and
                (int(args.gap) < 1)):
            raise RuntimeError("The gap size is out of range.")

        self.gap = args.gap

        # -- minimum signature size --
        # 1 <= size
        if (args.size is not None and
                (int(args.size) < 1)):
            raise RuntimeError("The signature size is out of range.")

        self.size = args.size

        # -- GC-content --
        # 0.0 <= gc <= 1.0
        if (args.gcContent is not None and
                (float(args.gcContent) < 0.0 or
                    float(args.gcContent) > 1.0)):
            raise RuntimeError("The GC-content is out of range.")

        self.gcContent = args.gcContent

        # -- statistical confidence --
        # 0.0 < confidence < 1.0
        if (args.confidence is not None and
                (float(args.confidence) <= 0.0 or
                    float(args.confidence) >= 1.0)):
            raise RuntimeError("The statistical confidence is out of range.")

        self.confidence = args.confidence

        # -- filter length --
        # 0.0 <= filterLength <= 1.0
        if (args.filterLength is not None and
                (float(args.filterLength) < 0.0 or
                    float(args.filterLength) > 1.0)):
            raise RuntimeError("The filter length is out of range.")

        self.filterLength = args.filterLength

        # -- filter percent --
        # 0.0 <= filterPercent <= 1.0
        if (args.filterPercent is not None and
                (float(args.filterPercent) < 0.0 or
                    float(args.filterPercent) > 1.0)):
            raise RuntimeError("The filter percent is out of range.")

        self.filterPercent = args.filterPercent

        # -- seed size --
        # 4 <= seedSize
        if (args.seedSize is not None and
                (int(args.seedSize) < 4)):
            raise RuntimeError("The seed size is out of range.")

        self.seedSize = args.seedSize

        # -- parallelization --
        # 1 <= parallelization
        if (args.parallelization is not None and
                (int(args.parallelization) < 0)):
            raise RuntimeError("The parallelization is out of range.")

        self.parallelization = args.parallelization

        # -- inclusion locations --
        # inclusion exists
        if args.inclusion is None:
            raise RuntimeError("Inclusion sequence(s) are missing.")

        self.inclusionLocations = []
        Utility.expandInput(args.inclusion, self.inclusionLocations)

        if len(args.inclusion) is 0:
            raise RuntimeError("Inclusion sequence(s) are missing.")

        # -- exclusion locations --
        # exclusion exists
        if args.exclusion is None:
            raise RuntimeError("Exclusion sequence(s) are missing.")

        self.exclusionLocations = []
        Utility.expandInput(args.exclusion, self.exclusionLocations)

        if len(args.exclusion) is 0:
            raise RuntimeError("exclusion sequence(s) are missing.")

        # -- reference locations --
        self.reference = args.reference

        # -- reference size --
        self.referenceSize = args.referenceSize

        # -- output locations --
        # output exists
        if args.output is None:
            raise RuntimeError("The output directory is missing.")

        # -- k-mer --
        # 1 <= k
        if (args.kmer is not None and
                (int(args.kmer) < 1)):
            raise RuntimeError("The k-mer size is out of range.")

        elif args.kmer is not None:
            self.k = int(args.kmer)

        else:
            self.estimateKMerSize()

        self.outputDirectoryLocation = os.path.abspath(args.output)

        if not os.path.exists(self.outputDirectoryLocation):
            os.makedirs(self.outputDirectoryLocation)

        self.candidatesDirectoryLocation = os.path.abspath(
            os.path.join(self.outputDirectoryLocation, Neptune.CANDIDATES))
        if not os.path.exists(self.candidatesDirectoryLocation):
            os.makedirs(self.candidatesDirectoryLocation)

        self.filteredDirectoryLocation = os.path.abspath(
            os.path.join(self.outputDirectoryLocation, Neptune.FILTERED))
        if not os.path.exists(self.filteredDirectoryLocation):
            os.makedirs(self.filteredDirectoryLocation)

        self.sortedDirectoryLocation = os.path.abspath(
            os.path.join(self.outputDirectoryLocation, Neptune.SORTED))
        if not os.path.exists(self.sortedDirectoryLocation):
            os.makedirs(self.sortedDirectoryLocation)

        self.consolidatedDirectoryLocation = os.path.abspath(
            os.path.join(self.outputDirectoryLocation, Neptune.CONSOLIDATED))
        if not os.path.exists(self.consolidatedDirectoryLocation):
            os.makedirs(self.consolidatedDirectoryLocation)

        self.databaseDirectoryLocation = os.path.abspath(
            os.path.join(self.outputDirectoryLocation, Neptune.DATABASE))
        if not os.path.exists(self.databaseDirectoryLocation):
            os.makedirs(self.databaseDirectoryLocation)

        self.kmersOutputDirectory = os.path.abspath(
            os.path.join(self.outputDirectoryLocation, Neptune.KMERS))
        self.inclusionOutputDirectory = os.path.abspath(
            os.path.join(self.kmersOutputDirectory, Neptune.INCLUSION))
        self.exclusionOutputDirectory = os.path.abspath(
            os.path.join(self.kmersOutputDirectory, Neptune.EXCLUSION))

        self.aggregateLocation = os.path.abspath(
            os.path.join(self.outputDirectoryLocation, Neptune.AGGREGATE))

        self.logDirectoryLocation = os.path.abspath(
            os.path.join(self.outputDirectoryLocation, Neptune.LOG))
        if not os.path.exists(self.logDirectoryLocation):
            os.makedirs(self.logDirectoryLocation)

        # -- job manager --
        self.jobManager = JobManager.JobManager(
            session, self.outputDirectoryLocation,
            self.logDirectoryLocation, args.defaultSpecification)

        # -- job specifications --
        if args.countSpecification:
            self.jobManager.setCountSpecification(
                args.countSpecification)

        if args.aggregateSpecification:
            self.jobManager.setAggregateSpecification(
                args.aggregateSpecification)

        if args.extractSpecification:
            self.jobManager.setExtractSpecification(
                args.extractSpecification)

        if args.databaseSpecification:
            self.jobManager.setDatabaseSpecification(
                args.databaseSpecification)

        if args.filterSpecification:
            self.jobManager.setFilterSpecification(
                args.filterSpecification)

        if args.consolidateSpecification:
            self.jobManager.setConsolidateSpecification(
                args.consolidateSpecification)

    """
    # =========================================================================

    CALCULATE EXPECTED K-MER HITS

    PURPOSE:
        Calculates the expected number of arbitrary k-mer matches with a
        genome.

        This is calculating P(k_x = k_y) * ((gs - k + 1) C (2)).

    INPUT:
        [0 <= FLOAT <= 1] [GC] - The GC content.
        [0 <= FLOAT <= 1] [GS] - The genome size.
        [1 <= INT] [K] - The k-mer size.

    RETURN:
        [FLOAT >= 0] [expected] - The expected number of arbitrary k-mer
            matches.

    # =========================================================================
    """
    @staticmethod
    def calculateExpectedKMerHits(gc, gs, k):

        # P(k_x = k_y) * ((gs - k + 1) C (2)) -- from manuscript
        a = 2.0 * math.pow((1.0 - gc) / 2.0, 2.0)
        b = 2.0 * math.pow(gc / 2.0, 2.0)
        c = math.pow(a + b, k)
        d = comb((gs - k + 1), (2), False)
        expected = c * d

        return expected

    """
    # =========================================================================

    ESTIMATE K-MER SIZE

    PURPOSE:
        Estimates the appropriate k-mer size.

    INPUT:
        [NONE]

    RETURN:
        [NONE]

    POST:
        The member variable [self.k] is assigned an appriorate value of k,
        determined by the most extreme GC-content of all references and the
        largest size of any reference.
    # =========================================================================
    """
    def estimateKMerSize(self):

        print "Estimating k-mer size ..."

        maxGenomeSize = 1
        maxGCContent = 0.5  # least extreme GC-content

        for inclusionLocation in self.inclusionLocations:

            inclusionFile = open(inclusionLocation, 'r')

            size = 0
            sumGC = 0
            sumAT = 0
            gcContent = 0

            for line in inclusionFile:

                if line[0] != ">":
                    line = line.strip()

                    size += len(line)
                    sumGC += line.count('G') + line.count('C')
                    sumAT += line.count('A') + line.count('T')

            gcContent = float(sumGC) / float(sumGC + sumAT)

            # SIZE
            if size > maxGenomeSize:
                maxGenomeSize = size

            # GC-CONTENT
            if gcContent > maxGCContent:
                maxGCContent = gcContent

            elif (1.0 - gcContent) > maxGCContent:
                maxGCContent = (1.0 - gcContent)

            inclusionFile.close()

        """
        NOTE:

        When GC-content is 1.0 (worst),
        and the length of genome is 10^80 (atoms in observable universe),
        and the threshold is 0.05,
        then the k required is 535.
        Therefore, k = 535 is the maximum we use.

        """
        for k in range(3, 535, 2):

            expected = self.calculateExpectedKMerHits(
                maxGCContent, maxGenomeSize, k)

            if expected < EXPECTED_HITS_THRESHOLD:

                self.k = k
                print "k = " + str(self.k)
                return

        # No suitable k estimated.
        raise RuntimeError("ERROR: No suitable value for k determined. \n")

    """
    # =========================================================================

    PRODUCE RECEIPT

    PURPOSE:
        Produces an execution receipt. This receipt provides information about
        the execution parameters.

    INPUT:
        [NONE]

    RETURN:
        [NONE]

    POST:
        An execution receipt is printed to standard output.

    # =========================================================================
    """
    def produceReceipt(self):

        receiptLocation = os.path.abspath(
            os.path.join(self.outputDirectoryLocation, Neptune.RECEIPT))
        receiptFile = open(receiptLocation, "w")

        receiptFile.write(
            "==============================================================\n")
        receiptFile.write("RUN RECEIPT\n")
        receiptFile.write(
            "==============================================================\n")
        receiptFile.write("\n")

        self.reportCommandLine(receiptFile)
        self.reportFiles(receiptFile)
        self.reportGeneralParameters(receiptFile)
        self.reportDRMAAParameters(receiptFile)

        receiptFile.close()

    """
    # =========================================================================

    REPORT COMMAND LINE

    PURPOSE:
        Reports the command line arguments to the execution receipt.

    INPUT:
        [FILE] [receiptFile] - The open and writable receipt file.

    RETURN:
        [NONE]

    POST:
        The command line arguments will be reported to the execution receipt.

    # =========================================================================
    """
    def reportCommandLine(self, receiptFile):

        receiptFile.write("-- Command Line -- \n")
        receiptFile.write("\n")

        for arg in sys.argv:
            receiptFile.write(str(arg) + " ")

        receiptFile.write("\n")
        receiptFile.write("\n")

    """
    # =========================================================================

    REPORT GENERAL PARAMETERS

    PURPOSE:
        Reports the general parameters used by Neptune.

    INPUT:
        [FILE] [receiptFile] - The open and writable receipt file.

    RETURN:
        [NONE]

    POST:
        The general parameters used by Neptune will be reported to the
        execution receipt.

    # =========================================================================
    """
    def reportGeneralParameters(self, receiptFile):

        receiptFile.write("-- General -- \n")
        receiptFile.write("\n")

        receiptFile.write(
            "k = "
            + str(self.k) + "\n")

        receiptFile.write(
            "SNV Rate = "
            + str(self.rate) + "\n")

        receiptFile.write(
            "Minimum Inclusion Observations = "
            + str(self.inhits) + "\n")

        receiptFile.write(
            "Minimum Exclusion Observations = "
            + str(self.exhits) + "\n")

        receiptFile.write(
            "Maximum Gap Size = "
            + str(self.gap) + "\n")

        receiptFile.write(
            "Minimum Signature Size = "
            + str(self.size) + "\n")

        receiptFile.write(
            "GC-Content = "
            + str(self.gcContent) + "\n")

        receiptFile.write(
            "Filter Length = "
            + str(self.filterLength) + "\n")

        receiptFile.write(
            "Filter Percent = "
            + str(self.filterPercent) + "\n")

        receiptFile.write(
            "Parallelization = "
            + str(self.parallelization) + "\n")

        receiptFile.write(
            "Reference Size = "
            + str(self.referenceSize) + "\n")

        receiptFile.write("\n")

    """
    # =========================================================================

    REPORT FILES

    PURPOSE:
        Reports the files used by Neptune.

    INPUT:
        [FILE] [receiptFile] - The open and writable receipt file.

    RETURN:
        [NONE]

    POST:
        The files used by Neptune will be reported to the execution receipt.

    # =========================================================================
    """
    def reportFiles(self, receiptFile):

        receiptFile.write("-- Files -- \n")
        receiptFile.write("\n")

        receiptFile.write("Inclusion Targets = \n")

        for location in self.inclusionLocations:
            receiptFile.write("\t" + str(location) + "\n")

        receiptFile.write("Exclusion Targets = \n")

        for location in self.exclusionLocations:
            receiptFile.write("\t" + str(location) + "\n")

        if self.reference:

            receiptFile.write("References = \n")

            for ref in self.reference:
                receiptFile.write("\t" + str(ref) + "\n")

        else:

            receiptFile.write("Reference = " + str(self.reference) + "\n")

        receiptFile.write(
            "Output = \n"
            + str("\t" + self.outputDirectoryLocation) + "\n")

        receiptFile.write("\n")

    """
    # =========================================================================

    REPORT DRMAA PARAMETERS

    PURPOSE:
        Reports the DRMAA parameters used by Neptune.

    INPUT:
        [FILE] [receiptFile] - The open and writable receipt file.

    RETURN:
        [NONE]

    POST:
        The DRMAA parameters used by Neptune will be reported to the execution
        receipt.

    # =========================================================================
    """
    def reportDRMAAParameters(self, receiptFile):

        receiptFile.write("-- DRMAA -- \n")
        receiptFile.write("\n")

        receiptFile.write(
            "CountKMers Specification = "
            + str(self.jobManager.countSpecification) + "\n")

        receiptFile.write(
            "AggregateKMers Specification = "
            + str(self.jobManager.aggregateSpecification) + "\n")

        receiptFile.write(
            "ExtractSignatures Specification = "
            + str(self.jobManager.extractSpecification) + "\n")

        receiptFile.write(
            "CreateDatabase Specification = "
            + str(self.jobManager.databaseSpecification) + "\n")

        receiptFile.write(
            "FilterSignatures Specification = "
            + str(self.jobManager.filterSpecification) + "\n")

        receiptFile.write(
            "ConsolidateSignatures Specification = "
            + str(self.jobManager.consolidateSpecification) + "\n")

        receiptFile.write("\n")
