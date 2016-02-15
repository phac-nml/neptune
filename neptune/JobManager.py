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

This file contains the JobManager class. JobManager is responsible for
managing the creation and execution of DRMAA jobs.

# =============================================================================
"""

import drmaa
import os
import sys
import inspect

import CountKMers
import AggregateKMers
import ExtractSignatures
import FilterSignatures
import ConsolidateSignatures

"""
# =============================================================================

JOB MANAGER

# =============================================================================
"""

class JobManager:

    NEPTUNE_JOB = "Neptune"
    COUNT_JOB = "Neptune-CountKMers"
    AGGREGATE_JOB = "Neptune-AggregateKMers"
    EXTRACT_JOB = "Neptune-ExtractSignatures"
    DATABASE_JOB = "Neptune-CreateDatabase"
    FILTER_JOB = "Neptune-FilterSignatures"
    CONSOLIDATE_JOB = "Neptune-ConsolidateSignatures"

    ID = 0

    """
    # =========================================================================

    CONSTRUCTOR

    INPUT:
        [FILE LOCATION] [outputDirectoryLocation] - The directory location to
            write DRMAA output.
        [FILE LOCATION] [logDirectoryLocation] - The directory location to
            write DRMAA output logs and error logs.

    # =========================================================================
    """
    def __init__(
            self, outputDirectoryLocation, logDirectoryLocation):

        self.outputDirectoryLocation = outputDirectoryLocation
        self.logDirectoryLocation = logDirectoryLocation

    """
    # =========================================================================

    GENERATE ID

    PURPOSE:
        Generates an ID unique to this JobManager. This unique ID DOES NOT
        correspond with the DRMAA job ID.

    INPUT:
        [NONE]

    RETURN:
        [INT] [ID] - The generated unique ID.

    POST:
        [NONE]

    # =========================================================================
    """
    def generateID(self):

        self.ID += 1

        return int(self.ID)
