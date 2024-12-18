#!/usr/bin/env python

"""
# =============================================================================

Copyright Government of Canada 2015-2024

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

def prepareSystemPath():

    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    print(root)

    if root not in sys.path:
         sys.path.insert(0, root)

def getPath(location):

    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    return os.path.join(root, location)
