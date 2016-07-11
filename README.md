# Neptune #
-----------

## Legal ##
-----------

Neptune

Copyright Government of Canada 2015-2016

Written by: Eric Marinier, Public Health Agency of Canada,
    National Microbiology Laboratory

Funded by the National Micriobiology Laboratory and the Genome Canada / Alberta
    Innovates Bio Solutions project "Listeria Detection and Surveillance
    using Next Generation Genomics"

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

## Description ##
-----------------

Neptune locates genomic signatures using an exact k-mer matching strategy while
accommodating k-mer mismatches. The software identifies sequences that are
sufficiently represented within inclusion targets and sufficiently absent from
exclusion targets. The signature discovery process is accomplished using
probabilistic models instead of heuristic strategies.

## Release ##
-------------

**Neptune 1.2.3**

2016 July 11

This release simplifies the installation process.

## Requirements ##
------------------

Neptune requires Python 2.7. You may check the installed version with the
following:

	$ python --version

If running a Debian distribution (ex: Ubuntu), dependencies may be installed
using the following command:

	$ sudo install/debian_dependencies.sh

Otherwise, the following dependencies must be installed manually:

- python-pip
- python-virtualenv
- build-essential
- python-dev
- NCBI BLAST+

## Installation ##
------------------

It is strongly recommended you refer to the manual for full installation
instructions. Neptune may be installed using the following command:

        $ INSTALL.sh

You may specify an install PREFIX location, and Neptune will install into
PREFIX/lib and PREFIX/bin. This only requires security privilages if the
install location requires them.

## Running Neptune ##
---------------------

Neptune's command line arguments can be found by running:

	$ neptune --help

A simple example of running Neptune:

	$ neptune --inclusion /path/to/inclusion/ --exclusion /path/to/exclusion/
		--output /path/to/output/

Please refer to the documentation for more details.

## Contact ##
-------------

**Eric Marinier**: eric.marinier@phac-aspc.gc.ca

