# Neptune #
-----------

## Legal ##
-----------

Neptune

Copyright Government of Canada 2015

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

**Neptune 1.1.1**

2016 January 19

This release of Neptune updates the installation instructions to be more
informative.

## Requirements ##
------------------

Neptune will attempt to automatically install several dependencies. However,
the following dependencies must be installed by the user:

- Python 2.7
- NumPy
- SciPy
- DRMAA-compliant scheduler
- Python DRMAA
- BLAST+
- pipsi

## Installation ##
------------------

**
It is strongly recommended you refer to the manual for full installation
instructions. The following description is an abridged version of the
installation process.
**

Install and configure a DRMAA-compliant scheduler and Python DRMAA bindings.
This process may require considerable time. Please refer to the
instructions available at:

	https://github.com/pygridtools/drmaa-python

The pipsi tool is installed using the following command:

	curl https://raw.githubusercontent.com/mitsuhiko/pipsi/master/get-pipsi.py | python

You may need to add the pipsi install location to your PATH variable. For BASH
shells, this can be done by adding the following line to your bashrc (usually
~/.bashrc) file:

	export PATH=$PATH:~/.local/bin

You will then need to source your bashrc file:

	source ~/.bashrc

You may use the following command to ensure the PATH was updated correctly. It
should display your PATH with "local/bin" highlighted. If nothing is displayed,
then you have not updated your PATH correctly.

	echo $PATH | grep local/bin

Install Neptune using the pipsi tool:

	pipsi install /path/to/neptune/download/location/


## Running Neptune ##
---------------------

Neptune's command line arguments can be found by running:

	neptune --help

A simple example of running Neptune:

	neptune --inclusion /path/to/inclusion/ --exclusion /path/to/exclusion/
		--output /path/to/output/ --parallelization 3

Please refer to the documentation for more details.

## Contact ##
-------------

**Eric Marinier**: eric.marinier@phac-aspc.gc.ca

