# Introduction #

A genomic signature is a genomic locus that is sufficiently represented in an 
inclusion group, and sufficiently absent from a background, or exclusion 
group. A signature might correlate genomic features with phenotypic traits, 
such as the presence of a gene with increased organism pathogenicity.

Neptune locates genomic signatures using an exact k-mer matching strategy while
accommodating k-mer mismatches. The software identifies sequences that are
sufficiently represented within inclusion targets and sufficiently absent from
exclusion targets. The signature discovery process is accomplished using
probabilistic models instead of heuristic strategies. Neptune may be leveraged 
to reveal discriminatory signature sequences to uniquely delineate one group of 
organisms, such as isolates associated with a disease cluster or event, from 
unrelated sporadic or environmental microbes.

# Resources #

* **Website**: [https://phac-nml.github.io/neptune/](https://phac-nml.github.io/neptune/)
* **Installation**: [https://phac-nml.github.io/neptune/install/](https://phac-nml.github.io/neptune/install/)
* **Walkthrough**: [https://phac-nml.github.io/neptune/walkthrough/](https://phac-nml.github.io/neptune/walkthrough/)

# Release #

**Neptune 1.2.4**

This release makes several small improvements, including: reducing the standard
output clutter, adding timings to stages, and updating the documentation.

# Requirements #

Neptune requires Python 2.7. You may check your installed Python version with 
the following:

        python --version

If running a Debian distribution (ex: Ubuntu), dependencies may be installed
using the following command:

        sudo install/debian_dependencies.sh

Otherwise, the following dependencies must be installed manually:

- python-pip
- python-virtualenv
- build-essential
- python-dev
- NCBI BLAST+

# Installation #

It is strongly recommended you refer to the
[documentation](https://phac-nml.github.io/neptune/install/) for full 
installation instructions. Neptune may be installed using the following 
command:

    INSTALL.sh

You may specify an install PREFIX location, and Neptune will install into
PREFIX/lib and PREFIX/bin. This only requires security privileges if the
install location requires them.

# Running Neptune #

Neptune's command line arguments can be found by running:

    neptune --help

A simple example of running Neptune:

    neptune --inclusion /path/to/inclusion/ --exclusion /path/to/exclusion/
            --output /path/to/output/

Please refer to the 
[documentation](https://phac-nml.github.io/neptune/parameters/) for more 
details.

# Contact #

**Eric Marinier**: eric.marinier@phac-aspc.gc.ca

# Legal #

Neptune

Copyright Government of Canada 2015-2017

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

