# Introduction #

A genomic signature is a genomic locus that is sufficiently represented in an inclusion group, and sufficiently absent from a background, or exclusion group. A signature might correlate genomic features with phenotypic traits, such as the presence of a gene with increased organism pathogenicity.

Neptune locates genomic signatures using an exact k-mer matching strategy while accommodating k-mer mismatches. The software identifies sequences that are sufficiently represented within inclusion targets and sufficiently absent from exclusion targets. The signature discovery process is accomplished using probabilistic models instead of heuristic strategies. Neptune may be leveraged to reveal discriminatory signature sequences to uniquely delineate one group of 
organisms, such as isolates associated with a disease cluster or event, from unrelated sporadic or environmental microbes.

# Resources #

* **Website**: [https://phac-nml.github.io/neptune/](https://phac-nml.github.io/neptune/)
* **Installation**: [https://phac-nml.github.io/neptune/install/](https://phac-nml.github.io/neptune/install/)
* **Walkthrough**: [https://phac-nml.github.io/neptune/walkthrough/](https://phac-nml.github.io/neptune/walkthrough/)

# Release #

**Neptune 1.2.5**

This release provides fixes for ambiguous crashes as a consequence of inputs containing no A, C, G, or T characters, and also makes improvements to the code quality.

# Installation #

It is strongly recommended you refer to the
[documentation](https://phac-nml.github.io/neptune/install/) for full 
installation instructions. Neptune may be installed on any 64-bit Linux system using Bioconda, preferably with [Mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) or with the Mamba resolver under Conda:

 1. Install [Bioconda](https://bioconda.github.io/)
 2. Create an environment for Neptune and install within it: `mamba create -n neptune bioconda::neptune -c conda-forge`
 3. Activate the environment: `mamba activate neptune`
 4. Test the install: `neptune --version`

If you run into problems installing Neptune with Conda / Mamba, you can try the following:

 - Modify your `~/.condarc` file to have `channel_priority: flexible`
 - Modify your conda solver within `~/.condarc` to use Mamba: `solver: libmamba`

Neptune may also be installed directly and instructions are available in the
[documentation](https://phac-nml.github.io/neptune/install/).

# Running Neptune #

Neptune's command line arguments can be found by running:

    neptune --help

A simple example of running Neptune:

    neptune --inclusion /path/to/inclusion/ --exclusion /path/to/exclusion/
            --output /path/to/output/

Please refer to the 
[documentation](https://phac-nml.github.io/neptune/parameters/) for more details.

# Contact #

**Eric Marinier**: eric.marinier@phac-aspc.gc.ca

# Legal #

Neptune

Copyright Government of Canada 2015-2017

Written by: Eric Marinier, Public Health Agency of Canada,
    National Microbiology Laboratory

Funded by the National Micriobiology Laboratory and the Genome Canada / Alberta
    Innovates Bio Solutions project "Listeria Detection and Surveillance using Next Generation Genomics"

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this work except in compliance with the License. You may obtain a copy of the License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
