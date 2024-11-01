# Introduction #

A genomic signature is a genomic locus that is sufficiently represented in an inclusion group, and sufficiently absent from a background, or exclusion group. A signature might correlate genomic features with phenotypic traits, such as the presence of a gene with increased organism pathogenicity.

Neptune locates genomic signatures using an exact k-mer matching strategy while accommodating k-mer mismatches. The software identifies sequences that are sufficiently represented within inclusion targets and sufficiently absent from exclusion targets. The signature discovery process is accomplished using probabilistic models instead of heuristic strategies. Neptune may be leveraged to reveal discriminatory signature sequences to uniquely delineate one group of organisms, such as isolates associated with a disease cluster or event, from unrelated sporadic or environmental microbes.

# Resources #

* **Website**: [https://phac-nml.github.io/neptune/](https://phac-nml.github.io/neptune/)
* **Installation**: [https://phac-nml.github.io/neptune/install/](https://phac-nml.github.io/neptune/install/)
* **Walkthrough**: [https://phac-nml.github.io/neptune/walkthrough/](https://phac-nml.github.io/neptune/walkthrough/)

# Release #

**Neptune 2.0.0 - 2024-10-21**

This release updates Neptune to Python3, removes DRMAA support, fixes a crash when no signatures are produced, and updates the installation process.

# Installation #

## Python 3

Ensure your version of Python is compatible (python>=3.10):

`python --version`

You may wish to use [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html) to create an environment specifically for this purpose:

`conda create --name neptune 'python>=3.10'`  
`conda activate neptune`

## pip

Ensure you can run pip:

`pip --version` or `python -m pip --version`

If pip is unavailable, please refer to [these instructions](https://packaging.python.org/en/latest/tutorials/installing-packages/) on how to install pip.

## BLAST

Neptune requires BLAST to be manually installed and made available as a command-line program:

`sudo apt-get install ncbi-blast+`

You can verify BLAST was installed by ensuring the follow commands are available:

`makeblastdb -h`  
`blastn -h`

## Neptune and Dependencies

After downloading Neptune's source files, you can install Neptune and all of its pip dependencies with the following:

`pip install /path/to/neptune_directory/` or `pip install .`

**CAUTION**: If you attempt `pip install neptune` (not interpreted as a file path), then you'll download a different package that's also named "neptune" that's available directly from pip.

The following packages and their dependencies will be installed:

- numpy
- scipy
- biopython
- neptune

You can verify the installation was successful with the following:

`neptune --version`

And you can test the installation with simple test inputs with the following:

`neptune -i tests/data/example/inclusion/ -e tests/data/example/exclusion/ -o output`

# Running Neptune #

Neptune's command line arguments can be found by running:

    neptune --help

A simple example of running Neptune:

    neptune --inclusion /path/to/inclusion/ --exclusion /path/to/exclusion/
            --output /path/to/output/

Please refer to the [documentation](https://phac-nml.github.io/neptune/parameters/) for more details.

# Contact #

**Eric Marinier**: eric.marinier@phac-aspc.gc.ca

# Legal #

Neptune

Copyright Government of Canada 2015-2024

Written by: Eric Marinier, Public Health Agency of Canada, National Microbiology Laboratory

Funded by the National Micriobiology Laboratory and the Genome Canada / Alberta Innovates Bio Solutions project "Listeria Detection and Surveillance using Next Generation Genomics"

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this work except in compliance with the License. You may obtain a copy of the License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
