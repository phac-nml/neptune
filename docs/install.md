# Installation

This installation guide assumes the use of the [BASH](https://en.wikipedia.org/wiki/Bash_(Unix_shell)) Unix shell and a 64-bit Linux system.

## Direct

The following instructions describe how to install Neptune directly. These instructions may require administrative privilages. Directly installing Neptune from the source files involves the following:

 1. Installing Python>=3.10
 2. Installing pip
 3. Installing BLAST (aptitude: `sudo apt-get install ncbi-blast+`)
 4. Installing Neptune (`pip install .`)

More detailed instructions are provided below.

### Python 3

Ensure your version of Python is compatible (python>=3.10):

`python --version`

You may wish to use [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html) to create an environment specifically for this purpose:

`conda create --name neptune 'python>=3.10'`  
`conda activate neptune`

### pip

Ensure you can run pip:

`pip --version` or `python -m pip --version`

If pip is unavailable, please refer to [these instructions](https://packaging.python.org/en/latest/tutorials/installing-packages/) on how to install pip.

### BLAST

Neptune requires BLAST to be manually installed and made available as a command-line program:

`sudo apt-get install ncbi-blast+`

You can verify BLAST was installed by ensuring the follow commands are available:

`makeblastdb -h`  
`blastn -h`

### Neptune and Dependencies

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

## Bioconda

Currently, there is no [Bioconda](https://bioconda.github.io/) package for Neptune v2.0.0. If you wish to install Neptune v1.2.5 with Bioconda, please use the following command to create an environment with Neptune installed:

```bash
conda create -n neptune bioconda::neptune -c conda-forge -c default
```

Please note that specifying `bioconda::neptune` is necessary, because otherwise Conda is likely to resolve the name to different software that's also named `neptune`.

You may also wish to review an [older version](https://github.com/phac-nml/neptune/blob/release/1.2/docs/install.md) of these installation instructions for installing Neptune v1.2.5.
