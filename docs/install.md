# Installation #

This installation guide assumes the use of the [BASH](https://en.wikipedia.org/wiki/Bash_(Unix_shell)) Unix shell. Neptune may be installed on most 64-bit Unix environments. However, the specifics of installations on all environments is beyond the scope of this documentation. Neptune may either be run on a single machine or a computing cluster. Neptune achieves maximum parallelization when submitting jobs through a DRMAA-compliant cluster computing scheduler. The installation and configuration of a DRMAA-compliant scheduler will require a significant understanding of Unix. However, it is possible to run Neptune in parallel on a single machine without DRMAA. Neptune is known to be compatible with the [SGE](http://gridscheduler.sourceforge.net/) and [Slurm](http://slurm.schedmd.com/) schedulers.

## Python ##

Neptune requires Python 2.7. Note that Python 2.7 is provided with many major distributions of Linux. The following may check your Python version:

    python --version

## Dependencies ##

### Debian-Based Installation ###

This section assumes the user has the [APT](https://help.ubuntu.com/community/AptGet/Howto) package manager. This is common to the [Ubuntu](https://en.wikipedia.org/wiki/Ubuntu_(operating_system)) operating system. However, this section should be compatible with any 64-bit Debian distribution. The following operation will automatically install Neptune's dependencies and require security privileges (sudo) to install the dependencies:

    sudo neptune/install/debian_dependencies.sh

### Manual Installation ###

If you cannot install the dependencies using the above script, the following dependencies must be manually installed, if necessary, by the user:

* pip
* virtualenv
* build-essential
* python-dev
* NCBI BLAST+

## Neptune ##

Neptune will be installed using pip into its own Python virtual environment. The following will install Neptune locally into the source directory and will not require security privileges:

    neptune/INSTALL.sh

Alternatively, you may specify an install location, PREFIX, such as /usr/local/. Neptune will create the directories PREFIX/lib and PREFIX/bin. This may require security privileges:

    neptune/INSTALL.sh PREFIX
