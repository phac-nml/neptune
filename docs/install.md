# Installation #

This installation guide assumes the use of the [BASH](https://en.wikipedia.org/wiki/Bash_(Unix_shell)) Unix shell. Neptune may be installed on most 64-bit Unix environments. However, the specifics of installations on all environments is beyond the scope of this documentation. Neptune may either be run on a single machine or a computing cluster. Neptune achieves maximum parallelization when submitting jobs through a DRMAA-compliant cluster computing scheduler. The installation and configuration of a DRMAA-compliant scheduler will require a significant understanding of Unix. However, it is possible to run Neptune in parallel on a single machine without DRMAA. Neptune is known to be compatible with the [SGE](http://gridscheduler.sourceforge.net/) and [Slurm](http://slurm.schedmd.com/) schedulers.

## Python ##

Neptune requires Python 2.7. Note that Python 2.7 is provided with many major distributions of Linux. The following may check your Python version:

    $ python --version

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

## DRMAA Requirements ##

The following is only necessary for execution of Neptune in DRMAA mode on a cluster computing environment. These instructions require a strong understanding of Unix and cluster computing configuration. The user will need to manually install and configure a DRMAA-compliant scheduler, such as [SGE](http://gridscheduler.sourceforge.net/) or [Slurm](http://slurm.schedmd.com/), on either a single machine or on a computing cluster. The user will additionally need to install and configure Python DRMAA bindings with considerations for the DRMAA-compliant scheduler. The following are required to operate Neptune in DRMAA mode:

* DRMAA-compliant scheduler
* Python DRMAA bindings

### DRMAA-Compliant Scheduler ###

Neptune has been tested using SGE installed on a single machine with the following instructions:

    https://scidom.wordpress.com/2012/01/18/sge-on-single-pc/

However, any DRMAA-compliant scheduler is expected to work. The instructions for installing and configuring such scheduling environments are beyond the scope of this resource.

### Python DRMAA Bindings ###

Neptune uses a Python DRMAA binding to schedule DRMAA jobs and communicate with the scheduler. The information necessary for installing and configuring the Python DRMAA bindings is available the following location:

    https://github.com/pygridtools/drmaa-python

## DRMAA Installation ##

It may be helpful to create a submission wrapper script for Neptune to avoid entering the same DRMAA native specification parameters for every submission. The following SGE and Slurm submission wrapper scripts automatically include native specification parameters, appropriate for the scheduling environment, which may be overwritten by the submitting user as necessary.

### Slurm Wrapper ###

    #!/usr/bin/env bash
    
    DRMAA_LIBRARY_PATH=/usr/local/lib/libdrmaa.so.1

    neptune --drmaa --default-specification "-n 1 --nodes=1 --ntasks-per-node=1 --mem=10240" $@

### Slurm Example ###

    neptune-slurm -i /path/to/inclusion/ -e /path/to/exclusion/ -o /path/to/output/

### SGE Wrapper ###

    #!/usr/bin/env bash

    DRMAA_LIBRARY_PATH=/opt/gridengine/lib/linux-x64/libdrmaa.so

    neptune --drmaa --default-specification "-l h_vmem=8G -pe smp 4" $@

### SGE Example ###

    neptune-sge -i /path/to/inclusion/ -e /path/to/exclusion/ -o /path/to/output/

