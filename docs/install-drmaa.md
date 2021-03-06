# DRMAA #

The following is only necessary for execution of Neptune in DRMAA mode on a cluster computing environment. These instructions require a strong understanding of Unix and cluster computing configuration. The user will need to manually install and configure a DRMAA-compliant scheduler, such as [SGE](http://gridscheduler.sourceforge.net/) or [Slurm](http://slurm.schedmd.com/), on either a single machine or on a computing cluster. The user will additionally need to install and configure Python DRMAA bindings with considerations for the DRMAA-compliant scheduler.

## DRMAA Requirements ##

The following are required to operate Neptune in DRMAA mode:

* DRMAA-Compliant Scheduler
* Python DRMAA Bindings

### DRMAA-Compliant Scheduler ###

Neptune has been tested using SGE installed on a single machine with the following instructions:

[https://scidom.wordpress.com/2012/01/18/sge-on-single-pc/](https://scidom.wordpress.com/2012/01/18/sge-on-single-pc/)

Furthermore, Neptune has been tested using Slurm on a large computing cluster. However, any DRMAA-compliant scheduler is expected to work. The instructions for installing and configuring such scheduling environments are beyond the scope of this resource.

### Python DRMAA Bindings ###

Neptune uses a Python DRMAA binding to schedule DRMAA jobs and communicate with the scheduler. The information necessary for installing and configuring the Python DRMAA bindings is available the following location:

[https://github.com/pygridtools/drmaa-python](https://github.com/pygridtools/drmaa-python)

## DRMAA Installation ##

It may be helpful to create a submission wrapper script for Neptune to avoid entering the same DRMAA native specification parameters for every submission. The following SGE and Slurm submission wrapper scripts automatically include native specification parameters, appropriate for the scheduling environment, which may be overwritten by the submitting user as necessary.

### Slurm ###

#### Wrapper ####

```bash
#!/usr/bin/env bash

DRMAA_LIBRARY_PATH=/usr/local/lib/libdrmaa.so.1

neptune --drmaa --default-specification "-n 1 --nodes=1 --ntasks-per-node=1 --mem=10240" $@
```

#### Example ####

```bash
neptune-slurm -i /path/to/inclusion/ -e /path/to/exclusion/ -o /path/to/output/
```

### SGE ###

#### Wrapper ####

```bash
#!/usr/bin/env bash

DRMAA_LIBRARY_PATH=/opt/gridengine/lib/linux-x64/libdrmaa.so

neptune --drmaa --default-specification "-l h_vmem=8G -pe smp 4" $@
```

#### Example ####

```bash
neptune-sge -i /path/to/inclusion/ -e /path/to/exclusion/ -o /path/to/output/
```
