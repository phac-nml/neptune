FROM debian:jessie

USER root
ARG DEBIAN_FRONTEND=noninteractive

#==============================================================================
# BASIC UTILITIES
#==============================================================================
RUN apt-get update -qqy \
  && apt-get -qqy install \
    apt-utils\
    dialog\
    less\
    vim\
    make\
    build-essential\
    git\
    wget\
    man-db\
  && rm -rf /var/lib/apt/lists/*

#==============================================================================
# SLURM
#==============================================================================
RUN apt-get update -qqy \
  && apt-get -qqy install \
    slurm-llnl\
    munge\
    heirloom-mailx\
    libslurm-dev\
    slurm-drmaa-dev\
  && rm -rf /var/lib/apt/lists/*

# Make necessary directories:
RUN mkdir /var/log/slurm/ /var/run/slurm-llnl/ \
  /var/lib/slurm/ /var/lib/slurm/slurmctld/

# Make necessary files:
RUN touch /var/log/slurm/slurmctld.log /var/run/slurm-llnl/slurmctld.pid

# Permissions:
RUN chown -R slurm /var/lib/slurm/ /var/log/slurm/ /var/run/slurm-llnl/

# Configuration file:
ADD slurm.conf /etc/slurm-llnl/slurm.conf

# Slurm setup:
ADD slurm.sh /root/

#==============================================================================
# Python / Python-DRMAA Bindings
#==============================================================================
RUN apt-get update -qqy \
  && apt-get -qqy install \
    python-dev\
    python-pip\
  && rm -rf /var/lib/apt/lists/*

RUN pip install drmaa

ENV SGE_ROOT=/var/lib/slurm
ENV SGE_CELL=default
ENV DRMAA_LIBRARY_PATH=/usr/lib/libdrmaa.so

#==============================================================================
# BLAST+
#==============================================================================
RUN apt-get update -qqy \
  && apt-get -qqy install \
    blast2\
    ncbi-blast+\
  && rm -rf /var/lib/apt/lists/*

#==============================================================================
# NEPTUNE
#==============================================================================
RUN pip install --upgrade pip

RUN pip install numpy
RUN pip install scipy
RUN pip install biopython

RUN pip install tox

#==============================================================================
# FINISH
#==============================================================================
WORKDIR /root/


