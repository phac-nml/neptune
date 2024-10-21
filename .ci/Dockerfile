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
# Python
#==============================================================================
RUN apt-get update -qqy \
  && apt-get -qqy install \
    python-dev\
    python-pip\
  && rm -rf /var/lib/apt/lists/*

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


