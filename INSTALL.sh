#!/bin/bash 

# Neptune directory:
DIR=`dirname $0`

# Create Python virtual environment:
VENV=$DIR/.venv
virtualenv $VENV
. $VENV/bin/activate

# PIP
pip install --upgrade pip

# NUMPY
pip install numpy

# SCIPY
pip install scipy

# BIOPYTHON
pip install biopython

# DRMAA
pip install drmaa

# NEPTUNE
pip install $DIR/
