#!/bin/bash 

# Create Python virtual environment:
virtualenv .venv
. .venv/bin/activate

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
pip install .

