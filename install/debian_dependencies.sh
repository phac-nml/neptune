#!/bin/bash 

# PIP
wget https://bootstrap.pypa.io/get-pip.py
python get-pip.py

# VIRTUALENV
pip install virtualenv

# BUILD ESSENTIAL, PYTHON DEV, BLAST+
apt-get --yes --force-yes install build-essential
apt-get --yes --force-yes install python-dev
apt-get --yes --force-yes install ncbi-blast+
