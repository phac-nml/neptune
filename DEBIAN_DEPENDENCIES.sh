#!/bin/bash 

# PIP
wget https://bootstrap.pypa.io/get-pip.py
python get-pip.py

# VIRTUALENV
pip install virtualenv

# BUILD ESSENTIAL, PYTHON DEV
apt-get --yes install build-essential
apt-get --yes install python-dev
