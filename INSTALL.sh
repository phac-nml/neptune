#!/bin/bash 

# Neptune directory:
DIR=`dirname $0`

# Install directory:
PREFIX=$DIR

# Any location provided?
if [ "$1" != "" ]; then
    PREFIX=$1
fi

BIN=$PREFIX/bin
LIB=$PREFIX/lib

# Create directories if needed:
if [ ! -d "$PREFIX" ]; then
    mkdir $PREFIX
fi

if [ ! -d "$BIN" ]; then
    mkdir $BIN
fi

if [ ! -d "$LIB" ]; then
    mkdir $LIB
fi

# Create Python virtual environment:
VENV=$LIB/neptune
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

# BIN
cp $DIR/install/neptune $BIN/neptune

