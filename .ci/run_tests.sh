#!/bin/bash -e

printf "setenv = \n    DRMAA_LIBRARY_PATH = /usr/lib/gridengine-drmaa/lib/libdrmaa.so.1.0\n" >> tox.ini
printf "    SGE_ROOT = /opt/gridengine\n" >> tox.ini
printf "    SGE_CELL=default\n" >> tox.ini
printf "    SGE_QMASTER_PORT=536\n" >> tox.ini
printf "    SGE_EXECD_PORT=537\n" >> tox.ini
printf "    SGE_ARCH=linux-x64\n" >> tox.ini
tox -r
