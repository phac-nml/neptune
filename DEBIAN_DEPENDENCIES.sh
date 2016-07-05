#!/bin/bash 

wget https://bootstrap.pypa.io/get-pip.py
python get-pip.py

pip install virtualenv

apt-get install build-essential
apt-get install python-dev
