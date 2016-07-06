#!/bin/bash 

# Neptune directory:
DIR=`dirname $0`

. $DIR/.venv/bin/activate
neptune $@
