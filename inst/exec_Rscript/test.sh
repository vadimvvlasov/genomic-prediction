#!/bin/bash

DIR=$(dirname $0)
cd $DIR
DIR=$(pwd)

echo $DIR
ls -lh ${DIR}/input
