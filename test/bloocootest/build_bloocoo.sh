#!/bin/sh

pwd=$PWD

cd

cd workspace/gatb-tools/gatb-tools/tools/bloocoo/build/
make
cp Bloocoo ~/workspace/scripts/bloocootest/

cd ${pwd}
