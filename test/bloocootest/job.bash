#!/bin/bash

#$ -S /bin/bash
#$ -M gaetan.benoit@inria.fr
#$ -m bea
#$ -cwd 

#for quake
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:~/workspace/gatb-tools/gatb-tools/tools/bloocoo/test/bloocootest/quake/jellyfish/lib
export PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:~/workspace/gatb-tools/gatb-tools/tools/bloocoo/test/bloocootest/quake/jellyfish/lib/pkgconfig

source /local/env/envR

cd jobs/job1
bash job.bash

