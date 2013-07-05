#!/bin/bash

#$ -S /bin/bash
#$ -M gaetan.benoit@inria.fr
#$ -m bea
#$ -cwd 

source /local/env/envR

cd jobs/job1
bash job.bash

