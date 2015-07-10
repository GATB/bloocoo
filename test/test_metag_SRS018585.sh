#!/bin/bash

#this dataset comes from http://2014-5-metagenomics-workshop.readthedocs.org/en/latest/assembly/qtrim.html
# human microbiome project

set -e # exit on error

################################################################################
# we download some banks
################################################################################
if [[ ! -f SRS018585/SRS018585.denovo_duplicates_marked.trimmed.1.fastq ]]
then
    wget http://downloads.hmpdacc.org/data/Illumina/anterior_nares/SRS018585.tar.bz2
    tar -xjf SRS018585.tar.bz2
fi

################################################################################
# we launch bloocoo
################################################################################
#../build/Bloocoo -file SRS018585/SRS018585.denovo_duplicates_marked.trimmed.1.fastq
../build/metaBloocoo count -file SRS018585/SRS018585.denovo_duplicates_marked.trimmed.1.fastq
../build/metaBloocoo correct -file SRS018585/SRS018585.denovo_duplicates_marked.trimmed.1.fastq -t1 2 -t3 5


if [ $? -eq 0 ]; then
   echo "TEST OK"
else
   echo "TEST KO"
fi

################################################################################
# clean up
################################################################################
rm -f  SRS018585_corrected* SRS018585_count*
