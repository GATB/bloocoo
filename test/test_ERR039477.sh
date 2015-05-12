#! /bin/bash

################################################################################
# we download some banks
################################################################################
wget "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR039/ERR039477/ERR039477.fastq.gz"

################################################################################
# we launch bloocoo
################################################################################
../build/Bloocoo -file ./ERR039477.fastq.gz  -nb-cores 1

################################################################################
# we check the result
################################################################################
md5sum ERR039477_corrected.fastq > ERR039477.check

diff ./ERR039477.md5 ./ERR039477.check

if [ $? -eq 0 ]; then
   echo "TEST OK"
else
   echo "TEST KO"
fi

################################################################################
# clean up
################################################################################
rm -f  ERR039477.fastq.gz  ERR039477.h5  ERR039477_corrected.fastq  ERR039477.check
