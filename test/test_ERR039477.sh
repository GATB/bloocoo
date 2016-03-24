#! /bin/bash

# look for Bloocoo binary. In devel mode, it's in ../build/bin directory.
# In production mode, it's in ../bin directory.
if [ -f "../bin/Bloocoo" ]
then
 bindir="../bin"
elif [ -f "../build/bin/Bloocoo" ]
then
 bindir="../build/bin"
else
 echo "could not find a compiled Bloocoo binary"
 exit 1
fi

################################################################################
# we download a sample bank from EBI
################################################################################
# if wget is not installed, you may use "curl -O ..."
DATA_SAMPLE="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR039/ERR039477/ERR039477.fastq.gz"
WGET_PATH=`which wget`
echo ">>> Retrieving data sample: ${DATA_SAMPLE}"
if [ ! -z "$WGET_PATH" ] ; then
  echo "    using '$WGET_PATH'..."
  if [ $SILENT_MODE=="true"  ] ; then
    wget --quiet ${DATA_SAMPLE}
  else
    wget ${DATA_SAMPLE}
  fi
else
   CURL_PATH=`which curl`
  if [ ! -z "$CURL_PATH" ] ; then
    echo "    using '$CURL_PATH'..."
    if [ $SILENT_MODE=="true"  ] ; then
      curl --silent -O ${DATA_SAMPLE}
    else
      curl -O ${DATA_SAMPLE}
    fi
  else
    echo "    /!\ error: unable to find 'wget' or 'curl'"
    exit 1
  fi
fi


################################################################################
# we launch bloocoo
################################################################################
$bindir/Bloocoo -file ./ERR039477.fastq.gz  -nb-cores 1

################################################################################
# we check the result
################################################################################
md5sum ERR039477.fastq_corrected.fastq > ERR039477.check

diff ./ERR039477.md5 ./ERR039477.check

if [ $? -eq 0 ]; then
   echo "TEST OK"
else
   echo "TEST KO"
fi

################################################################################
# clean up
################################################################################
rm -f  ERR039477.fastq.gz  ERR039477.fastq.h5  ERR039477.fastq_corrected.fastq  ERR039477.check
