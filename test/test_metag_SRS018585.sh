#!/bin/bash

#this dataset comes from http://2014-5-metagenomics-workshop.readthedocs.org/en/latest/assembly/qtrim.html
# human microbiome project

#set -e # exit on error

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
# we download some banks
################################################################################
#if [[ ! -f SRS018585/SRS018585.denovo_duplicates_marked.trimmed.1.fastq ]]
#then
#    wget http://downloads.hmpdacc.org/data/Illumina/anterior_nares/SRS018585.tar.bz2
#    tar -xjf SRS018585.tar.bz2
#fi

################################################################################
# we download a sample bank from EBI
################################################################################
# if wget is not installed, you may use "curl -O ..."
DATA_SAMPLE="http://downloads.hmpdacc.org/data/Illumina/anterior_nares/SRS018585.tar.bz2"
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

echo ">>> Unarchiving data sample: SRS018585.tar.bz2"

tar -xjf SRS018585.tar.bz2

################################################################################
# we launch bloocoo
################################################################################
#$bindir/Bloocoo -file SRS018585/SRS018585.denovo_duplicates_marked.trimmed.1.fastq
$bindir/metaBloocoo count -file SRS018585/SRS018585.denovo_duplicates_marked.trimmed.1.fastq #-out SRS018585 #-alt_ascii_graph
$bindir/metaBloocoo correct -file SRS018585/SRS018585.denovo_duplicates_marked.trimmed.1.fastq -t1 2 -t2 5 -t3 10 #-out SRS018585

################################################################################
# clean up
################################################################################
rm -rf  SRS018585*

if [ $? -eq 0 ]; then
   echo "TEST OK"
else
   echo "TEST KO"
   exit 1
fi
