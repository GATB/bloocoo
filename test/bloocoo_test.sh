#!/bin/bash
#simple  test with synthetic data

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

echo "Testing two close errors..."

$bindir/Bloocoo -file datatest/errclose.fasta -kmer-size 31 -abundance-min 5 -err-tab &> /dev/null

diff errclose_bloocoo_corr_errs.tab ./datatest/true_res_test4 > /dev/null

var=$?

# SOME CLEANUP
rm  *_corrected.fasta  *.tab  *.h5

if [ $var -eq 0 ]
then
echo  PASSED
#    exit 0
else
echo  FAILED
exit 1
fi

exit 0
