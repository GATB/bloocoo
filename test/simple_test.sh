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

echo -n "Testing isolated error, middle of read..."

$bindir/Bloocoo -file datatest/errm.fasta -kmer-size 31 -abundance-min 5  -err-tab &> /dev/null

diff errm_bloocoo_corr_errs.tab ./datatest/true_res_test1 > /dev/null

var=$?
if [ $var -eq 0 ]
then
    echo  PASSED
#    exit 0
else
    echo  FAILED
    exit 1
fi


echo -n "Testing isolated error, left side..."

$bindir/Bloocoo -file datatest/errleft.fasta -kmer-size 31 -abundance-min 5 -err-tab  &> /dev/null


diff errleft_bloocoo_corr_errs.tab ./datatest/true_res_test2 > /dev/null

var=$?
if [ $var -eq 0 ]
then
    echo  PASSED
#    exit 0
else
    echo  FAILED
    exit 1
fi



echo -n "Testing isolated error, right side..."

$bindir/Bloocoo -file datatest/errright.fasta -kmer-size 31 -abundance-min 5 -err-tab  &> /dev/null


diff errright_bloocoo_corr_errs.tab ./datatest/true_res_test3 > /dev/null

var=$?
if [ $var -eq 0 ]
then
echo  PASSED
#    exit 0
else
echo  FAILED
exit 1
fi


echo -n "Testing two close errors..."

$bindir/Bloocoo -file datatest/errclose.fasta -kmer-size 31 -abundance-min 5 -err-tab  &> /dev/null


diff errclose_bloocoo_corr_errs.tab ./datatest/true_res_test4 > /dev/null

var=$?
if [ $var -eq 0 ]
then
echo  PASSED
#    exit 0
else
echo  FAILED
exit 1
fi

# SOME CLEANUP
rm  *_corrected.fasta  *.tab  *.h5

exit 0
