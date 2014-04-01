#!/bin/bash
#simple  test with synthetic data

echo -n "Testing two close errors"

../build/Bloocoo -file datatest/errclose.fasta -kmer-size 31 -nks 5 -err-tab &> /dev/null


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
