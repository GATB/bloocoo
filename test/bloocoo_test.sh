#!/bin/bash
#simple  test with synthetic data

echo -n "Testing two close errors"

../build/Bloocoo -db datatest/errclose.fasta -kmer-size 31 -nks 5


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

rm errclose_bloocoo_corr_errs.tab errclose_corrected.fasta




exit 0
