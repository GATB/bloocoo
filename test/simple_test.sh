#!/bin/bash
#simple  test with synthetic data

echo -n "Testing isolated error, middle of read"

../build/Bloocoo -db datatest/errm.fasta -kmer-size 31 -nks 5 &> /dev/null

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

rm errm_bloocoo_corr_errs.tab errm_corrected.fasta



echo -n "Testing isolated error, left side"

../build/Bloocoo -db datatest/errleft.fasta -kmer-size 31 -nks 5 &> /dev/null


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

rm errleft_bloocoo_corr_errs.tab errleft_corrected.fasta






echo -n "Testing isolated error, right side"

../build/Bloocoo -db datatest/errright.fasta -kmer-size 31 -nks 5 &> /dev/null


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

rm errright_bloocoo_corr_errs.tab errright_corrected.fasta


echo -n "Testing two close errors"

../build/Bloocoo -db datatest/errclose.fasta -kmer-size 31 -nks 5 &> /dev/null


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
