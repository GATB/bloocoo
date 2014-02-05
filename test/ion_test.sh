#!/bin/bash
#simple  test with synthetic data

echo -n "Testing insertion error, middle of read"

../build/Bloocoo -db datatest/errins.fasta -kmer-size 31 -nks 5   -ion  &> /dev/null

diff errins_corrected.fasta ./datatest/true_errins1 > /dev/null


var=$?
if [ $var -eq 0 ] 
then
    echo  PASSED
#    exit 0
else
    echo  FAILED
    exit 1
fi

rm errins_corrected.fasta


echo -n "Testing deletion error, middle of read"

../build/Bloocoo -db datatest/errdel.fasta -kmer-size 31 -nks 5   -ion  &> /dev/null

diff errdel_corrected.fasta ./datatest/true_errdel1 > /dev/null


var=$?
if [ $var -eq 0 ] 
then
    echo  PASSED
#    exit 0
else
    echo  FAILED
    exit 1
fi

rm errdel_corrected.fasta



echo -n "Testing 2nt-insertion error, middle of read"

../build/Bloocoo -db datatest/err2ins.fasta -kmer-size 31 -nks 5   -ion  &> /dev/null

diff err2ins_corrected.fasta ./datatest/true_errins2 > /dev/null


var=$?
if [ $var -eq 0 ] 
then
    echo  PASSED
#    exit 0
else
    echo  FAILED
    exit 1
fi

rm err2ins_corrected.fasta



echo -n "Testing insertion error, end of read"

../build/Bloocoo -db datatest/errinsfin.fasta -kmer-size 31 -nks 5   -ion  &> /dev/null

diff errinsfin_corrected.fasta ./datatest/true_errinsfin > /dev/null


var=$?
if [ $var -eq 0 ] 
then
    echo  PASSED
#    exit 0
else
    echo  FAILED
    exit 1
fi

rm errinsfin_corrected.fasta



echo -n "Testing deletion error, end of read"

../build/Bloocoo -db datatest/errdelfin.fasta -kmer-size 31 -nks 5   -ion  &> /dev/null

diff errdelfin_corrected.fasta ./datatest/true_errdelfin > /dev/null


var=$?
if [ $var -eq 0 ] 
then
    echo  PASSED
#    exit 0
else
    echo  FAILED
    exit 1
fi

rm errdelfin_corrected.fasta



exit 0
