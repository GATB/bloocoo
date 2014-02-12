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



echo -n "Testing 2nt-deletion error, middle of read"

../build/Bloocoo -db datatest/err2del.fasta -kmer-size 31 -nks 5   -ion  &> /dev/null

diff err2del_corrected.fasta ./datatest/true_del2 > /dev/null


var=$?
if [ $var -eq 0 ] 
then
    echo  PASSED
#    exit 0
else
    echo  FAILED
    exit 1
fi

rm err2del_corrected.fasta


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



echo -n "Testing 2-nt deletion error, end of read"

../build/Bloocoo -db datatest/err2delfin.fasta -kmer-size 31 -nks 5   -ion  &> /dev/null

diff err2delfin_corrected.fasta ./datatest/true2errdelfin > /dev/null


var=$?
if [ $var -eq 0 ] 
then
    echo  PASSED
#    exit 0
else
    echo  FAILED
    exit 1
fi

rm err2delfin_corrected.fasta



echo -n "Testing 2-nt insertion error, end of read"

../build/Bloocoo -db datatest/err2insfin.fasta -kmer-size 31 -nks 5   -ion  &> /dev/null

diff err2insfin_corrected.fasta ./datatest/true_err2insfin > /dev/null


var=$?
if [ $var -eq 0 ] 
then
    echo  PASSED
#    exit 0
else
    echo  FAILED
    exit 1
fi

rm err2insfin_corrected.fasta





exit 0
