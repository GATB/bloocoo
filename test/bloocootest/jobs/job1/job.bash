#!/bin/bash

#$ -S /bin/bash
#$ -M gaetan.benoit@inria.fr
#$ -m bea
#$ -cwd 

#usage: 
#	Faut que les parametre soit dans l'ordre, le -g est facultatif mais doit etre le dernier param
#	python main.py <genome_size> <reads_size> <reads_count> <kmer_size> <bloom_threshold> <error_rate>
#	-g: generate a random genome

#Creation du fichier ou seront mis les resultats
rm -rf test_result
mkdir -p test_result

python test.py
