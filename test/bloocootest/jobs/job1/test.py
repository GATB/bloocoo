
#Add the root tests folder to the pythonpath
import sys; sys.path.append("../..")
#import all tests methods
from libs.TestExec import *


#variables parametrables (valeur par defaut)

p = {
	"genome_filename": "ecoli.fasta",
	"genome_size": 10000,
	"reads_size": 100,
	"cover": 30,
	"kmer_size": 24,
	"bloom_threshold": 4,
	"error_rate": 0.01,
	"result_filename_prefix": "test",
	"result_filename_suffix": "",
	"R_script_filename": "rec_pre_vs_err.R",
	"nb_kmer_checked": 8,
	"min_nb_kmer_checked": 8,
	"regenerate_reads": True,
	"correction_software": 1, 
	"read_format":"fasta"
}



#Always call setup() before executing tests to prepare the new genome
setup_test(p)

p["result_filename_prefix"] = "test"

execute_test(p)


