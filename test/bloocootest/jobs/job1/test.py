
#Add the root tests folder to the pythonpath
import sys; sys.path.append("../..")
#import all tests methods
from libs.TestExec import *


#variables parametrables (valeur par defaut)
#
#Correction software:
#  - 0: Bloocoo
#  - 1: Musket
#  - 2: Racer

p = {
	"genome_filename": "",
	"genome_size": 5000,
	"reads_size": 100,
	"cover": 30,
	"kmer_size": 30,
	"bloom_threshold": 4,
	"error_rate": 0.01,
	"result_filename_prefix": "test",
	"result_filename_suffix": "",
	"R_script_filename": "rec_pre_vs_err.R",
	"nb_kmer_checked": 4,
	"min_nb_kmer_checked": 4,
	"regenerate_reads": True,
	"correction_software": 0, 
	
}


#["ecoli.fasta", "humch22.fasta", ""]

for kmer_size in [30]:
	for cover in [30]:
		for genome_filename in [""]:
		
			p["genome_filename"] = genome_filename
			p["kmer_size"] = kmer_size
			p["cover"] = cover

			#Always call setup() before executing tests to prepare the new genome
			setup_test(p)
			
			#-----------------------------------------------------
			# recall/precision en fonction du taux d'erreur
			#-----------------------------------------------------
			p["result_filename_prefix"] = "test_error"
			p["result_filename_suffix"] = "_c"+str(p["cover"]) + "_k"+str(p["kmer_size"])
			p["bloom_threshold"] = 4
			p["R_script_filename"] = "rec_pre_vs_err.R"
			#error_rates = [0.001, 0.01, 0.02, 0.03, 0.05]
			error_rates = [0.01]

			for error_rate in error_rates:
				p["error_rate"] = error_rate
				
				p["nb_kmer_checked"] = 4
				p["min_nb_kmer_checked"] = 3
				execute_test(p)
				p["nb_kmer_checked"] = 4
				p["min_nb_kmer_checked"] = 4
				execute_test(p)
				p["nb_kmer_checked"] = 6
				p["min_nb_kmer_checked"] = 5
				execute_test(p)
				p["nb_kmer_checked"] = 6
				p["min_nb_kmer_checked"] = 6
				execute_test(p)
				p["nb_kmer_checked"] = 8
				p["min_nb_kmer_checked"] = 6
				execute_test(p)
				p["nb_kmer_checked"] = 8
				p["min_nb_kmer_checked"] = 7
				execute_test(p)
				p["nb_kmer_checked"] = 8
				p["min_nb_kmer_checked"] = 8
				execute_test(p)
				p["nb_kmer_checked"] = 10
				p["min_nb_kmer_checked"] = 9
				execute_test(p)
				p["nb_kmer_checked"] = 10
				p["min_nb_kmer_checked"] = 10
				execute_test(p)
				p["nb_kmer_checked"] = 12
				p["min_nb_kmer_checked"] = 8
				execute_test(p)
				p["nb_kmer_checked"] = 12
				p["min_nb_kmer_checked"] = 10
				execute_test(p)
				p["nb_kmer_checked"] = 12
				p["min_nb_kmer_checked"] = 12
				execute_test(p)
				
			#execute_graph(p)

			#-----------------------------------------------------
			# recall/precision avec taux d'erreur fixe et seuil variable
			#-----------------------------------------------------
			"""
			p["result_filename_prefix"] = "bloocootest_threshold"
			p["error_rate"] = 0.01
			p["R_script_filename"] = "rec_pre_vs_threshold.R"
			cov_ts = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 35, 40]

			for bloom_threshold in cov_ts:
				p["bloom_threshold"] = coverage_threshold
				execute_test(p)

			execute_graph(p)
			"""


