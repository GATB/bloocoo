
#importe toute les methodes de tests
import sys; sys.path.insert(0, "../.."); from libs.TestExec import *


#variables parametrables (valeur par defaut)
p = {
	"genome_filename": "",
	"genome_size": 5000,
	"reads_size": 100,
	"cover": 30,
	"kmer_size": 30,
	"coverage_threshold": 4,
	"error_rate": 0.01,
	"result_filename_prefix": "bloocoo_test",
	"R_script_filename": "rec_pre_vs_err.R",
	"nb_iter": 2,
	"nb_kmer_checked": 4,
}


#["ecoli.fasta", "humch22.fasta", ""]

for kmer_size in [30]:
	for cover in [30]:
		for genome_filename in [""]:
		
			p["genome_filename"] = genome_filename
			p["kmer_size"] = kmer_size
			p["cover"] = cover

			#Always call setup() before executing tests to prepare the new genome
			setup(p)
			#-----------------------------------------------------
			# recall/precision en fonction du taux d'erreur
			#-----------------------------------------------------
			p["result_filename_prefix"] = "bloocootest_error"
			p["coverage_threshold"] = 4
			p["R_script_filename"] = "rec_pre_vs_err.R"
			#error_rates = [0.001, 0.01, 0.02, 0.03, 0.05]
			error_rates = [0.01]

			for error_rate in error_rates:
				p["error_rate"] = error_rate
				execute_test(p)

			execute_graph(p)

			#-----------------------------------------------------
			# recall/precision avec taux d'erreur fixe et seuil variable
			#-----------------------------------------------------
			"""
			p["result_filename_prefix"] = "bloocootest_threshold"
			p["error_rate"] = 0.01
			p["R_script_filename"] = "rec_pre_vs_threshold.R"
			cov_ts = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 35, 40]

			for coverage_threshold in cov_ts:
				p["coverage_threshold"] = coverage_threshold
				execute_test(p)

			execute_graph(p)
			"""


