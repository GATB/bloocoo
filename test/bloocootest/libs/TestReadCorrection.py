
import sys, os
import cPickle, commands, copy, time
from os.path import exists, basename, splitext, join


class TestReadCorrection:

	TIMER = {
		"bloocoo": 0,
		"test": 0,
	}


	@staticmethod
	def main(params):
		TestReadCorrection.remove_file("reads.fasta.bin")

		#args = sys.argv
		genome_size = str(params["genome_size"])
		reads_size = str(params["reads_size"])
		reads_count = str(params["reads_count"])
		kmer_size = str(params["kmer_size"])
		coverage_threshold = str(params["coverage_threshold"])
		error_rate = str(params["error_rate"])
		nb_iter = str(params["nb_iter"])
		nb_kmer_checked = str(params["nb_kmer_checked"])
		genome_filename = str(params["genome_filename"])
		
		
		#Mutaread
		if genome_filename == "":
			#Generation d'un genome aleatoirement, le genome est contenu dans le fichier alea.seq
			os.system("../../gener_alea " + genome_size + " 1")
			#Decoupage du genome contenu dans le fichier en multiple reads
			#Params: file_in, file_out, reads_count, reads_length, error_rate, error_rate, error_rate
			os.system("../../mutareads alea.seq  reads " + reads_count + " " + reads_size + " " + error_rate + " 0 0 -errfile")
		else:
			os.system("../../mutareads " + join("../../genomes", genome_filename) + " reads " + reads_count + " " + reads_size + " " + error_rate + " 0 0 -errfile")
		
		#Bloocoo correction
		t = time.time()
		#command = ./Bloocoo  datatest/errm.fasta   31 -nks 3 -nb-iter 2 -nkmer-checked 0
		os.system("../../Bloocoo -db reads.fasta -kmer-size " + kmer_size + " -nks " + coverage_threshold + " -nb-iter " + nb_iter + " -nkmer-checked " + nb_kmer_checked)
		TestReadCorrection.TIMER["bloocoo"] = time.time() - t
		TestReadCorrection.remove_file("reads.fasta.bin")
		TestReadCorrection.remove_file("foo_partition0.out")
		TestReadCorrection.remove_file("solids.bin")

		err_tab_filename = "reads_errs.tab"
		corrected_tab_filename = "reads_bloocoo_corr_errs.tab"

		TestReadCorrection.print_params(params)
		TestReadCorrection.execute(params, err_tab_filename, corrected_tab_filename)

		TestReadCorrection.remove_file("alea.seq")
		TestReadCorrection.remove_file("reads_corrected.fasta")
		TestReadCorrection.remove_file("sorted_"+err_tab_filename)
		TestReadCorrection.remove_file("sorted_"+corrected_tab_filename)
		TestReadCorrection.remove_file("test_diff_result.temp")
		#TestReadCorrection.remove_file("reads_bloocoo_corr_errs.tab")
		#TestReadCorrection.remove_file("reads_errs.tab")
		TestReadCorrection.remove_file("comm.temp")
		TestReadCorrection.remove_file("tmp.binary")
		TestReadCorrection.remove_file("tmp.solid")
		#TestReadCorrection.remove_file("reads.fasta")
		
	@staticmethod
	def remove_file(filename):
		try:
			os.remove(filename)
		except:
			pass

	@staticmethod
	def execute(params, err_tab_filename, corrected_tab_filename):
		err_tab_file = open(err_tab_filename, "r")
		corrected_tab_file = open(corrected_tab_filename, "r")
		#---
		total_error_count = int(commands.getstatusoutput("wc -l " + err_tab_filename)[1].split(" ")[0])
		#---
		t = time.time()
		os.system("sort " + err_tab_filename + "> sorted_"+err_tab_filename)
		os.system("sort " + corrected_tab_filename + "> sorted_"+corrected_tab_filename)
		print "Time for sorting both tabs:", time.time() - t
		#---
		t = time.time()
		command = "comm -1 " + "sorted_"+err_tab_filename + " " + "sorted_"+corrected_tab_filename + " > comm.temp"
		os.system(command)
		predit = int(commands.getstatusoutput("wc -l comm.temp")[1].split(" ")[0])
		print "Time to compare (predit):", time.time() - t
		#---
		command = "comm -12 " + "sorted_"+err_tab_filename + " " + "sorted_"+corrected_tab_filename + " > comm.temp"
		os.system(command)
		error_corrected_count = int(commands.getstatusoutput("wc -l comm.temp")[1].split(" ")[0])
		#t = time.time()
		#command = "comm -13 " + "sorted_"+err_tab_filename + " " + "sorted_"+corrected_tab_filename + " > comm.temp"
		#os.system(command)
		#wrong_error_corrected_count = int(commands.getstatusoutput("wc -l comm.temp")[1].split(" ")[0])
		print "Time to compare (vp):", time.time() - t
		#---
		TestReadCorrection.print_result(params, total_error_count, predit, error_corrected_count)

	@staticmethod
	def print_params(params):
		print "----------------------------------------"
		print "- Test Start"
		for param_name, param_value in params.items():
			print "-\t" + param_name + ": " + str(param_value)
		print "----------------------------------------"

	@staticmethod
	def print_result(params, total_error_count, predit, error_corrected_count):
		vrai = total_error_count
		vp = error_corrected_count
		fp = predit - error_corrected_count
		print "----------------------------------------"
		print "- Test Result"
		print "-\tVrai (Nombre total d\'erreur):", vrai
		print "-\tVP (Nombre total d\'erreur corrigees):", vp
		if total_error_count == 0:
			recall = 100
		else:
			recall = round((float(vp)/vrai)*100,1)
		print "-\tRecall (Taux Erreur corrigee): " + str(recall) + "%"
		print "-\tPredit: " + str(predit)
		print "-\tFP (Erreur ajoutees):", fp
		if predit == 0:
			precision = 100
		else:
			precision = round((float(vp)/predit)*100,1)
		
		print "-\tPrecision: " + str(precision) + "%"
		print "----------------------------------------"
		print "- Temps Execution"
		print "-\tBloocoo: " + str(TestReadCorrection.TIMER["bloocoo"])
		print "----------------------------------------"
		#---
		if not exists("test_result"):
			os.mkdir("test_result")
		coverage = params["reads_size"] * params["reads_count"] / float(params["genome_size"])
		real_coverage = int(coverage * (params["reads_size"] - params["kmer_size"] + 1) / params["reads_size"])
		
		filename = "test_result/"
		filename += params["result_filename_prefix"] + "_"
		filename += splitext(basename(params["genome_filename"]))[0]
		filename += "_c"+str(int(real_coverage)) + "_k"+str(params["kmer_size"])
		filename += ".tab"
		if not exists(filename):
			file = open(filename, "w")
			file.write("err\tcov_t\trec\tpre\n")
		else:
			file = open(filename, "a")
		#---
		file.write(str(params["error_rate"]) + "\t" + str(params["coverage_threshold"]) + "\t" + str(recall) + "\t" + str(precision) + "\n")


