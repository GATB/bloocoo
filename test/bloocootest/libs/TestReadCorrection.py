 
import sys, os
import cPickle, commands, copy, time
from os.path import exists, basename, splitext, join
from itertools import izip

#-------------------------------------------------------------
# ** TestReadCorrection
#-------------------------------------------------------------
class TestReadCorrection:

	#-------------------------------------------------------------
	# * class variables
	#-------------------------------------------------------------


	#-------------------------------------------------------------
	# * main
	#-------------------------------------------------------------
	@staticmethod
	def main(params):
		
		#Print all test params before starting the tests
		TestReadCorrection.print_params(params)
			
		#Correction will fail if this file exists
		TestReadCorrection.remove_file("reads.fasta.bin")

		#Convert params to variable
		genome_size = str(params["genome_size"])
		reads_size = str(params["reads_size"])
		reads_count = str(params["reads_count"])
		kmer_size = str(params["kmer_size"])
		bloom_threshold = str(params["bloom_threshold"])
		error_rate = str(params["error_rate"])
		nb_kmer_checked = str(params["nb_kmer_checked"])
		min_nb_kmer_checked = str(params["min_nb_kmer_checked"])
		regenerate_reads = bool(params["regenerate_reads"])
		genome_filename = str(params["genome_filename"])
		correction_software = int(params["correction_software"])
		
		#Mutaread
		if regenerate_reads:
			if genome_filename == "":
				#Generation d'un genome aleatoirement, le genome est contenu dans le fichier alea.seq
				os.system("../../gener_alea " + genome_size + " 1")
				#Decoupage du genome contenu dans le fichier en multiple reads
				#Params: file_in, file_out, reads_count, reads_length, error_rate, error_rate, error_rate
				os.system("../../mutareads alea.seq  reads " + reads_count + " " + reads_size + " " + error_rate + " 0 0 -errfile")
			else:
				os.system("../../mutareads " + join("../../genomes", genome_filename) + " reads " + reads_count + " " + reads_size + " " + error_rate + " 0 0 -errfile ")
		
		#Reads correction
		t = time.time()
		
		#Bloocoo correction
		if correction_software == 0:
			os.system("../../Bloocoo -verbose -db reads.fasta -kmer-size " + kmer_size + " -nks " + bloom_threshold + " -nbmin-valid " + min_nb_kmer_checked + " -nkmer-checked " + nb_kmer_checked)
			correction_duration = time.time() - t
			#Remove bloocoo temp files
			TestReadCorrection.remove_file("reads.fasta.bin")
			TestReadCorrection.remove_file("foo_partition0.out")
			TestReadCorrection.remove_file("solids.bin")
			#Start tests specific for bloocoo
			TestReadCorrection.print_params(params)
			true, tp, fp, FP0, FP1, FP2 ,FP3, FP4, TC0, TC1, TC2 ,TC3, TC4 = BloocooTest.execute(params, "reads_errs.tab", "reads_bloocoo_corr_errs.tab", "reads_bloocoo_corr_errs_full.tab")
		#Musket correction
		elif correction_software == 1:
			#os.system("../../musket reads.fasta -o reads_corrected.fasta")
			os.system("../../musket -k " + kmer_size + " " + genome_size + " -maxiter " + nb_kmer_checked + " -p 12 -inorder reads.fasta -o reads_corrected.fasta")
			correction_duration = time.time() - t
		#Racer correction
		elif correction_software == 2:
			os.system("../../racer reads.fasta reads_corrected.fasta " + str(genome_size))
			correction_duration = time.time() - t
		#
		elif correction_software == 3: # bloocoo sans errtab
			os.system("../../Bloocoo -verbose -db reads.fasta -kmer-size " + kmer_size + " -nks " + bloom_threshold + " -nbmin-valid " + min_nb_kmer_checked + " -nkmer-checked " + nb_kmer_checked)
			correction_duration = time.time() - t
		elif correction_software == 4:
			pass
			
		#Use generic tests if bloocoo is not used
		if correction_software != 0:
			TestReadCorrection.print_params(params)
			true, fp, tp = GenericTest.execute(params, "reads.fasta", "reads_corrected.fasta", "reads_errs.tab")
			
		#Recall
		predicted = tp + fp
		fn = true - tp
		if true == 0:
			recall = 100
		else:
			recall = float('%.2f' % ((float(tp)/true)*100))
			
		#Precision
		if predicted == 0:
			precision = 0
		else:
			precision = float('%.2f' % ((float(tp)/predicted)*100))
			
		#Fscore
		if precision+recall == 0:
			Fscore = 0
		else:
			Fscore = (2*precision*recall) / (precision+recall)
		Fscore = '%.2f' % Fscore
		
		#Gain
		if true == 0:
			gain = 0
		else:
			gain = ((tp-fp) / float(true))*100
		gain = '%.2f' % gain
		
		if correction_software == 0: #Bloocoo
			BloocooTest.print_stats(recall, precision, FP0, FP1, FP2 ,FP3, FP4, TC0, TC1, TC2 ,TC3, TC4)
				
		#Print result
		TestReadCorrection.print_result(params, true, tp, fp, recall, precision, Fscore, gain, correction_duration)

		#Write result to tabs
		if not exists("test_result"):
			os.mkdir("test_result")
		#coverage = params["reads_size"] * params["reads_count"] / float(params["genome_size"])
		#real_coverage = int(coverage * (params["reads_size"] - params["kmer_size"] + 1) / params["reads_size"])
		
		#tab_filename = "test_result/"
		tab_filename = params["result_filename_prefix"] + "_"
		tab_filename += splitext(basename(params["genome_filename"]))[0]
		#filename += "_c"+str(int(real_coverage)) + "_k"+str(params["kmer_size"])
		tab_filename += params["result_filename_suffix"]
		tab_filename += ".tab"
		
		#Create tabs file and write the column names
		if not exists(join("test_result", tab_filename)):
			tab_file = open(join("test_result", tab_filename), "w")
			formatted_tab_file = open(join("test_result", "formatted_" + tab_filename), "w")
			tab_column_names = ["err", "read_cover", "bloom_t", "kmer_check", "kmer_min", "recall", "precision", "Fscore", "gain", "time"]
			for column_name in tab_column_names:
				tab_file.write(column_name + "\t")
				formatted_tab_file.write(("{:<16}").format(column_name))
			tab_file.write("\n")
			formatted_tab_file.write("\n")
		#Tab file exists so we open them at the end of the file
		else:
			tab_file = open(join("test_result", tab_filename), "a")
			formatted_tab_file = open(join("test_result", "formatted_" + tab_filename), "a")
		
		#Write the values corresponding to each column
		tab_column_values = [error_rate, params["cover"], bloom_threshold, nb_kmer_checked, min_nb_kmer_checked, recall, precision, Fscore, gain, correction_duration]
		
		for value in tab_column_values:
			tab_file.write(str(value) + "\t")
			formatted_tab_file.write(("{:<16}").format(str(value)))
		tab_file.write("\n")
		formatted_tab_file.write("\n")
		
		#Remove temp files
		TestReadCorrection.remove_file("alea.seq")
		#TestReadCorrection.remove_file("reads_corrected.fasta")
		TestReadCorrection.remove_file("sorted_reads_errs.tab")
		TestReadCorrection.remove_file("sorted_reads_bloocoo_corr_errs.tab")
		TestReadCorrection.remove_file("test_diff_result.temp")
		#TestReadCorrection.remove_file("reads_bloocoo_corr_errs.tab")
		#TestReadCorrection.remove_file("reads_errs.tab")
		#TestReadCorrection.remove_file("comm.temp")
		TestReadCorrection.remove_file("tmp.binary")
		TestReadCorrection.remove_file("tmp.solid")
		#TestReadCorrection.remove_file("reads.fasta")
		TestReadCorrection.remove_file("reads_bloocoo_corr_errs_full.tab")
		TestReadCorrection.remove_file("sorted_reads_bloocoo_corr_errs_full.tab")
		TestReadCorrection.remove_file("FN")
		TestReadCorrection.remove_file("FP")
		TestReadCorrection.remove_file("TP")
		TestReadCorrection.remove_file("FPfull")
		TestReadCorrection.remove_file("tmp.solid")
		

		
	#-------------------------------------------------------------
	# * remove_file
	#-------------------------------------------------------------
	@staticmethod
	def remove_file(filename):
		try:
			os.remove(filename)
		except:
			pass
			
	#-------------------------------------------------------------
	# * print_params
	#-------------------------------------------------------------
	@staticmethod
	def print_params(params):
		print "----------------------------------------"
		print "- Params:"
		for param_name, param_value in params.items():
			print "-\t" + param_name + ": " + str(param_value)
		print "----------------------------------------"

	#-------------------------------------------------------------
	# * print_params
	#-------------------------------------------------------------
	@staticmethod
	def print_result(params, true, tp, fp, recall, precision, Fscore, gain, correction_duration):
		
		print "----------------------------------------"
		print "- Test Result"
		print "\ttrue (" + str(true) + "), tp (" + str(tp) + "), fp (" + str(fp) + ")"
		print "-\tRecall: " + str(recall) + "%"	
		print "-\tPrecision: " + str(precision) + "%"
		print "-\tFscore: " + str(Fscore) + "%"
		print "-\tGain: " + str(gain) + "%"
		print "----------------------------------------"
		print "- Temps execution: " + str(correction_duration)
		print "----------------------------------------"

		
#-------------------------------------------------------------
# ** GenericCorrection
#-------------------------------------------------------------
class GenericTest:
	
	#-------------------------------------------------------------
	# ** ReadError
	#-------------------------------------------------------------
	class ReadError:

		def __init__(self, pos, noerr_letter, err_letter):
			self.pos = pos
			self.noerr_letter = noerr_letter
			self.err_letter = err_letter

		def __str__(self):
			return "<ReadError> [pos " + str(self.pos) + "] [" + self.noerr_letter + "->" + self.err_letter + "]"


	#-------------------------------------------------------------
	# * execute
	#-------------------------------------------------------------
	@staticmethod
	def execute(params, reads_filename, corrected_reads_filename, err_tab_filename):
		reads_file = open(reads_filename, "r")
		#reads_file.seek(1)
		corrected_reads_file = open(corrected_reads_filename, "r")
		#corrected_reads_file.seek(1)
		err_tab_file = open(err_tab_filename, "r")
		#---
		#n = 0
		read_index = 0
		true = 0
		tp = 0
		fp = 0
		
		#---
		while(True):
			#print "------------------", read_index
			err_read = GenericTest.next_read(reads_file)
			corrected_read = GenericTest.next_read(corrected_reads_file)
				
			
			if err_read == "":
				break
			#print corrected_read
			#if n == 0:
			#	n += 1
			#	continue
			#else:
			#	n = 0
			#err_read = err_line.strip()
			#corrected_read = corrected_line.strip()
			#---
			#print corrected_read
			errors = GenericTest.search_errors_at_read_index(err_tab_file, read_index)
			true += len(errors)
			tp += GenericTest.check_true_correction(corrected_read, errors)
			fp += GenericTest.check_false_correction(err_read, corrected_read, errors)
			read_index += 1
			#if read_index == 3: break
			
		return true, fp, tp
		
	@staticmethod	
	def next_read(reads_file):
		#print "-----------------------------------------"
		read = ""
		
		for line in reads_file:
			#print "lilneLALA:     ", line.strip()
			if line[0] == ">" and len(read) != 0:
				break
			elif line[0] != ">" :
				read += line.strip()
			
		#print "-----------------------------------------"	
		return read
			
	#-------------------------------------------------------------
	# * search_errors_at_read_index
	#-------------------------------------------------------------
	@staticmethod
	def search_errors_at_read_index(err_tab_file, read_index):
		errors = []
		last_pos = None
		
		while(True):
			last_pos = err_tab_file.tell()
			line = err_tab_file.readline()
			if not line: #eof
				return errors
			#print "lililouloulala:", line.strip()
			#print line
			#print line.strip()
			columns = line.strip().split("\t")
			#print int(columns[0]) == read_index
			if int(columns[0]) == read_index:
				#print "added"
				errors.append(GenericTest.ReadError(int(columns[1]), columns[2], columns[3]))
			else:
				#print "octect", sys.getsizeof(line)
				break
			
		#print read_index, errors
		#if last_pos:
		err_tab_file.seek(last_pos)
		return errors
		"""
		errors = []
		for i in range(0, len(noerr_read)):
			if noerr_read[i] != err_read[i]:
				error = ReadError(i, noerr_read[i], err_read[i])
				errors.append(error)
		return errors
		"""

	#-------------------------------------------------------------
	# * check_true_correction
	#-------------------------------------------------------------
	@staticmethod
	def check_true_correction(corrected_read, errors):
		error_corrected = 0
		for error in errors:
			if corrected_read[error.pos] == error.noerr_letter:
				error_corrected += 1
		return error_corrected

	#-------------------------------------------------------------
	# * check_false_correction
	#-------------------------------------------------------------
	@staticmethod
	def check_false_correction(err_read, corrected_read, errors):
		false_error_corrected = 0;
		errors_pos = []
		for error in errors:
			errors_pos.append(error.pos)
		#---
		for i in range(0, len(corrected_read)):
			if (err_read[i] != corrected_read[i]):
				if not (i in errors_pos):
					false_error_corrected += 1
					
		return false_error_corrected
		

#-------------------------------------------------------------
# ** BloocooCorrection
#-------------------------------------------------------------
class BloocooTest:
	
	#-------------------------------------------------------------
	# * execute
	#-------------------------------------------------------------
	@staticmethod
	def execute(params, err_tab_filename, corrected_tab_filename, corrected_tab_full_filename):
		err_tab_file = open(err_tab_filename, "r")
		corrected_tab_file = open(corrected_tab_filename, "r")
		#---
		true = int(commands.getstatusoutput("wc -l " + err_tab_filename + "| awk '{print $1}'")[1])#int(commands.getstatusoutput("wc -l " + err_tab_filename)[1].split(" ")[0])
		#---
		t = time.time()
		os.system("sort " + err_tab_filename + "> sorted_"+err_tab_filename)
		os.system("sort " + corrected_tab_filename + "> sorted_"+corrected_tab_filename)
		os.system("sort " + corrected_tab_full_filename + "> sorted_"+ corrected_tab_full_filename)
		print "Time for sorting both tabs:", time.time() - t
		#---
		t = time.time()
		command = "comm -1 " + "sorted_"+err_tab_filename + " " + "sorted_"+corrected_tab_filename + " > comm.temp"
		os.system(command)
		predicted = int(commands.getstatusoutput("wc -l comm.temp | awk '{print $1}'")[1])
		print "Time to compare (predicted):", time.time() - t
		#---
		command = "comm -12 " + "sorted_"+err_tab_filename + " " + "sorted_"+corrected_tab_filename + " > TP"
		os.system(command)
		tp = int(commands.getstatusoutput("wc -l TP| awk '{print $1}'")[1])
		command = "comm -13 " + "sorted_"+err_tab_filename + " " + "sorted_"+corrected_tab_filename + " > FP"
		os.system(command)
		command = "comm -23 " + "sorted_"+err_tab_filename + " " + "sorted_"+corrected_tab_filename + " > FN"
		os.system(command)
		command = "join -t \";\" -1 1 -2 1 "   + "sorted_"+ corrected_tab_full_filename  + " FP > FPfull"
		os.system(command)

		command = "grep twosided FPfull | wc -l | awk '{print $1}'"
		FP0 = int(commands.getstatusoutput(command)[1])
		command = "grep aggressive FPfull | wc -l | awk '{print $1}'"
		FP1 = int(commands.getstatusoutput(command)[1])
		command = "grep vote FPfull | wc -l | awk '{print $1}'"
		FP2 = int(commands.getstatusoutput(command)[1])
		command = "grep multi FPfull | wc -l | awk '{print $1}'"
		FP3 = int(commands.getstatusoutput(command)[1])
		command = "grep \;side FPfull | wc -l | awk '{print $1}'"
		FP4 = int(commands.getstatusoutput(command)[1])

		command = "grep twosided sorted_"+corrected_tab_full_filename + " | wc -l | awk '{print $1}'"
		TC0 = int(commands.getstatusoutput(command)[1])
		command = "grep aggressive sorted_"+corrected_tab_full_filename + " | wc -l | awk '{print $1}'"
		TC1 = int(commands.getstatusoutput(command)[1])
		command = "grep vote sorted_"+corrected_tab_full_filename + " | wc -l | awk '{print $1}'"
		TC2 = int(commands.getstatusoutput(command)[1])
		command = "grep multi sorted_"+corrected_tab_full_filename + " | wc -l | awk '{print $1}'"
		TC3 = int(commands.getstatusoutput(command)[1])
		command = "grep \;side sorted_"+corrected_tab_full_filename + " | wc -l | awk '{print $1}'"
		TC4 = int(commands.getstatusoutput(command)[1])
		
		return true, tp, predicted-tp, FP0, FP1, FP2 ,FP3, FP4, TC0, TC1, TC2 ,TC3, TC4

	@staticmethod
	def print_stats(recall, precision, FP0, FP1, FP2 ,FP3, FP4, TC0, TC1, TC2 ,TC3, TC4):
		print "----------------------------------------"
		print "- Bloocoo stats"
		print "-\tRecall: " + str(recall) + "%"
		print "-\t\tFP twosided :", FP0, "\t/ ", TC0
		print "-\t\tFP aggressive :", FP1, "\t/ ", TC1
		print "-\t\tFP vote :", FP2, "\t/ ", TC2
		print "-\t\tFP multivote :", FP3, "\t/ ", TC3
		print "-\t\tFP side :", FP4, "\t/ ", TC4
		
		print "-\tPrecision: " + str(precision) + "%"

		if TC0 == 0:
			precision = 100
		else:
			precision = '%.2f' % ((float(TC0-FP0)/TC0)*100)
		
		print "-\t\t Prec twosided: " + str(precision) + "%"

		if TC1 == 0:
			precision = 100
		else:
			precision = '%.2f' % ((float(TC1-FP1)/TC1)*100)
		
		print "-\t\t Prec aggressive: " + str(precision) + "%"

		if TC2 == 0:
			precision = 100
		else:
			precision = '%.2f' % ((float(TC2-FP2)/TC2)*100)
		
		print "-\t\t Prec vote: " + str(precision) + "%"


		if TC3 == 0:
			precision = 100
		else:
			precision = '%.2f' % ((float(TC3-FP3)/TC3)*100)
		
		print "-\t\t Prec multivote: " + str(precision) + "%"


		if TC4 == 0:
			precision = 100
		else:
			precision = '%.2f' % ((float(TC4-FP4)/TC4)*100)
		
		print "-\t\t Prec side: " + str(precision) + "%"

