
import os
from os import listdir
from os.path import exists, isfile, join, basename, splitext, getsize
import shutil
from math import ceil


#fonction permettant d'executer un test avec les variables parametrables comme arguments
def execute_test(p):
	os.system("python ../../main.py " + p["result_filename_prefix"] + " " + str(p["genome_size"]) + " " + str(p["reads_size"]) +\
	" " + str(p["reads_count"]) + " " + str(p["kmer_size"]) + " " + str(p["coverage_threshold"]) + " " +\
	str(p["error_rate"]) + " " + str(p["nb_iter"]) + " " + str(p["nb_kmer_checked"]) + " " + p["genome_filename"]) 


#fonction permettant de creer un graphe en prennant comme donnee
#le dernier tabbed file cree par les appels successif a execute_test()
def execute_graph(p):
	if not exists("test_result/tabs"):
		os.mkdir("test_result/tabs")
	if not exists("test_result/graphs"):
		os.mkdir("test_result/graphs")
	#--- get the last result file added
	tab_filename = None
	filenames = listdir("test_result")
	for filename in filenames:
		complete_filename = join("test_result", filename)
		if isfile(complete_filename):
			tab_filename = complete_filename
			break
	graph_filename = splitext(basename(tab_filename))[0] + ".png"
	graph_filename = join("test_result/graphs", graph_filename)
	#---
	command = "Rscript " + join("../../Rscripts", p["R_script_filename"]) + " " + tab_filename + " " + graph_filename
	os.system(command)
	#--- Move the last result file to tabs dir
	shutil.move(tab_filename, "test_result/tabs")


#Prepare the tests for a new genome
#Currently this method:
#	- determine the genome size
#	- determine the number of reads needed depending of cover
def setup(p):
	if p["genome_filename"] != "":
		p["genome_size"] = getsize(join("../../genomes", p["genome_filename"]))
	
	#Calc the reads_count depending of cover
	offset = (p["reads_size"]-p["kmer_size"]+1) / float(p["reads_size"])
	p["reads_count"] = int(ceil((p["cover"]*p["genome_size"]) / (p["reads_size"]*offset)))



