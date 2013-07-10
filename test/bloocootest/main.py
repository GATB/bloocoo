
import sys, os
from libs.TestReadCorrection import TestReadCorrection

params = {
	"result_filename_prefix": sys.argv[1],
	"genome_size": int(sys.argv[2]),
	"reads_size": int(sys.argv[3]),
	"reads_count": int(sys.argv[4]),
	"kmer_size": int(sys.argv[5]),
	"coverage_threshold": int(sys.argv[6]),
	"error_rate": float(sys.argv[7]),
	"nb_kmer_checked": int(sys.argv[8]),
	"regenerate_reads": sys.argv[9],
	"genome_filename": "",
}

try:
	params["genome_filename"] = sys.argv[10]
except:
	pass
	
TestReadCorrection.main(params)

