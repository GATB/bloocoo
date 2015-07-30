#!/usr/bin/env python
# requires: python >=2.7, pysam and pyfasta (e.g. pip install --user pysam pyfasta) if it fails, read this: http://stackoverflow.com/questions/20585218/install-python-package-without-root-access

import pysam
import pyfasta
import string
from collections import Counter
import sys

if len(sys.argv) < 4:
    exit("args: cor.sam uncor.sam ref.fa")

def open_sambam(file):
    return pysam.Samfile(file, "rb" if file.endswith("bam") else "r")#, check_header = False, check_sq = False)

mapped_reads_uncorrected = open_sambam(sys.argv[1])
mapped_reads_corrected = open_sambam(sys.argv[2])
reference = pyfasta.Fasta(sys.argv[3], key_fn=lambda key: key.split()[0]) #get before the first space

print "reference number of sequences",len(reference)
print "initialization done"

# -----------
# parse SAM file of reads aligned to the reference

def parse_cigar(cigar,f):
    correspondence = { 'M': 0, 'I': 1, 'D': 2}
    return sum([y for x,y in cigar if correspondence[f] == x] + [0])

def get_mismatches_tuples(r, refname):
    # return M tuples where the ref nucleotide doesnt match the read nucleotide (i.e. sequencing errors or snps)
    return set(filter(lambda (x,y): r.query_sequence[x] != reference[refname][y], r.get_aligned_pairs(matches_only=True)))

cor_read_counter = 0
nb_bases_in_cor_reads, nb_bases_in_uncor_reads = 0, 0
good_corrections, bad_corrections = 0,0
uncor_iterator = mapped_reads_uncorrected


def report():
    print "---"
    print "done, %d corrected reads processed" % cor_read_counter
    print "number of bases in corrected/uncorrected reads:",nb_bases_in_cor_reads,"/",nb_bases_in_uncor_reads
    print "good corrections", good_corrections
    print "bad corrections", bad_corrections
    print "---"


for cor_read in mapped_reads_corrected:
    if cor_read.alen is None:
        continue
    nb_bases_in_cor_reads += cor_read.alen 
    cor_read_counter += 1

    if cor_read_counter % 500000 == 0:
        print "%d corrected reads processed" % cor_read_counter
        print "so far",
        report()

    # synchronine with uncorrected reads
    uncor_read = uncor_iterator.next()
    while uncor_read.alen is None:
        uncor_read = uncor_iterator.next()
    nb_bases_in_uncor_reads += uncor_read.alen 
    while cor_read.qname != uncor_read.qname:
        uncor_read = uncor_iterator.next()
        while uncor_read.alen is None:
            uncor_read = uncor_iterator.next()
        nb_bases_in_uncor_reads += uncor_read.alen 
    #assert(cor_read.qname == uncor_read.qname)

    ref_end = cor_read.aend # last base of the read + 1
    if ref_end is None:
        continue
    ref_end -= 1
    ref_start = ref_end - cor_read.alen + 1 # first base of the read
    #assert(ref_start == cor_read.pos)

    refname = mapped_reads_corrected.getrname(cor_read.reference_id)
    cor_mismatches = get_mismatches_tuples(cor_read, refname)
    uncor_mismatches = get_mismatches_tuples(uncor_read, refname)
    good_corrections += len(uncor_mismatches - cor_mismatches)
    bad_corrections += len(cor_mismatches - uncor_mismatches)

report()
