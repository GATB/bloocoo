#!/usr/bin/env python
# requires: python >=2.7, pysam and pyfasta (e.g. pip install --user pysam pyfasta) if it fails, read this: http://stackoverflow.com/questions/20585218/install-python-package-without-root-access

import pysam
import pyfasta
import string
from collections import Counter
import sys

if len(sys.argv) < 4:
    exit("args: cor.sam uncor.sam ref.fa")

mapped_reads_uncorrected =  pysam.Samfile(sys.argv[1],"r")
mapped_reads_corrected =  pysam.Samfile(sys.argv[2],"r")
reference = pyfasta.Fasta(sys.argv[3])

print "reference number of sequences",len(reference)
print "initialization done"

# -----------
# parse SAM file of reads aligned to the reference

def parse_cigar(cigar,f):
    correspondence = { 'M': 0, 'I': 1, 'D': 2}
    return sum([y for x,y in cigar if correspondence[f] == x] + [0])

def get_mismatches_tuples(r, refname):
    print reference[refname]
    print r.query_sequence
    print r.get_aligned_pairs(matches_only=True)
    print  set(filter(lambda (x,y): r.query_sequence[x] != reference[refname][y], r.get_aligned_pairs(matches_only=True)))
    return set(filter(lambda (x,y): r.query_sequence[x] != reference[refname][y], r.get_aligned_pairs(matches_only=True)))

cor_read_counter = 0
nb_bases_in_cor_reads, nb_bases_in_uncor_reads = 0, 0
good_corrections, bad_corrections = 0,0
uncor_iterator = mapped_reads_uncorrected.fetch()

for cor_read in mapped_reads_corrected.fetch():
    if cor_read.alen is None:
        continue
    nb_bases_in_cor_reads += cor_read.alen 
    cor_read_counter += 1

    if cor_read_counter % 500000 == 0:
        print "%d corrected reads processed" % cor_read_counter

    # synchronine with uncorrected reads
    uncor_read = uncor_iterator.next()
    nb_bases_in_uncor_reads += uncor_read.alen 
    while cor_read.qname != uncor_read.qname:
        print uncor_read.qname, cor_read.qname
        uncor_read = uncor_iterator.next()
        nb_bases_in_uncor_reads += uncor_read.alen 
    assert(cor_read.qname == uncor_read.qname)

    ref_end = cor_read.aend # last base of the read + 1
    if ref_end is None:
        continue
    ref_end -= 1
    ref_start = ref_end - cor_read.alen + 1 # first base of the read
    assert(ref_start == cor_read.pos)

    refname = mapped_reads_corrected.getrname(cor_read.reference_id)
    cor_mismatches = get_mismatches_tuples(cor_read, refname)
    uncor_mismatches = get_mismatches_tuples(uncor_read, refname)
    print cor_mismatches
    good_corrections += len(uncor_mismatches - cor_mismatches)
    bad_corrections = len(cor_mismatches - uncor_mismatches)

print "done, %d corrected reads processed" % cor_read_counter
print "number of bases in corrected/uncorrected reads:",nb_bases_in_cor_reads,"/",nb_bases_in_uncor_reads
print "---"
print "good corrections", good_corrections
print "bad corrections", bad_corrections
