#!/usr/bin/env python
# requires: python >=2.7, pysam and pyfasta>=0.8 (e.g. pip install --user pysam pyfasta) if it fails, read this: http://stackoverflow.com/questions/20585218/install-python-package-without-root-access
# also on genocluster, run with "python -S" else it will take the machine pyfasta (older version, not compatible)

import pysam
import pyfasta
import string
from collections import defaultdict 
import sys

if len(sys.argv) < 4:
    exit("args: cor.sam uncor.sam ref.fa")

def open_sambam(file):
    return pysam.Samfile(file, "rb" if file.endswith("bam") else "r")#, check_header = False, check_sq = False)

mapped_reads_corrected = open_sambam(sys.argv[1])
mapped_reads_uncorrected = open_sambam(sys.argv[2])
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
good_corrections_per_ref, bad_corrections_per_ref = defaultdict(int), defaultdict(int)
better_reads, worse_reads = 0,0
better_reads_per_ref, worse_reads_per_ref = defaultdict(int), defaultdict(int)
uncor_iterator = mapped_reads_uncorrected


def report():
    print "---"
    print "done, %d corrected reads processed" % cor_read_counter
    print "number of bases in corrected/uncorrected reads:",nb_bases_in_cor_reads,"/", nb_bases_in_uncor_reads
    print "good/bad base corrections", good_corrections, "/", bad_corrections
    print "better/worse reads", better_reads, "/", worse_reads
    print "per reference:"
    print "refname\tb_good\tb_bad\tr_better\tr_worse"
    for refname in good_corrections_per_ref.keys():
        print refname,"\t",good_corrections_per_ref[refname], "\t", bad_corrections_per_ref[refname], "\t",better_reads_per_ref[refname], "\t", worse_reads_per_ref[refname]
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
    try:
        uncor_read = uncor_iterator.next()
        if uncor_read.alen is not None:
            nb_bases_in_uncor_reads += uncor_read.alen 
        skip=0
        while cor_read.qname != uncor_read.qname:
            uncor_read = uncor_iterator.next()
            if uncor_read.alen is not None:
                nb_bases_in_uncor_reads += uncor_read.alen 
            skip += 1
        if skip:
            print "skipped",skip,"uncorrected reads"
    except StopIteration:
        print "end of uncorr reads (skip %d, nb cor reads processed %d, was looking for cor read id %s)" % (skip, cor_read_counter, cor_read.qname)
        sys.exit(1) 
    #assert(cor_read.qname == uncor_read.qname)
    
    if uncor_read.alen is None:
        continue

    ref_end = cor_read.aend # last base of the read + 1
    if ref_end is None:
        continue
    ref_end -= 1
    ref_start = ref_end - cor_read.alen + 1 # first base of the read
    #assert(ref_start == cor_read.pos)

    # try to make sure both the corrected and uncorrected reads map at the same place
    if cor_read.reference_id != uncor_read.reference_id:
        continue
    if abs(cor_read.aend - uncor_read.aend) > cor_read.alen:
        continue

    refname = mapped_reads_corrected.getrname(cor_read.reference_id)
    cor_mismatches = get_mismatches_tuples(cor_read, refname)
    uncor_mismatches = get_mismatches_tuples(uncor_read, refname)
    good_corrections += len(uncor_mismatches - cor_mismatches)
    bad_corrections += len(cor_mismatches - uncor_mismatches)    
    good_corrections_per_ref[refname] += len(uncor_mismatches - cor_mismatches)
    bad_corrections_per_ref[refname] += len(cor_mismatches - uncor_mismatches)    
    #print cor_mismatches - uncor_mismatches, cor_mismatches, uncor_mismatches, cor_read.aend, uncor_read.aend

    if len(cor_mismatches) < len(uncor_mismatches):
        better_reads += 1
        better_reads_per_ref[refname] += 1
    elif len(cor_mismatches) > len(uncor_mismatches):
        worse_reads += 1
        worse_reads_per_ref[refname] += 1

report()
