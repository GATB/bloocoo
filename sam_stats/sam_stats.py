#!/usr/bin/env python
# requires: python >=2.7, pysam and pyfasta (e.g. pip install pysam pyfasta) don't have root? read this: http://stackoverflow.com/questions/20585218/install-python-package-without-root-access

k = 21
mapped_reads = "merged_tags_mapped.sam"
#mapped_reads = "bloocoo-corrected_mapped.sam"
reference = 'cp31.fasta'
assembly = '../minia/minia-k21-m4.contigs.fa'

import pysam
import pyfasta
import string
from collections import Counter

samfile = pysam.Samfile(mapped_reads, "r" )
ref_file = pyfasta.Fasta(reference)
ref_seq = str(ref_file[ref_file.keys()[0]])
ref_len = len(ref_seq)
print "reference length",ref_len
mapped_kmer_coverage = [0]*(ref_len+k-1)
ref_kmers = set([ref_seq[i:i+k] for i in xrange(len(ref_seq)-k+1)])
print "how many distinct kmers:",len(ref_kmers)
mapped_kmer_count = Counter(ref_kmers)
read_kmer_count = Counter()
revcomp_trans = string.maketrans('ACTGN10', 'TGACN10')

def revcomp(s):
    complement = string.translate(s, revcomp_trans)
    return complement[::-1]

def normalize_kmer(s):
    return min(revcomp(s),s)

print "initialization done"

# -----------
# parse SAM file of reads aligned to the reference, to compute
# 1) the k-mer coverage at each reference position 
# 2) a k-mer count for each k-mers in read (we actually didn't need an alignment for that..)
# 3) some stats about errors in reads

def parse_cigar(cigar,f):
    correspondence = { 'M': 0, 'I': 1, 'D': 2}
    return sum([y for x,y in cigar if correspondence[f] == x] + [0])

read_counter = 0
nb_bases_in_reads, nb_M, nb_D, nb_I = 0,0,0,0
for read in samfile.fetch():
    ref_end = read.aend # last base of the read + 1
    if ref_end is None:
        continue
    ref_end -= 1
    ref_start = ref_end - read.alen + 1 # first base of the read
    assert(ref_start == read.pos)
    # here we count read kmers (not using any alignment data)
    for pos in xrange(len(read.seq)-k+1):
        true_kmer = normalize_kmer(read.seq[pos:pos+k])
        read_kmer_count[true_kmer] += 1
    # here we only take reference k-mers into account, as long as they're
    # withing the bounday of an aligned reads
    for pos in xrange(ref_start,ref_end-k+1):
        mapped_kmer_coverage[pos] += 1
        ref_kmer = normalize_kmer(ref_seq[pos:pos+k])
        mapped_kmer_count[ref_kmer] += 1
    read_counter += 1
    nb_M += parse_cigar(read.cigar,'M')
    nb_I += parse_cigar(read.cigar,'I')
    nb_D += parse_cigar(read.cigar,'D')
    nb_bases_in_reads += read.alen # len(positions) would give exactly nb_M
    if read_counter % 500000 == 0:
        print "%d reads processed" % read_counter

samfile.close()

print "most common kmers in reads",read_kmer_count.most_common(5)

import numpy
mean_mapped_kmer_coverage = numpy.mean(mapped_kmer_coverage)
print "--- statistics about all k-mers from reads aligned to the reference"
print "--- the __read k-mer count__ is defined as the number of time a kmer exactly appears in the reads"
print "--- the __mapped k-mer count__ is defined as the number of time a reference kmer has a read kmer that aligns to it (possibly with mismatches)"
print "--- for a given reference position, the __mapped k-mer coverage__ is the number read kmers that align (possibly with mismatches) to the reference kmer at this position"
print "mean mapped kmer-coverage %.2f" % mean_mapped_kmer_coverage
nb_locations_below_5 = len([x for x in mapped_kmer_coverage if x < 5])
nb_locations_below_2 = len([x for x in mapped_kmer_coverage if x < 2])
nb_locations_at_0 = len([x for x in mapped_kmer_coverage if x == 0])
print "number of reference positions where mapped kmer coverage < 5:",nb_locations_below_5
print "number of reference positions where mapped kmer coverage < 2:",nb_locations_below_2
print "number of reference positions where mapped kmer coverage = 0:",nb_locations_at_0
print "number of kmers from the reference that are not seen in the reads:",len([x for x in ref_kmers if x not in read_kmer_count ])
print "number of kmers from the reference that have mapped kmer count < 2:",len([x for x in ref_kmers if mapped_kmer_count[x] < 2])

def percent(v):
    return " (%.2f %%)" % (100.0*v/nb_bases_in_reads)

print "number of bases in reads:",nb_bases_in_reads
print "among them, bases matched:",nb_M,percent(nb_M)
print "            deletions:",nb_D,percent(nb_D)
print "            insertions:",nb_I,percent(nb_I)
# -----------
# if an assembly is provided, parse contigs to see how covered are the 4 possible extensions at the extremities 

if not assembly:
    exit("Finished. Provide contigs for extended analysis.")

def right_extension(kmer):
    return [ normalize_kmer(kmer[1:] + y) for y in "ACTG" ]

def left_extension(kmer):
    return [ normalize_kmer(y + kmer[:-1]) for y in "ACTG" ]

contigs = pyfasta.Fasta(assembly)
could_have_been_extended = 0
no_chance_to_be_extended = 0
f = open("contigs_stats.txt","w")
g = open("critical_kmers.txt","w")
for contig_name in contigs:
    contig = str(contigs[contig_name])
    start_kmer = contig[:k]
    end_kmer = contig[-k:]
    f.write("contig %s:\n" % contig_name)
    for dir_extension, dir, kmer in [(left_extension, "left", start_kmer), (right_extension, "right", end_kmer)]:
        max_covs = []
        for count_method,count_dict in [("mapped",mapped_kmer_count),("read",read_kmer_count)]:
            max_cov, max_cov_kmer = max([(count_dict[x],x) for x in dir_extension(start_kmer) ])
            if count_method == "read":
                max_covs += [(max_cov, max_cov_kmer)]
            f.write("\t max %s kmer count in %s extension: %d\n" % (count_method, dir, max_cov))
        m = 4
        cov, kmer = max(max_covs)
        if cov >= m:
            could_have_been_extended += 1
        else:
            g.write("%s %d\n" % (kmer,cov))
            no_chance_to_be_extended += 1

print "parsed %d contigs, for %d extremities there was enough coverage to extend, for %d extremities there weren't" % (len(contigs), could_have_been_extended, no_chance_to_be_extended)
