#!/usr/bin/env python
#
# gene_predict.py
#
# Aniket Schneider
# April 2012
#

from Bio import SeqIO
from collections import defaultdict
import math
from fractions import Fraction
import sys
from orf_finder import find_all_orfs
import re
import random

model_order = 5
MIN_ORF_LENGTH = 90
MAX_OVERLAP = 30

def read_genes(filename):
    '''
    Generate a list of sequence objects from a FASTA file containing
    one or more genes.  All header information is discarded.
    '''
    genes = []

    try:
        with open(filename) as infile:
            for seq_record in SeqIO.parse(infile, "fasta"):
                genes.append(seq_record)
    except:
        print "Error reading file:", filename

    return genes

def count_patterns(sequence_list, pattern_length):
    '''
    Return dictionaries mapping sequence fragments of fixed length to number of
    occurrences. Return value is a pair of two triples of dicts, the first of
    which is the counts for each reading frame in the forward direction, and the
    second of which is the counts for each reading frame in the reverse
    direction.

    Expects a list of sequence objects.
    '''
    fwd_counts = (defaultdict((lambda: 0)),
                  defaultdict((lambda: 0)),
                  defaultdict((lambda: 0)))
    rev_counts = (defaultdict((lambda: 0)),
                  defaultdict((lambda: 0)),
                  defaultdict((lambda: 0)))
    
    sequence_list = [rec.seq for rec in sequence_list]

    for gene in sequence_list:
        rev_gene = gene.reverse_complement()
        for i in xrange(len(gene) - pattern_length + 1):
            frame = (i + pattern_length) % 3
            fragment = str(gene[i:i + pattern_length])
            fwd_counts[frame][fragment] += 1

            rev_fragment = str(rev_gene[i:i + pattern_length])
            rev_counts[frame][rev_fragment] += 1

    return (fwd_counts, rev_counts)

def counts_to_mm(counts):
    '''
    Convert a dictionary of counts for various string fragments into a Markov
    Model mapping contexts to the probability of each emission given that 
    context.
    '''
    # break each pattern into prefix and last character, merge prefixes
    count_tree = defaultdict( (lambda: defaultdict(lambda: 0)) )
    for key in counts:
        count_tree[key[:-1]][key[-1]] += counts[key]

    # convert counts into log probabilities
    model = defaultdict( (lambda: defaultdict(lambda: 0.) ))

    for key in count_tree:
        total_count = sum(count_tree[key][lastch] for lastch in count_tree[key])
        assert total_count > 0

        for lastch in count_tree[key]:
            assert count_tree[key][lastch] > 0
            model[key][lastch] = math.log(
                    float(count_tree[key][lastch])/float(total_count))
            
    return model


def build_fmm(genes, order):
    '''
    Build a fixed-order markov model based on a set of gene sequences. The 
    model consists of six dicts, three mapping contexts to emissions for 
    forward reading frames and three for the reverse reading frames.
    '''
    (fwd_counts, rev_counts) = count_patterns(genes, order + 1)
    fwd_mms = tuple(counts_to_mm(d) for d in fwd_counts)
    rev_mms = tuple(counts_to_mm(d) for d in rev_counts)

    return (fwd_mms, rev_mms)


def score_sequence_fmm(seq, model, order):
    '''
    Predict whether a sequence is likely to be a gene or not.

    Score the sequence using a frame offset of 0, 1, and 2 in each direction,
    and predict a gene if the 0 offset in the forward direction scores the
    highest.
    '''
    length = len(seq)
    score = [Fraction(1) for i in xrange(6)]

    # fo = 'frame offset'
    for fo in [0,1,2,3,4,5]:
        for i in xrange(length - order):
            # fo % 2 will cycle through both the forward and reverse models
            # the order is a bit odd, but the correct forward frame is always
            # first
            frame_model = model[fo % 2][(i + order + 1 + fo) % 3]

            # the model contains strings, so seq must be explicitly converted 
            # to strings in case another sequence container is used.
            prefix = str(seq[i:i+order])
            emission = str(seq[i+order])

            score[fo] += frame_model[prefix][emission]
            if frame_model[prefix][emission] == float('-inf'):
                print prefix, emission, frame_model[prefix][emission]

    if score[0] == max(score):
        return (True, score[0])
    else:
        return (False, None)

def filter_overlaps(orfs):
    '''
    For a set of predicted coding regions, filter out regions that overlap by
    more than 30 base pairs.  For each overlap conflict, the longer predicted 
    ORF is always chosen.  Room for improvement here - could choose overlapping
    regions based on scores, or scores for overlapping regions.

    BUG: A chain of overlapping regions of increasing length will throw away
    earlier short regions despite the fact that they do not overlap with the
    long region two or more iterations later.  This may or may not be
    significant.
    '''
    # sort by increasing start position
    orfs.sort(key=lambda orf: orf['start'])

    orfs_no_overlap = []

    current = orfs[0]
    for orf in orfs[1:]:
        if overlap( (current['start'],current['end']), 
                    (orf['start'],orf['end']) ) > MAX_OVERLAP:
            if orf['score']/orf['length'] > current['score']/current['length']:
                current = orf
        else:
            orfs_no_overlap.append(current)
            current = orf

    return orfs_no_overlap
        
def overlap(intv1, intv2):
    '''
    Determine whether the end of intv1 overlaps the start of intv2
    '''
    (s1,e1) = intv1
    (s2,e2) = intv2

    if e1 > s2 and e2 > s1:
        return e1 - s2
    else:
        return 0

def orf_intervals(orfs):
    '''
    Transform a list of orfs into a sorted list of intervals corresponding to
    the positions of the genes.
    '''
    pred_genes = filter(lambda orf: orf['predicted'], orfs)
    intervals = map(lambda orf: (orf['start'],orf['end']), pred_genes)

    # sort by increasing start
    intervals.sort(key=lambda (a,b): a)

    return intervals

def parse_location(desc):
    '''
    Searches for the location field in the sequence record description and
    extracts the start and end point.
    '''
    match = re.search(r'\[location=[^\d]*([\d]+)\.\.([\d]+)\)?\]', 
            desc)
    if match:
        start = int(match.group(1))
        end = int(match.group(2))
        return (start,end)
    else:
        return None

def gene_intervals(genes):
    '''
    Turn a set of gene records in biopython sequence record format into a set of
    intervals.
    '''
    # genes are expected to be in biopython sequence record format, with
    # annotation as to their location.
    locs = map(lambda rec: parse_location(rec.description), genes)
    locs = filter(lambda x: x != None, locs)

    # sort by increasing start position
    locs.sort(key=lambda (a,b): a)

    return locs

def evaluate_results(orfs, ref_genes):
    '''
    Compare a set of predicted genes to a reference set of curated coding
    regions.
    '''
    pred_intervals = set(orf_intervals(orfs))
    ref_intervals = set(gene_intervals(ref_genes))

    pred_stops = set(map(lambda (a,b): b, pred_intervals))
    ref_stops = set(map(lambda (a,b): b, ref_intervals))
        
    print "Predicted", len(pred_intervals), "coding regions"
    print "Reference data contained", len(ref_intervals), "coding regions"
    print "Got", len(pred_intervals.intersection(ref_intervals)), "exactly correct"
    print "Predicted", len(pred_stops.intersection(ref_stops)), "correct ORFs (including incorrect start positions)"
    print "Predicted", len(pred_stops.difference(ref_stops)), "false positives"
    print "Missed", len(ref_stops.difference(pred_stops)), "ORFs"

if __name__ == "__main__":
    def usage():
        print "{0} -t <training datafile> [-r <reference data>] -i <input datafile>".format(sys.argv[0])
        print "\tTraining datafile should be a set of coding DNA sequences in FASTA format"
        print "\tInput datafile should be a DNA sequence in FASTA format containing one or more ORFs"
        print "\tReference data should be a curated set of coding sequences found in the input datafile to compare the results to"

    if '-h' in sys.argv or '--help' in sys.argv or '-help' in sys.argv:
        usage()
        exit(0)

    if len(sys.argv) < 4 or "-t" not in sys.argv or "-i" not in sys.argv:
        usage()
        exit(1)

    if "-r" in sys.argv:
        ropt = sys.argv.index("-r")
        ref_filename = sys.argv[ropt + 1]
        ref_genes = read_genes(ref_filename)
        use_ref = True
    else:
        use_ref = False

    topt = sys.argv.index("-t")
    training_filename = sys.argv[topt + 1]
    training_genes_all = read_genes(training_filename)

    # randomly throw out some percentage of the training set
    random.seed(1)
    training_genes = training_genes_all
    #for i in xrange(len(training_genes_all)):
        #if random.random() > 1.0:
            #training_genes.append(training_genes_all[i])

    iopt = sys.argv.index("-i")
    input_filename = sys.argv[iopt + 1]
    input_sequences = read_genes(input_filename)

    model = build_fmm(training_genes, model_order)

    orfs = []
    for rec in input_sequences:
        seq = str(rec.seq)
        orfs.extend(find_all_orfs(seq, minlen=MIN_ORF_LENGTH))
        if len(orfs) == 0:
            print "No ORFs found in: {0}".format(rec.name)
        for orf in orfs:
            (prediction, score) = score_sequence_fmm(
                    orf['seq'], model, model_order)
            if prediction:
                orf['predicted'] = True
                orf['score'] = score
            else:
                orf['predicted'] = False

        raw_predictions = filter(lambda orf: orf['predicted'], orfs)
        predicted_genes = filter_overlaps(raw_predictions)
        '''
        for orf in predicted_genes:
            print orf['start'], '-', orf['end'], 
            print ('reverse strand' if orf['frame'] < 0 
                    else 'forward strand')
        '''


    if use_ref:
        evaluate_results(predicted_genes, ref_genes)
