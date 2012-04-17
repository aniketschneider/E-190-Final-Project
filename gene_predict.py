#!/usr/bin/env python
#
# gene_predict.py
#
# Aniket Schneider
# April 2012
#

from Bio import SeqIO
from collections import defaultdict

model_order = 5

def read_genes(training_filename):
    '''
    Generate a list of coding sequences from a FASTA file containing
    one or more genes.  All header information is discarded.
    '''
    genes = []

    for seq_record in SeqIO.parse(training_filename, "fasta"):
        genes.append(seq_record.seq)

    return genes

def count_patterns(sequence_list):
    fwd_counts = (defaultdict((lambda: 0)),
                  defaultdict((lambda: 0)),
                  defaultdict((lambda: 0)))
    rev_counts = (defaultdict((lambda: 0)),
                  defaultdict((lambda: 0)),
                  defaultdict((lambda: 0)))
    
    for gene in sequence_list:
        rev_gene = gene.reverse_complement()
        for i in xrange(len(gene) - model_order + 2):
            frame = (i + model_order) % 3
            fragment = str(gene[i:i + model_order + 1])
            fwd_counts[frame][fragment] += 1

            rev_fragment = str(rev_gene[i:i + model_order + 1])
            rev_counts[frame][rev_fragment] += 1

    print fwd_counts[0]


if __name__ == "__main__":
    genes = read_genes("Data/short_gene.fasta")
    build_model(genes)

