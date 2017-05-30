#!/usr/bin/env python
__author__ = 'Jessica_Bryant'
__email__ = 'jessawbryant@gmail.com'


"""
This program, for each amino acid sequence in the HOT assemblies dataset, maps the taxa composition and calculates the
GC, average nitrogen and carbon content of encoded protein side chains and average aa molecular weight for each
amino acid sequence.


Need to activate a virtual environment: source ~/DEEP_SEQUENCING/ENV/bin/activate
usage: python /mnt/lysine/jbryant/DEEP_SEQUENCING/codon_bias/GC_and_ARSC_from_abundances.py all_mg_genes.fna
"""

import sys
import re
import pandas as pd
from Bio import SeqIO
import codon_table as ct
import ntpath
from Bio import SeqUtils


def combine_samples(sample):
    '''
    This script takes a raw sample name and returns the final name if the sample came from a combined run.
    '''


    #samples that should be combined
    if(sample == 'HOT229_3_0500m'): sample_out = 'HOT229_0500m'
    elif(sample == 'HOT229_2_0500m'): sample_out = 'HOT229_0500m'
    elif(sample == 'HOT229_1_0500m'): sample_out = 'HOT229_0500m'
    elif(sample == 'HOT232_1_0500m'): sample_out = 'HOT232_0500m'
    elif(sample == 'HOT232_2_0500m'): sample_out = 'HOT232_0500m'
    elif(sample == 'HOT232_3_0500m'): sample_out = 'HOT232_0500m'
    elif(sample == 'HOT237_3_0500m'): sample_out = 'HOT237_0500m'
    elif(sample == 'HOT237_2_0500m'): sample_out = 'HOT237_0500m'
    elif(sample == 'HOT237_1_0500m'): sample_out = 'HOT237_0500m'
    elif(sample == 'HOT237_3_0025m'): sample_out = 'HOT237_0025m'
    elif(sample == 'HOT237_1_0025m'): sample_out = 'HOT237_0025m'
    elif(sample == 'HOT236_1_0500m'): sample_out = 'HOT236_0500m'
    elif(sample == 'HOT236_2_0500m'): sample_out = 'HOT236_0500m'
    elif(sample == 'HOT234_1_0200m'): sample_out = 'HOT234_0200m'
    elif(sample == 'HOT234_2_0200m'): sample_out = 'HOT234_0200m'
    elif(sample == 'HOT229_2_0025m'): sample_out = 'HOT229_0025m'
    elif(sample == 'HOT229_1_0025m'): sample_out = 'HOT229_0025m'
    elif(sample == 'HOT237_2_0770m'): sample_out = 'HOT237_0770m'
    elif(sample == 'HOT237_1_0770m'): sample_out = 'HOT237_0770m'
    elif(sample == 'HOT238_2_0500m'): sample_out = 'HOT238_0500m'
    elif(sample == 'HOT238_1_0500m'): sample_out = 'HOT238_0500m'

    else: sample_out = sample

    return(sample_out)


if __name__ == '__main__':

    fna_handle = sys.argv[1]
    errorfile_handle = open('~/jbryant/DEEP_SEQUENCING/codon_bias/' + ntpath.basename(fna_handle) +
                            '.errors', 'w')


    outfile = open(fna_handle + '_GC_ARSC.txt', 'w')
    outfile.writelines('\t'.join(['gene.id', 'sample', 'GC', 'SCU', 'N_ARSC', 'C_ARSC', 'N/C', 'av_MW',
                                  'codon_choice_rank', 'gene_length']) + '\n')


    # run through header of each unique DNA sequence in the .fna file, grab GC and calculate average amino
    #  acid usage

    fna_handle_open = open(fna_handle, "rU")
    gene_generator = SeqIO.parse(fna_handle_open, "fasta")

    timer = 0
    while True:
        try:
            gene = gene_generator.next()
            timer+=1
            print timer, gene.id

            #skip genes that are < 100bps long
            if len(gene.seq)<100:
                error2 = gene.id.strip() + '\t' + 'less than 100bps' + '\n'
                errorfile_handle.writelines(error2)
                continue

            # skip gene if it contains an internal stop codon
            if ct.scan_for_stop_codons(gene.seq) == 'True':
                error2 = gene.id.strip() + '\t' + 'internal stop codon' + '\n'
                errorfile_handle.writelines(error2)
                continue

            # calculate values
            sample = combine_samples(re.search('(HOT.*?m)', gene.id).group(1))
            GC = str(SeqUtils.GC(gene.seq[0:-3]))
            av_ARSC, av_MW, av_C_ARSC, N_C_ratio = ct.ARSC_MW_from_nucleotides(ct.screen_out_ambiguous_codons(gene.seq))
            SCU, codon_choice_rank = ct.calculate_SCU(gene, errorfile_handle)
            gene_length = str(len(gene.seq[0:-3]))
            outfile.writelines('\t'.join([gene.id, sample, GC, SCU, av_ARSC, av_C_ARSC, N_C_ratio, av_MW, codon_choice_rank, gene_length]) + '\n')


        except StopIteration as e:
            break
