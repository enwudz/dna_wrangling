#!/usr/bin/python
import pandas as pd
#blast_results = pd.read_csv('')

def main():
    
    ### build up a dictionary of sequences from a fasta file    
    # seqfile = 'human_protein_sequences_199k.fa'
    seqfile = 'h_exemplaris_ensembl_peps_Feb2025.fa'
    sequences, headers = build_seq_dictionary(seqfile)

    #### get a list of accession numbers for which we need sequences
    
    # OPTION 1: ... from a BLAST table
    blast_results = 'hum_v_exemplaris_prots.blastp'
    accessions = get_accessions_from_blast_table(blast_results, 'Tardigrade Hit')

    # OPTION 2:  ... from a FILE containing a list of accession numbers
    # accessions = get_accessions_from_file('human_movement_disorder_prot_accessions_Feb2025.txt')

    print_sequences(accessions, sequences, headers)

def get_accessions_from_file(fname):
    return [x.rstrip() for x in open(fname,'r').readlines()]

def print_sequences(accessions, sequences, headers):
    for accession in accessions:
        print(headers[accession] + '\n' + sequences[accession])

def get_accessions_from_blast_table(blast_results, colname='tardigrade'):
    results = pd.read_csv(blast_results, delimiter='\t')
    accessions = results[colname].values
    accessions = [x.split('.')[0] for x in accessions]
    accessions = list(set(accessions))
    return accessions

def build_seq_dictionary(seqfile):
    # build dictionaries of sequences and headers, keyed by accession number
    sequences = {}
    headers = {}
    with open(seqfile,'r') as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if '|' in line:
                    accession = line.split()[0].split('|')[1]
                elif 'gene:' in line:
                    accession = line.split('gene:')[1].split()[0]
                else:
                    accession = line.split()[0][1:].split('.')[0]
                sequences[accession] = ''
                headers[accession] = line
            elif len(line) > 0:
                sequences[accession] += line
    return sequences, headers

if __name__ == '__main__':
    main()
