#!/usr/bin/python
# compare blast results, check reciprocal blast
# these blast results are in tables.
# see https://docs.google.com/document/d/1N1RexllnvWfXaJ5gfw8klEegDt2vkMlyam2jSC7XKvM/edit#bookmark=id.svttxn2t0wxe
import pandas as pd

def remove_end(acc_list):
    return [x.split('.')[0] for x in acc_list]


def main():
    
    # Important Files
    accessions_file = 'human_movement_disorder_prot_accessions_Feb2025.txt'
    human_protein_sequences = 'human_protein_sequences_199k.fa'
    blast_v_tardigrade = 'hum_v_exemplaris_prots.blastp' # hum_v_exemplaris_transcripts.txt
    blast_back = 'exemplaris_hits_v_human_all.blastp' # tard_transcript_hits_v_human.txt
    
    # get dictionary of tardigrade_seq ==> human top hit
    blastback_tophits = get_blastback_tophits(blast_back)
    
    # get dictionary of human accession# ==> gene ID
    gene_ids = get_geneids(human_protein_sequences)
    
    # compare the two blast resuults
    header, blast_info = reciprocal_blast_compare(blast_v_tardigrade, blastback_tophits, gene_ids)
    
    accessions = get_accession_list(accessions_file)
    
    print(header)
    for acc in accessions:
        print(blast_info[acc])


def get_accession_list(accessions_file):
    f = open(accessions_file,'r')
    accessions = [x.rstrip() for x in f.readlines()]
    f.close()
    return accessions

def get_blastback_tophits(blast_back):
    # make dictionary of tardigrade accession => backhit
    blast_back_results = pd.read_csv(blast_back, delimiter='\t')
    tardicolumn = blast_back_results.columns[0]
    humcolumn =  blast_back_results.columns[1]
    tard_accessions = blast_back_results[tardicolumn].values
    back_hits = blast_back_results[humcolumn].values 
    tard_accessions = remove_end(tard_accessions)
    back_hits = remove_end(back_hits)
    blastback_tophits = dict(zip(tard_accessions,back_hits))
    return blastback_tophits

def get_geneids(protein_seqs):
    # make dictionary of accession ==> GeneIDs
    f = open(protein_seqs,'r')
    fasta = f.readlines()
    f.close()
    gene_ids = {}
    for line in fasta:
        if line.startswith('>'):
            accession = line.split()[0].split('.')[0][1:]
    #         print(accession) # testing OK
            geneID = int(line.split('GeneID=')[1].split(']')[0]) 
    #         print(geneID) # testing OK
            gene_ids[accession] = geneID
    return gene_ids

def reciprocal_blast_compare(blast_v_tardigrade, blastback_tophits, gene_ids):
    
    # go through original file (human vs. tardigrade BLAST)
    # add a column of GENE_IDS
    # and add column of top blastback top_hits
    # add a column of gene IDS for top blastback hits
    # and another column of whether the geneID matches.
    
    f = open(blast_v_tardigrade,'r')
    tardigrade_hits = f.readlines()
    f.close()
    header = tardigrade_hits[0].rstrip() + '\t'
    header += '\t'.join(['human_geneID','top_backhit','top_back_geneID','match'])
    # print(header)
    
    blast_info = {} # key = human protein; val = info from BLAST searches
    
    for line in tardigrade_hits[1:]:
        line = line.rstrip()
        hum_prot = line.split()[0].split('.')[0]
        hum_geneID = gene_ids[hum_prot]
        tardihit = line.split()[1].split('.')[0]
        evalue = line.split()[2]
        
        if tardihit in blastback_tophits.keys():
            blastbacktop = blastback_tophits[tardihit]
            blastbacktopID = gene_ids[blastbacktop]
        else:
            blastbacktop = 'no blastback hit for ' + tardihit
            blastbacktopID = 'NA'
        
        if hum_geneID == blastbacktopID:
            match = 'yes'
        else:
            match = 'no'
        blast_info[hum_prot] = '\t'.join([hum_prot,tardihit,evalue,str(hum_geneID),
                         blastbacktop,str(blastbacktopID),match ])
        
    return header, blast_info
        
main()
