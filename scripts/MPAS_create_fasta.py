#!/usr/bin/env python

import os
import pandas as pd
import csv
import sys


# setup command two line arguments:
# 1. sequencing_report.xlsx
# 2. targets.tsv

sequencing_report: sys.argv[1]
targets: sys.argv[2]


# 1. read in all_data tab and metadata tabs from sequening_report.xlsx file as two separate dataframes and read in the targets.tsv file

val = pd.read_excel(sys.argv[1], sheet_name='all_data')
meta = pd.read_excel(sys.argv[1], sheet_name='metadata')
targets = pd.read_table(sys.argv[2])

# 2. clean up each dataframe to only include specified fields

#all_data DF: only indluce Index, abundance, sequence number, subject, percent identity, and sequence
seqs = val[['Index', 'abundance', 'primer_name', 'sequence_number', 'subject','percent_identity', 'observed_sequence']]

#metadata DF: only include CSID, Index, Tick_Genus_species, Lifestage, state
meta_sub = meta[['Index', 'CSID', 'Morphological_Ectoparasite_Genus_Species', 'Lifestage', 'State']]

#targets DF:
targets_sub = targets[['primer_name', 'ref_sequence_name', 'sequence']]


# 3. Join seq and meta_sub into one DF based on Index
meta_sub_seqs = pd.merge(seqs, meta_sub, on='Index').reindex(columns=['Index', 'CSID', 'abundance', 'primer_name', 'sequence_number', 'subject', 'percent_identity', 'Morphological_Ectoparasite_Genus_Species', 'Lifestage', 'State', 'observed_sequence'])

# 4. Keep all sequences with abundance greter than 50
seq_50 = meta_sub_seqs[meta_sub_seqs['abundance'] >=50]



# 6. Group representative sequences found in each state and tick species and create CSV
seq_50_group = seq_50.groupby(['State','sequence_number', 'Morphological_Ectoparasite_Genus_Species', 'Lifestage', 'subject', 'primer_name', 'observed_sequence']).size().reset_index(name='n')
seq_50_group_reorder = seq_50_group[['State', 'Morphological_Ectoparasite_Genus_Species', 'Lifestage', 'sequence_number', 'subject', 'primer_name', 'n', 'observed_sequence']]



# 7. Funtions to Convert sample CSV and ref seq CSV files to fasta files

def write_fasta(in_seq, out_seq):

    try:
        fin = open(in_seq, 'r')
        fout = open(out_seq,'w')
    except:
        return -1

    with fin:
        reader = csv.DictReader(fin, skipinitialspace=True)
        for line in reader:
            state = line['State']
            genus = line['Morphological_Ectoparasite_Genus_Species']
            lifestage = line['Lifestage']
            number = line['sequence_number']
            tot = line['n']
            seq = line['observed_sequence']
            print(f'>{state}_{genus}_{lifestage}_seq{number}_n={tot}\n{seq}', file=fout)
    return

def write_ref_fasta(infile, outfile):
    try:
        f_in = open(infile, 'r')
        f_out = open(outfile,'w')
    except:
        return -1

    with f_in:
        read = csv.DictReader(f_in, skipinitialspace=True)
        for line in read:
            primer = line['primer_name']
            name = line['ref_sequence_name']
            ref_seq = line['sequence']
            print(f'>{primer}_{name}\n{ref_seq}', file=f_out) 
    return


# 8.Create csv files and fasta files for samples and ref sequences based on primers

for primers, seqs in seq_50_group_reorder.groupby(['primer_name']):
    seqs.to_csv(f'{primers}.csv', index=False)
    write_fasta(f'{primers}.csv', f'{primers}_fasta.fasta')

for primer, seq in targets_sub.groupby(['primer_name']):
    seq.to_csv(f'{primer}_ref.csv', index=False)
    write_ref_fasta(f'{primer}_ref.csv', f'{primer}_ref_fasta.fasta')


# 9. Concatenate the reference seq fasta files with the sample fasta files based on primer

for primer, seq in targets_sub.groupby(['primer_name']): #put primer name into variable
    primer_name = primer
    with open(f'{primer_name}_all.fasta', 'w') as outfile: #open .fasta files for each primer name
        for file in os.listdir(os.curdir): #files in the curren direcotry
            if file.startswith(f'{primer_name}') and file.endswith('.fasta'): #if files start with the primer name and ends with end with .fasta
                with open(file) as infile: #open the selected files
                    contents = infile.read() #read the files and put into variables
                    outfile.write(contents) #write the contents of those files into the outfiles
