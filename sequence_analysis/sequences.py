""" Importing sequences from Entrez directly """
from Bio import Entrez, SeqIO
from keys import Email, NCBI

Entrez.email = Email # Change to you NCBI Email
Entrez.api_key = NCBI # Change to your NCBI API

# If you do not have credentials, use the sample FASTA file in this folder

gene_name = "Sirt1" # Gene name
organisms = "Mus musculus" # Organisms to search
output_file = f"{gene_name}_{organisms}.fasta" # Output file

handle = Entrez.esearch( # Get conection for search
    db="nucleotide",
    term=f"{gene_name}[Gene Name] AND {organisms}[Organism] AND mRNA[Filter]",
    retmax=10)

result = Entrez.read(handle) # Get the XML / Use just once
handle.close()

print(f"Total IDs: {result}") # Show IDs

handle = Entrez.efetch( # Get conection to download sequences
    db="nucleotide",
    id=",".join(result["IdList"]),
    rettype="fasta",
    retmode="text"
)

records = list(SeqIO.parse(handle, "fasta")) # Parse sequences
handle.close() # Close handle

SeqIO.write(records, output_file, "fasta") # Save sequences


""" Genetic sequences manipulation """

def read_FASTA(filename):
    with open(filename) as file:
        return [(part[0].split('|'),
        part[2].replace('\n', ''))
        for part in
        [entry.partition('\n')
        for entry in file.read().split('>')[1:]]]

def print_FASTA_headers(filename): 
    with open(filename) as file:
        for line in file:
            if line[0] == '>':
                print(line[1:-1])

print_FASTA_headers(output_file)

print_FASTA_headers(output_file) # Print headers
seq = read_FASTA(output_file)[0][1] # Get the sequence


seq.count('A') # count number of occurrences of 'A'
len(seq) # length of the sequence
seq.find('T') # position of the first occurrence of 'T'
seq.upper() # convert to uppercase
seq.lower() # convert to lowercase
seq.strip() # remove whitespace
seq.split() # split the string into a list of strings
seq.join() # join a list of strings into a string
seq.startswith() # check if the string starts with a substring
seq.endswith() # check if the string ends with a substring
seq.find() # find the position of a substring
seq.index() # find the position of a substring
seq.count() # count the number of occurrences of a substring
seq.replace() # replace a substring with another
seq.strip() # remove whitespace
seq.split() # split the string into a list of strings
seq.join() # join a list of strings into a string

def validate_base_sequence(base_sequence, RNAflag = False):
    """Return True if the string base_sequence contains only upper- or lowercase
    T (or U, if RNAflag), C, A, and G characters, otherwise False"""
    seq = base_sequence.upper()
    return len(seq) == (seq.count('U' if RNAflag else 'T') +
    seq.count('C') +
    seq.count('A') +
    seq.count('G'))

help(validate_base_sequence)
print(validate_base_sequence(seq))

def gc_content(base_seq):
    """Return the percentage of G and C characters in base_seq"""
    assert validate_base_sequence(base_seq), \
    'argument has invalid characters'
    seq = base_seq.upper()
    return ((seq.count('G') + seq.count('C')) /
    len(seq))

print(gc_content(seq))

""" Translation """
import pprint as pp

nucleotide_names = dict((('A', 'adenine'),
    ('C', 'cytosine'),
    ('G', 'guanine'),
    ('T', 'thymine')
    ))

print(nucleotide_names)

RNA_codon_table = {
# Second Base
# U C A G
# U
'UUU': 'Phe', 'UCU': 'Ser', 'UAU': 'Tyr', 'UGU': 'Cys', # UxU
'UUC': 'Phe', 'UCC': 'Ser', 'UAC': 'Tyr', 'UGC': 'Cys', # UxC
'UUA': 'Leu', 'UCA': 'Ser', 'UAA': '---', 'UGA': '---', # UxA
'UUG': 'Leu', 'UCG': 'Ser', 'UAG': '---', 'UGG': 'Urp', # UxG
# C
'CUU': 'Leu', 'CCU': 'Pro', 'CAU': 'His', 'CGU': 'Arg', # CxU
'CUC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg', # CxC
'CUA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg', # CxA
'CUG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg', # CxG
# A
'AUU': 'Ile', 'ACU': 'Thr', 'AAU': 'Asn', 'AGU': 'Ser', # AxU
'AUC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser', # AxC
'AUA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg', # AxA
'AUG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg', # AxG
# G
'GUU': 'Val', 'GCU': 'Ala', 'GAU': 'Asp', 'GGU': 'Gly', # GxU
'GUC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly', # GxC
'GUA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly', # GxA
'GUG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'  # GxG
}

pp.pprint(RNA_codon_table)

def translate_RNA_codon(codon):
    return RNA_codon_table[codon]

translate_RNA_codon('UUU')

def translate_DNA_codon(codon):
    return translate_RNA_codon(codon.replace('T', 'U'))

seq_RNA = seq.replace('T', 'U') # Replace 'T' with 'U'
RNA_start = seq_RNA.find('AUG'); print(RNA_start)
processed_RNA = seq_RNA[RNA_start:]

def amino_acid_sequence(nucleotide_sequence):
    """Return the amino acid sequence of a nucleotide sequence, starting from the first AUG"""
    sequence = []
    for i in range(0, len(nucleotide_sequence), 3):
        codon = nucleotide_sequence[i:i+3]
        if len(codon) == 3:
            sequence.append(translate_RNA_codon(codon))
    return ''.join(sequence)

amino_acid = amino_acid_sequence(processed_RNA)
print(amino_acid)

#################################################
# Testes

# Sequence IDs from multiple files
def extract_gi_id(description):
    """Given a FASTA file description line, return its GenInfo ID if it has one"""
    if line[0] != '>':
        return None
    fields = description[1:].split('|')
    if 'gi' not in fields:
        return None
    return fields[1 + fields.index('gi')]

def get_gi_ids(filename):
    """Return a list of the GenInfo IDs of all sequences found in the file named filename"""
    with open(filename) as file:
        return [extract_gi_id(line) for line in file if line[0] == '>']

print(get_gi_ids('sirt1_Mus musculus.fasta'))

def get_gi_ids_from_files(filenames):
    """Return a list of the GenInfo IDs of all sequences found in the
    files whose names are contained in the collection filenames"""
    idlst = []
    for filename in filenames:
        idlst += get_gi_ids(filename)
    return idlst

print(get_gi_ids_from_files(['sirt1_Mus musculus.fasta']))

def rna_sequence_is_valid(seq):
    """Return True if the RNA sequence is valid"""
    for base in seq:
        if base not in 'UCAGucag':
            return False
    return True

def dna_sequence_contains_N(seq):
    """Return True if the DNA sequence contains N"""
    for base in seq:
        if base == 'N':
            return True
    return False

def count_gc_bases(seq):
    """Return the number of G and C bases in the sequence"""
    count = 0
    for base in seq:
        if base in 'GCgc':
            count += 1
    return count

def count_gc_percentage(seq):
    """Return the percentage of G and C bases in the sequence"""
    return count_gc_bases(seq) / len(seq) * 100

# DNA sequence manipulation
def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join([complement[base] for base in seq[::-1]])

print(reverse_complement('ATGC'))

def print_FASTA_headers(filename):
    """Print the headers of a FASTA file"""
    with open(filename) as file:
        for line in file:
            if line[0] == '>':
                print(line[1:-1])

print_FASTA_headers('sirt1_Mus musculus.fasta')

def extract_matching_sequences(filename, string):
    """From a FASTA file named filename, extract all sequences whose descriptions contain string"""
    sequences = []
    seq = ''
    with open(filename) as file:
        for line in file:
            if line[0] == '>':
                if seq: # not first time through
                    sequences.append(seq)
                seq = '' # next sequence detected
                includeflag = string in line # flag for later iterations
            else:
                if includeflag:
                    seq += line[:-1]
        if seq: # last sequence in file is included
            sequences.append(seq)
    return sequences

print(extract_matching_sequences('sirt1_Mus musculus.fasta', 'Mus musculus'))

def list_sequences_in_files(filelist):
    """For each file whose name is contained in filelist,
    list the description of each sequence it contains"""
    for filename in filelist:
        print(filename)
        with open(filename) as file:
            for line in file:
                if line[0] == '>':
                    print('\t', line[1:-1])

print(list_sequences_in_files(['sirt1_Mus musculus.fasta']))

DNA_bases = ['A', 'T', 'C', 'G']

def print_codon_table():
    """Print the DNA codon table in a nice, but simple, arrangement"""
    for base1 in DNA_bases: # horizontal section (or "group")
        for base3 in DNA_bases: # line (or "row")
            for base2 in DNA_bases: # vertical section (or "column")
                # the base2 loop is inside the base3 loop!
                print(base1+base2+base3,
                translate_DNA_codon(base1+base2+base3),
                end=' ')
            print()
        print()

print_codon_table()

