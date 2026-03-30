""" Importing sequences from Entrez directly """
from Bio import Entrez, SeqIO
from keys import Email, NCBI

Entrez.email = Email
Entrez.api_key = NCBI


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

