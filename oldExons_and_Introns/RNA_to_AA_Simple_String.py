#!/usr/bin/env python

'''
Rosalind Translating RNA into Protein
Sample Dataset
  AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA
Sample Output
  MAMAPRTEINSTRING
http://rosalind.info/problems/prot/
'''
# read in dataset, create output file
file = open('rosalind_prot.txt', 'r')
output = open('./rosalind_prot_answer.txt', 'w')
# it wasn't divisible by 3 until I used strip ... probably new line character!
RNA = file.readline().strip()
print(RNA)
RNA_length = len(RNA)
print("RNA length is: " + str(RNA_length))
# RNA can only be translated to protein in groups of 3
# this'll alert the user if there's a problem
if RNA_length % 3 == 0:
    print("It's disivble by 3!")
else:
    print("PROBLEM! Not divisible by 3!")

# good idea from StackOverflow!
# https://stackoverflow.com/questions/9475241/split-string-every-nth-character
'''
n = 3
line = RNA
print([line[i:i+n] for i in range(0, len(line), n)])
'''
# dictionary of RNA to Protein that I wrote out by hand :S
RNA_dict = {'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'UCU': 'S', 'UCC': 'S', \
'UCA': 'S', 'UCG': 'S', 'UAU': 'Y', 'UAC': 'Y', 'UAA': 'Stop', 'UAG': 'Stop', 'UGU': 'C', \
'UGC': 'C', 'UGA': 'Stop', 'UGG': 'W', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L', \
'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', \
'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M', \
'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', \
'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V', \
'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', \
'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

# this gets each trinucleotide all the way to the end!
three = 3
# initalize an empty string for the translated protein
protein_string = ''
# this for loop jumps up i by threes
for i in range(0, len(RNA), three):
    print(i)
    # it's i+3 INSTEAD of i+2 because the end of a python slice is NOT inclusive
    current_trinucleotide = RNA[i:i+3]
    # use the dictionary to translate RNA into protein
    current_amino_acid = RNA_dict[current_trinucleotide]
    print("Current trinucleotide is: ")
    print(current_trinucleotide)
    print("Current amino acid is: ")
    print(current_amino_acid)
    # break the loop if a Stop codon is encountered
    if current_amino_acid == 'Stop':
        print("STOP codon!")
        break
    protein_string += current_amino_acid
print("the protein string is: ")
print(protein_string)
output.write(protein_string)
