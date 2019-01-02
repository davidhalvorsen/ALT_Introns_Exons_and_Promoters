#This program takes in a sequence of DNA that contains introns and exons. It
#returns a sequence of coding bases upper case and non-coding bases as lowercase.
my_dna = "ATCGATCGATCGATCGACTGACTAGTCATAGCTATGCATGTAGCTACTCGATCGATCGATCGATCGATCGATCGATCGATCGATCATGCTATCATCGATCGATATCGATGCATCGACTACTAT"
#The instructions state that exon 1 starts at CHARACTER 1 and ends at CHARACTER 63
#Slicing starts at 0, so I need to subtract the character positions by 1 to use with slicing.
exon1 = my_dna[0:63]
#The instructions state that the intron 1 starts at the 64th character and runs to the 90th character.
#Subtracting 1 from the chracter positions results in the slicing positions.
intron1 = my_dna[63:90]
#The instructions state that exon 2 starts at CHARAcTER 91 and goes to the end of the sequence.
#The same -1 slicing rules apply AND recall that, to go to end of sequence, leave the end slice blank.
exon2 = my_dna[90: ]
#Here the upper and lower methods are used to make exons upper case and introns lower case.
upper_exon1 = exon1.upper()
lower_intron1 = intron1.lower()
upper_exon2 = exon2.upper()
#The case-transitioned string after concatenation is the following:
upperExon_lowerIntron = upper_exon1 + lower_intron1 + upper_exon2

print("The provided DNA sequence is: ")
print(my_dna)
print("Exon 1 is: ")
print(exon1)
print("Intron 1 is: ")
print(intron1)
print("Exon 2 is: ")
print(exon2)
print("The coding regions of the DNA sequence were made upper case. The non-coding regions were made lower case.")
print("The string with the case transitions is as follows: ")
#Here's the final print call for the completed string:
print(upperExon_lowerIntron)
