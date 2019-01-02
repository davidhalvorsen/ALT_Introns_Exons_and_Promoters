#!/usr/bin/env python

# this prints the human ATRX gene
file = "/media/david/Linux/Introns_Exons_and_Promoters/Data/Homo_sapiens_ATRX_Gene.fasta"
from Bio import SeqIO
for seq_record in SeqIO.parse(file, "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
