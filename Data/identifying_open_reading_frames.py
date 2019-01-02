#!/usr/bin/env python

# http://biopython.org/DIST/docs/tutorial/Tutorial.html
# 20.1.13. Identifying open reading frames
# https://biopython.readthedocs.io/en/latest/Tutorial/chapter_cookbook.html

from Bio import SeqIO
# record = SeqIO.read("NC_005816.fna", "fasta")
file = "Homo_sapiens_ATRX_Gene.fasta"

record = SeqIO.read(file, "fasta")
table = 11
min_pro_len = 100

for strand, nuc in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
    for frame in range(3):
        length = 3 * ((len(record)-frame) // 3) #Multiple of three
        for pro in nuc[frame:frame+length].translate(table).split("*"):
            if len(pro) >= min_pro_len:
                print("%s...%s - length %i, strand %i, frame %i" \
                % (pro[:30], pro[-3:], len(pro), strand, frame))

    for frame in range(3):
        length = 3 * ((len(record)-frame) // 3) #Multiple of three
        for pro in nuc[frame:frame+length].translate(table).split("*"):
            if len(pro) >= min_pro_len:
                print("%s...%s - length %i, strand %i, frame %i" \
                % (pro[:30], pro[-3:], len(pro), strand, frame))

