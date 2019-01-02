#!/usr/bin/env python
from Bio import SeqIO

# this is the location of the human ATRX gene
gene_fasta = "/media/david/Linux/Introns_Exons_and_Promoters/Data/Homo_sapiens_ATRX_Gene.fasta"
# this is the location of the human ATRX protein
protein_fasta = "/media/david/Linux/Introns_Exons_and_Promoters/Data/Homo_sapiens_ATRX_Protein.fasta"

# making a record of the gene_fasta
gene_record = SeqIO.read(gene_fasta, "fasta")
# makign a copy of the gene_fasta in another .fasta file
SeqIO.write(gene_record, "hATRX_Copied_Gene.fasta", "fasta")

# making a record of the protein protein_fasta
protein_record = SeqIO.read(protein_fasta, "fasta")

# write protein
SeqIO.write(protein_record, "hATRX_Copied_Protein.fasta", "fasta")
