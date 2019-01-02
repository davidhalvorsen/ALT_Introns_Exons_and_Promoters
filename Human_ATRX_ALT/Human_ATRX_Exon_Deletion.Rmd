---
title: "Human ATRX Exon Deletion"
author: "Dave_Halvorsen"
date: "January 2, 2019"
output: html_document
---

```{r}
# who's in my personal space?
setwd("/media/david/Linux/Introns_Exons_and_Promoters/Human_ATRX_ALT")
getwd()
list.files()
```

```{r}
library(seqinr)
WT_hATRX_Gene <- read.fasta("NM_000489.5_homo_sapiens_ATRX_Gene.fasta")
WT_hATRX_Gene_Nucleotides <- WT_hATRX_Gene[[1]]
length(WT_hATRX_Gene_Nucleotides)
typeof(WT_hATRX_Gene_Nucleotides)
WT_hATRX_Gene_Nucleotides
```

```{r}
# U2OS is missing ATRX exons 2-29
# NCBI says exon 2 is [236:348]
# exon 29 is 6542..6719
# https://www.ncbi.nlm.nih.gov/nuccore/NM_000489

U2OS_hATRX_Gene_Nucleotide_FIRST <- WT_hATRX_Gene_Nucleotides[1:235]
U2OS_hATRX_Gene_Nucleotide_SECOND <- WT_hATRX_Gene_Nucleotides[6720:length(WT_hATRX_Gene_Nucleotides)]
# U2OS_hATRX_Gene_Nucleotide <- paste(U2OS_hATRX_Gene_Nucleotide_FIRST, U2OS_hATRX_Gene_Nucleotide_SECOND, sep=" ")
U2OS_ATRX_Characters <- c(U2OS_hATRX_Gene_Nucleotide_FIRST, U2OS_hATRX_Gene_Nucleotide_SECOND)
U2OS_ATRX_DNAstring <- DNAString(paste(toupper(U2OS_ATRX_Characters), collapse = ""))
#U2OS_hATRX_Gene_Nucleotide <- paste(toupper(concat_these), collapse = "")
#writeXStringSet(paste(toupper(concat_these), "./U2OS_hATRX_Gene_Nucleotide.fasta"))
#typeof(U2OS_hATRX_Gene_Nucleotide)
#U2OS_hATRX_Gene_Nucleotide


```

```{r}
# write fasta of U2OS ATRX
write.fasta(sequences = U2OS_ATRX_DNAstring, names = "U2OS hATRX Missing Exons 2-29", file.out = "U2OS_hATRX_Gene_Nucleotide.fasta", open = "w", nbchar = 70, as.string = FALSE)



```

```{r}
library(msa)
# both_ATRX_Sequences <- read.fasta("WT_and_U2OS_hATRX.fasta")
both_ATRX_Sequences <- "WT_and_U2OS_hATRX.fasta"
typeof(both_ATRX_Sequences)
#both_ATRX_DNAStringSet <- readDNAStringSet(both_ATRX_Sequences)
#both_ATRX_Sequences_Alignment <- msa(both_ATRX_DNAStringSet)
both_ATRX_Sequences_StringSet <- readDNAStringSet(both_ATRX_Sequences)
both_ATRX_Sequences_Alignment <- msa(both_ATRX_Sequences_StringSet)
both_ATRX_Sequences_Alignment
```

```{r}
print(both_ATRX_Sequences_Alignment, show="complete")

```

```{r}
# THIS IS TOO LONG
#install.packages("pdflatex")
#library(pdflatex)
#install.packages("texi2pdf")
#library(texi2pdf)
msaPrettyPrint(both_ATRX_Sequences_Alignment, output="pdf", showNames="none",
showLogo="none", askForOverwrite=FALSE, verbose=TRUE)
#texi2pdf("both_ATRX_Sequences_Alignment.tex", clean=TRUE)
```

```{r}
# limit it to 400 ... that's more than enough to see exon absence
U2OS_ATRX_DNA_Short <- U2OS_ATRX_DNAstring[1:400]
WT_hATRX_DNA_Short <- toupper(WT_hATRX_Gene_Nucleotides[1:400])
write.fasta(sequences = U2OS_ATRX_DNA_Short, names = "U2OS_ATRX_DNA_Short", file.out = "U2OS_ATRX_DNA_Short.fasta", open = "w", nbchar = 70, as.string = FALSE)
write.fasta(sequences = WT_hATRX_DNA_Short, names = "WT_hATRX_DNA_Short", file.out = "WT_hATRX_DNA_Short.fasta", open = "w", nbchar = 70, as.string = FALSE)
```


```{r}
library(msa)
# both_ATRX_Sequences <- read.fasta("WT_and_U2OS_hATRX.fasta")
both_ATRX_Sequences_SHORT <- "both_ATRX_Sequences_SHORT.fasta"
# typeof(both_ATRX_Sequences)
#both_ATRX_DNAStringSet <- readDNAStringSet(both_ATRX_Sequences)
#both_ATRX_Sequences_Alignment <- msa(both_ATRX_DNAStringSet)
both_ATRX_Sequences_SHORT_StringSet <- readDNAStringSet(both_ATRX_Sequences_SHORT)
both_ATRX_Sequences_Alignment_SHORT <- msa(both_ATRX_Sequences_SHORT_StringSet)
both_ATRX_Sequences_Alignment_SHORT
```

```{r}
print(both_ATRX_Sequences_Alignment_SHORT, show="complete")

```

```{r}
# THIS IS TOO LONG
#install.packages("pdflatex")
#library(pdflatex)
#install.packages("texi2pdf")
#library(texi2pdf)
msaPrettyPrint(both_ATRX_Sequences_Alignment_SHORT, output="pdf", showNames="none",
showLogo="none", askForOverwrite=FALSE, verbose=TRUE)
#texi2pdf("both_ATRX_Sequences_Alignment.tex", clean=TRUE)
```