# Telomere Maintenance Mechanisms
When DNA is getting copied, small end pieces aren't able to be fully copied. This is known as the End Replication Problem and is a result of Okazaki fragments incompletely covering the end of the chromosome. This is an issue for cell types that need to divide a lot, like Hematopoietic stem cells (HSC). HSCs have plenty of division to do so they can keep up with all the differentiation to make blood cells. This is why they express telomerase, which is a reverse transcriptase that adds repetitive telomeric DNA (TTAGGG)n to the ends of the chromosomes (Allsopp 2001). Did you notice that I'm citing a paper from the Weissman lab? ;)  
 
![telomeres-shorten-with-age](/Assets/telomeres-shorten-with-age.jpg "telomeres-shorten-with-age")

(Finkel 2007)

Cells that divide a lot will senesce in the absence of an active telomere maintenance mechanism (Shay 2012). That's why Telomerase is great for stem cells and other cell types that need to divide a lot. HOWEVER, it's also part of how the majority of cancers become immortal (Cesare 2010). This is why several companies are working on telomere-based anti-cancer therapies:

* Telomerase enzyme inhibitor: 
	* GRN163L (Geron)
* Telomerase active immunotherapy:	
	* GRNVAC1 (Geron)
	* GV1001 (Pharmexa)
	* P540-548 (Gemvax)
	* Vx01 (Vaxon Biotech)
	* TLI (Cosmo Bioscience) 

(Shay 2012, HARLEY 2008)

One problem with telomerase inhibitors is that they have been reported to cause problems with blood stem cells (Hu 2017). Another problem with the telomerase inhibition approach is that approximately 10-15% of cancers use the Alterantive Lengthening of Telomeres (ALT) to extend telomeres, some cancers won't be treated with these anti-telomerase therapies. But, it's way worse of a problem than that! Tumors have been reported to use both ALT and TEL simultaneously (Gocha 2013) AND in vitro inhibition of telomerase selects for ALT activity (Sahin 2012). 

![TEL_ALT_Reversible](/Assets/TEL_ALT_Reversible.jpg "TEL_ALT_Reversible")

(Shay 2012)

On an interesting note, there is a high frequency of cancers with a Mesenchymal Stem Cell origin using the ALT pathway (77% malignant fibrous histocytomas, 47-66% of osteosarcomas (Lafferty-Whyte 2009), 21.4% of liposarcomas (Venturini APB ALT LIPOSARCOMA 2008). This may be a result of the telomerase promoter being repressed with chromatin compaction (Atkinson 2005). Telomerase inhibition MIGHT not be alone in potentially causing problems for stem cells. There is some evidence to suggest that stem cells use ALT (Kalmbach 2014, Huang 2014). To make things worse, some cancers appear to extensively divide without a telomere maintenance mechanism (Dagg 2017).  

![stem_cell_ALT.jpg](/Assets/stem_cell_ALT.jpg "stem_cell_ALT.jpg")

(Kalmbach 2014)

# ATRX Exon Deletion is Common in ALT
This project can be found in the Human_ATRX_ALT folder. ATRX gene mutations are found in a range of cancers. 10-15% of cancers are estimated to use ALT. ALT involves homologous recombination-based telomere elongation. Inactivating mutations in either ATRX or DAXX are found in many cancers. Depletion of ATRX seems insufficient to trigger ALT, but it does seem to play a key role in the ALT pathway. The absence of ATRX might lead to the failure of stalled replication forks to get resolved. The required fork restart would require homologus recombination and could jumpstart the ALT pathway (Clynes 2013). ALT involves a template-based lengthening of telomeres with homologous recombination. The genetic and epigenetic changes are not full understood. Lovejoy 2012 reported that ATRX gene mutations are a common feature of ALT. Specifically 19/22 ALT+ cell lines had an issue with the expression of ATRX or DAXX (Lovejoy 2012). See the Lovejoy 2012 supplementary information for the Excel table of Exon deletions in ALT cell lines. 
![ATRX_Prevents_Fork_Collapse](/Assets/ATRX_Prevents_Fork_Collapse.jpg "ATRX_Prevents_Fork_Collapse")

(Clynes 2013)

#### Getting ATRX DNA
Searching Ensembl for human ATRX yielded ATRX-201 and ATRX-202. I picked ATRX-201 cause it has 35 exons (which matches the Lovejoy 2012 paper). It was Ensembly ENST00000373344.10. Ensembl refseq switch to NCBI Reference Sequence yielded NM_000489.5 for the gene. I saved it as NM_000489.5_homo_sapiens_ATRX_Gene.fasta.

#### Removing ATRX Exons 2-29 
See the Lovejoy 2012 supplementary Excel table for a list of commonly missing ATRX Exons. I decided to play with the U2OS variant because that is a cell line that I used to grow :) U2OS is missing ATRX exons 2-29. NCBI says exon 2 is [236:348] and exon 29 is 6542..6719 https://www.ncbi.nlm.nih.gov/nuccore/NM_000489. I used R to remove those exons.

```r
library(seqinr)
WT_hATRX_Gene <- read.fasta("NM_000489.5_homo_sapiens_ATRX_Gene.fasta")
WT_hATRX_Gene_Nucleotides <- WT_hATRX_Gene[[1]]
length(WT_hATRX_Gene_Nucleotides)
typeof(WT_hATRX_Gene_Nucleotides)
WT_hATRX_Gene_Nucleotides
U2OS_hATRX_Gene_Nucleotide_FIRST <- WT_hATRX_Gene_Nucleotides[1:235]
U2OS_hATRX_Gene_Nucleotide_SECOND <- WT_hATRX_Gene_Nucleotides[6720:length(WT_hATRX_Gene_Nucleotides)]
U2OS_ATRX_Characters <- c(U2OS_hATRX_Gene_Nucleotide_FIRST, U2OS_hATRX_Gene_Nucleotide_SECOND)
U2OS_ATRX_DNAstring <- DNAString(paste(toupper(U2OS_ATRX_Characters), collapse = ""))
```

#### Sequence Alignment of WT ATRX to Mutant ATRX
The R msa package can't handle the full length of the ATRX gene, so I shortened it down to 400 nucleotides.
```{r}
# limit it to 400 ... that's more than enough to see exon absence
U2OS_ATRX_DNA_Short <- U2OS_ATRX_DNAstring[1:400]
WT_hATRX_DNA_Short <- toupper(WT_hATRX_Gene_Nucleotides[1:400])
write.fasta(sequences = U2OS_ATRX_DNA_Short, names = "U2OS_ATRX_DNA_Short", file.out = "U2OS_ATRX_DNA_Short.fasta", open = "w", nbchar = 70, as.string = FALSE)
write.fasta(sequences = WT_hATRX_DNA_Short, names = "WT_hATRX_DNA_Short", file.out = "WT_hATRX_DNA_Short.fasta", open = "w", nbchar = 70, as.string = FALSE)
library(msa)
# both_ATRX_Sequences <- read.fasta("WT_and_U2OS_hATRX.fasta")
both_ATRX_Sequences_SHORT <- "both_ATRX_Sequences_SHORT.fasta"
# typeof(both_ATRX_Sequences)
#both_ATRX_DNAStringSet <- readDNAStringSet(both_ATRX_Sequences)
#both_ATRX_Sequences_Alignment <- msa(both_ATRX_DNAStringSet)
both_ATRX_Sequences_SHORT_StringSet <- readDNAStringSet(both_ATRX_Sequences_SHORT)
both_ATRX_Sequences_Alignment_SHORT <- msa(both_ATRX_Sequences_SHORT_StringSet)
both_ATRX_Sequences_Alignment_SHORT
msaPrettyPrint(both_ATRX_Sequences_Alignment_SHORT, output="pdf", showNames="none",
showLogo="none", askForOverwrite=FALSE, verbose=TRUE)
#texi2pdf("both_ATRX_Sequences_Alignment.tex", clean=TRUE)
```
You can see that the sequences are the same until postion 236. That is where the Exon deletion for U2OS starts! 
![ATRX_Exon_Deletion_Alignment](/Assets/ATRX_Exon_Deletion_Alignment.jpg "ATRX_Exon_Deletion_Alignment")

# TERT Promoter Compaction is Found in ALT
ALT cells commonly have long, heterogenous telomere lengths (Kumakura 2005). Mouse embryonic stem cells deficient for DNMT have HUGE telomeres. Under normal conditions, mouse subtelomeres are heavily methylated, BUT that is not the case in mESC deficient for DNMT. The lack of DNMT increased the rate of telomeric sister chromatid exchanges (T-SCE), and ALT-associated Promyelocytic Nuclear Bodies (APBs). T-SCE and APBs are both common features of ALT activity. The authors concluded that the increased telomeric recombination MIGHT lead to telomere length changes, BUT they do not exclude the involvement of telomerase in the weirdly long telomeres that were seen (Gonzalo 2006). Luckily, I found these two other papers that go into more detail about TERT chromatin compaction in ALT!

Atkinson 2005 found that chromatin modifications of hTR and hTERT promoters were commonly found in ALT activity. Treatment of ALT+ cells with 5-AZC or Trichostatin A lead to chromatin remodeling of hTR and hTERT. This induced telomerase expression. Interestingly enough they found that mehtylated Lys20 Histon H4 was not associated with gene expression, BUT does seem to be ALT specific (Atkinson 2005). This might be a new marker of ALT activity! Acetylation of H3K9 and methylation of H3K4 is known to be associated with an open chromatin conformation. In Kumakura 2005, the authors found that ALT+ cells had H3K9 methylation and low levels of H3K4 methylation and H3K9+H3K14 acetylations. The ratio of H3K9 methylation / H3K4 methylation was different across ALT+ and TEL+ cell lines. They found that treating ALT+ cells with TSC or 5-AZC caused a reversion from complete to partial methylation of the CpG islands on the hTERT promoter. They switched an E6CL TEL+ line to TEL- and it was able to grow for well over 240 population doublings (Kumakura 2005). That's some ALT activity right there!
![H3K9_H3K4_Methylation_Ratio](/Assets/H3K9_H3K4_Methylation_Ratio.jpg "H3K9_H3K4_Methylation_Ratio")

(Kumakura 2005)

#### Getting The hTERT Sequence
The UniProt entry O14746 is for hTERT https://www.uniprot.org/uniprot/O14746. Following GeneID 7015 gets TERT telomerase reverse transcriptase for humans https://www.ncbi.nlm.nih.gov/gene/7015. Note that the reverse arrow on TERT indicates that the sequence is on the reverse strand (this will become important later). I downloaded the FASTA as "NC_000005.10_hChrom5_TERT_CpG_Start.fasta". A quick text search shows that the start codon is at position 59. Take care with this sequence cause it's the reverse complement of the actual sequence!

![hTERT_NCBI_Reverse_Strand](/Assets/hTERT_NCBI_Reverse_Strand.jpg "hTERT_NCBI_Reverse_Strand")

#### Obtaining hTERT WITH the CpG Island Region
Stay with me ... we're about to dig a bit into the literature! The hTERT sequence that I grabbed from NCBI DOES NOT contain the CpG island for hTERT. It doesn't even contain the normal promoter region for hTERT! Cong 1999 reports that the core hTERT promoter region is from -330 to +361 bp of the ATG start codon. HOWEVER, Kumakura 2005 found the hTERT CpG island to be from 654 bp upstream to 510 bp downstream of the ATG start codon, so this is the actual region that I need to grab! 

Grabbing the hTERT FASTA sequence from https://www.ncbi.nlm.nih.gov/gene/7015 INITIALLY is from: 1253167 to: 1295047. Checking Cong 1999 and The FASTA file, I can see that the hTERT start codon, AND a bit more of that region, of "ATGCCGCGCGCT" is at the end of the first FASTA line, which is 59 in from the left (59 is A of ATG). CpG is 654 bp upstream of the transcriptoin start site, SO going 595 back from current start site 1253167-595 = 1252572. 1252572 should be the start of the CpG island, RIGHT?!? WRONG!!! ... Why aren't I getting any more nucleotides before the current start read?!?!? OH!!! I'm looking at the reverse complement, lol! ;) 

It should be  1295047 + 595 = 1295642 YESSSSSSS, that's right :) Now I have more at the beginniig, so "atgccgcgcgctccccgct" is fully searchable! This is the new range:
https://www.ncbi.nlm.nih.gov/nuccore/NC_000005.10?report=fasta&from=1253167&to=1295642&strand=true. I saved the FASTA as NC_000005.10_hChrom5_TERT_CpG_Start.fasta. 

#### Analyzing the Alleged CpG Promoter Region
Dessain 2000 reported that the hTERT CpG island has a GC content of 74% and a CG:GC ratio of 0.87. Is that what I get for the same region?!? I wrote code in R to get the GC content and CG:GC ratio of the hTERT CpG promoter region that I identified. I didn't comment my code ... I am sorry. Note that I picked i=1164 cause Kumakura 2005 says that the hTERT CpG island to be from 654 bp upstream to 510 bp downstream of the ATG start codon, which is 654 + 510 = 1164 :) Here's R code that I didn't bother commenting :( 

```r
dna_file <- read.fasta("/media/david/Linux/Introns_Exons_and_Promoters/hTERT_CpG_Island/NC_000005.10_hChrom5_TERT_CpG_Start.fasta")
individual_characters <- dna_file[[1]]
i <- 1
CG_count <- 0
GC_count <- 0
g_count <- 0
c_count <- 0
for (letter in individual_characters) {
  #print(letter)
  i <- i + 1
  if (letter == "g") {
    g_count <- g_count + 1
  }
  if (letter == "c") {
    c_count <- c_count + 1
  }
  if ((letter == "c") & (last_letter == "g")) {
    GC_count <- GC_count + 1
  }
  if ((letter == "g") & (last_letter == "c")) {
    CG_count <- CG_count + 1
  }
  last_letter <- letter
  if (i >= 1164) {
    break
  }
}
print(100*(c_count+g_count)/1164)
print(GC_count)
print(CG_count)
print(CG_count/GC_count)
```
The lazily unlabeled output is:

[1] 76.03093
[1] 167
[1] 141
[1] 0.8443114

Cool-ness! My GC content is 76 % and the ratio of CpG/GpC is 0.84. Recall that Dessain 2000 reported that the hTERT CpG island has a GC content of 74% and a CG:GC ratio of 0.87. I'm 2 % off of the GC content and 0.03 off of the CpG/GpC. THAT'S PRETTY GOOD FOR REPLICATING DATA FROM A PAPER THAT IS ALMOST TWO DECADES OLD :D But, what if that was just random luck? I re-ran that same R code on the region AFTER the CpG island and I got wildly different data. Here's that code (yes, I should've used a function; I am lazy, lol):

```r
dna_file <- read.fasta("/media/david/Linux/Introns_Exons_and_Promoters/hTERT_CpG_Island/NC_000005.10_hChrom5_TERT_CpG_Start.fasta")
individual_characters <- dna_file[[1]]
#individual_characters[5]
i <- 1164
CG_count <- 0
GC_count <- 0
g_count <- 0
c_count <- 0
for (letter in individual_characters) {
  if (i <1164) {
    next
  }
  #print(letter)
  i <- i + 1
  if (letter == "g") {
    g_count <- g_count + 1
  }
  if (letter == "c") {
    c_count <- c_count + 1
  }
  if ((letter == "c") & (last_letter == "g")) {
    GC_count <- GC_count + 1
  }
  if ((letter == "g") & (last_letter == "c")) {
    CG_count <- CG_count + 1
  }
  last_letter <- letter
  
  
  if (i >= 42476) {
    break
  }
}
print(100*(c_count+g_count)/41312)
print(GC_count)
print(CG_count)
print(CG_count/GC_count)
```

The lazily unlabeled output is:

[1] 58.55926
[1] 3047
[1] 1654
[1] 0.542829

My GC content is 58.5 % and the ratio of CpG/GpC is 0.54. Recall that Dessain 2000 reported that the hTERT CpG island has a GC content of 74% and a CG:GC ratio of 0.87. THE REGION THAT IS NOT A CpG ISLAND IS 15.5% off of the GC content and 0.33 off of the CpG/GpC. I could dig into this with more statistical rigor, but I think you get the idea. I'M EXCITED!!! This was a really cool biological programming exercise!!! :D

 



# POT-1 Deficiency Creates ALT+ C. Elegans Strains
Telomeres cap linear chromosomes because DNA polymerase can't completely copy chromosomes. Telomerase adds telomeric repeats to the ends of linear chromosomes with reverse transcription. Most human cancers have long, heterogenous telomeres. Telomere shortening leads to senescence and potentially crisis. Cancer emerges as part of massive cell death and genomic rearrangements after crisis. 10-15% of cancers are estimated to use ALT (Cheng 2012). 

ALT can happen in Caenorhabditis elegans! Mammalian POT1 has homologs in C. elegans as pot-1 (CeOB2) and pot-2 (CeOB1). What's the deal with the reversing of 1 and 2? That's how it's reported in the paper ... it's odd. pot-1 mutant C. elegans have HUGE telomere lengths while pot-2 mutants have normal telomere lengths. The authors of Cheng 2012 created a variety of mutants in C. elegans.  The trt-1 C. elegans mutant has a deletion in telomerase reverse transcriptase. trt-1 & pot-2 absence led to ALT+ Caenorhabditis elegans with normal telomere lengths. trt-1 and pot-1 mutants were found to have long, heterogenous telomere lengths like those seen in human ALT. Here is the survival figure showing that C. elegans can survive in the absence of telomerase reverse transcriptase.

![Celegans_ALT_Generation_Survival](/Assets/Celegans_ALT_Generation_Survival.jpg "Celegans_ALT_Generation_Survival")

(Cheng 2012)

#### Multiple Sequence Alignment of pot-1 Genes
YES, pot-2 was the central point of the paper, but it won't be as fun to play with because it only has one isoform. I picked pot-1 cause there is a lot of cool stuff to play with. There were a lot of workup steps to get all of the sequences ... It would take a long while to review them. Essentially, I looked up the proteins on UniProt and then grabbed the DNA files from NCBI GenBank and WormBase. Check out the Celegans_POT1_ALT folder for the file names of everything. The file containing all the C. elegans genes is Celegans_POT1_genes.fasta. I used the R package "msa" for multiple sequence alignment with this code:

```r
library(msa)
Celegans_POT1_genes <- "/media/david/Linux/Introns_Exons_and_Promoters/Celegans_POT1_ALT/DNA/Celegans_POT1_genes.fasta"
Celegans_POT1_genes_DNA <- readDNAStringSet(Celegans_POT1_genes)
Celegans_POT1_gene_alignment <- msa(Celegans_POT1_genes_DNA)
msaPrettyPrint(Celegans_POT1_gene_alignment, output="pdf", showNames="none",
showLogo="none", askForOverwrite=FALSE, verbose=FALSE)
```

The aligned sequences aren't very pretty ... I decided not to include sequence labels cause it shortened the available nucleotide space for each new line. Here's part of the output for you to get the idea of the work:

![Celegans_POT1_gene_alignment](/Assets/Celegans_POT1_gene_alignment.jpg "Celegans_POT1_gene_alignment")

#### Multiple Sequence Alignment of pot-1 Proteins
I grabbed all the C elegans pot-1 isoform sequences from UniProt. You can check them out in Celegans_POT1_ALT/Protein. The file containing all of the sequences is Celegans_POT1_Proteins.fasta. I aligned all of the proteins with code that is similar to the DNA alignment code:

```r
Celegans_POT1_proteins <- "/media/david/Linux/Introns_Exons_and_Promoters/Celegans_POT1_ALT/Protein/Celegans_POT1_Proteins.fasta"
Celegans_POT1_proteins_AA <- readAAStringSet(Celegans_POT1_proteins)
Celegans_POT1_protein_alignment <- msa(Celegans_POT1_proteins_AA)
Celegans_POT1_protein_alignment
msaPrettyPrint(Celegans_POT1_protein_alignment, output="pdf", showNames="none",
showLogo="none", askForOverwrite=FALSE, verbose=FALSE)
```

This alignment looks great! You can see all of the alignments between the different pot-1 isoforms :)

![Celegans_POT1_protein_alignment](/Assets/Celegans_POT1_protein_alignment.jpg "Celegans_POT1_protein_alignment")

#### Displaying pot-1 Open Reading Frames
I used code from the BioPython Tutorial to identify the C. elegans pot-1 gene Open Reading Frames AND to report the translated proteins! Compare this to the last section to see that the the DNA -> Protein translations all have the correct lengths! The code is from http://biopython.org/DIST/docs/tutorial/Tutorial.html in the section titled "20.1.13. Identifying open reading frames". You should check this code out! It can identify Open Reading Frames in the +/- strand AND in three different reading frames, SO IT DOES ALL 6 FRAMES!!! I re-used the code four times (instead of making a function, haha). Here's part of the code that I used:

```python
#!/usr/bin/env python

from Bio import SeqIO
# record = SeqIO.read("NC_005816.fna", "fasta")
file_0 = "NM_001361730.1_Caenorhabditis_elegans_pot-1_gene_homolog.fasta"
file_1 = "NM_001361731.1_Caenorhabditis_elegans_pot-1_gene_homolog.fasta"
file_2 = "NM_001361732.1_Caenorhabditis_elegans_pot-1_gene_homolog.fasta"
file_3 = "NM_066157.3_pot-1_NCBI_DNA_matchesP42001_Celegans.fasta"

print("")
print("NM_001361730.1_Caenorhabditis_elegans_pot-1_gene_homolog.fasta")
record = SeqIO.read(file_0, "fasta")
table = 11
min_pro_len = 150

for strand, nuc in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
    for frame in range(3):
        length = 3 * ((len(record)-frame) // 3) #Multiple of three
        for pro in nuc[frame:frame+length].translate(table).split("*"):
            if len(pro) >= min_pro_len:
                print("%s...%s - length %i, strand %i, frame %i" \
                % (pro[:30], pro[-3:], len(pro), strand, frame))
```

![Displaying_pot-1_Open_Reading_Frames](/Assets/Displaying_pot-1_Open_Reading_Frames.jpg "Displaying_pot-1_Open_Reading_Frames")

#### Discussing C. elegans pot-1 Alternative Splicing
WormBase has three isoforms for pot-1 in C. elegans https://wormbase.org/species/c_elegans/gene/WBGene00015105#0-9g-3

![WormBase_pot-1_Celegans_Isoforms](/Assets/WormBase_pot-1_Celegans_Isoforms.jpg "WormBase_pot-1_Celegans_Isoforms")

Transcript B0280.10a.1 is 1216 nucleotides in length and codes for a 400 amino acid protien. The WormBase curators used RNA-seq data from Boeck 2016 to alter the original WormBase entry to include the published alternate intron splicing. You can see from the table that Exons 1-10 are part of this pot-1 isoform.

![B0280.10a.1_Celegans_pot-1_isoformA_1203NT_400AA](/Assets/B0280.10a.1_Celegans_pot-1_isoformA_1203NT_400AA.jpg "B0280.10a.1_Celegans_pot-1_isoformA_1203NT_400AA")

Transcript B0280.10b.1 is 462 nucleotides in length and it codes for a 153 amino acid protein. You can see from the table that Exons 1-4 are part of this pot-1 isoform. Boeck 2016 goes into more detail about the alternative intron that causes alternative splicing here. 

![B0280.10b.1_Celegans_pot-1_isoformB_462NT_153AA](/Assets/B0280.10b.1_Celegans_pot-1_isoformB_462NT_153AA.jpg "B0280.10b.1_Celegans_pot-1_isoformB_462NT_153AA")

Transcript B0280.10c.1 is 1140 nucleotides long and it codes for a 379 amino acid protein. You can see from the table that Exons 1-10 make it into this protein isoform.

![B0280.10c.1_Celegans_pot-1_isoformC_1140NT_379AA](/Assets/B0280.10c.1_Celegans_pot-1_isoformC_1140NT_379AA.jpg "B0280.10c.1_Celegans_pot-1_isoformC_1140NT_379AA")


# STN1 Mutation Triggers ALT in Yeast
ALT is a recombination-based telomere maintenace mechanism used by some human cancers to maintain cellular immortality. ALT cells typically have widly long, heterogenous telomeres. The exact molecular mechanism involved in this pathway are unknown. Iyer 2005 found that a STN1 gene mutation can initiate ALT in the yeast Kluyveromyces lactis. These ALT-like yeast cells experience a rapid telomere shortening when WT STN1 is re-introduced (Iyer 2005). There aren't any neat figures in this paper to talk about, so I'm going to try to replicate the protein multiple sequence alignment from Iyer 2005 figure 4. STN1 is aligned from K. lactis S. cerevisiae (Sc) and Candida glabrata (Cgl).

![Yeast_Protein_Alignment](/Assets/Yeast_Protein_Alignment.jpg "Yeast_Protein_Alignment")

#### Obtaining STN1 From 3 Yeast Organisms
The paper states that the S. cerevisiae (Sc) and Candida glabrata (Cgl) GenBank accession numbers are P_38960 and XP_448655. HOWEVER, it wasn't that easy to find the sequences cause those accession numbers are from back in 2005. NCBI says "The following term was not found in Nucleotide: P_38960." The XP_448655 is here: https://www.ncbi.nlm.nih.gov/protein/XP_448655. I had trouble finding the sequence for K. lactis they were talking about. Here's the crazy search term that I used on NCBI:

and it's the only result (other than whole chromosome chunks) w/ this crazy search term I made:
(((stn1) NOT "Pyrenophora tritici-repentis"[porgn:__txid45151] NOT "Fusarium fujikuroi"[porgn:__txid5127] NOT "[Candida] glabrata"[porgn:__txid5478] NOT "Hortaea werneckii"[porgn:__txid91943] NOT "Saccharomyces cerevisiae"[porgn:__txid4932]) NOT "Metarhizium robertsii"[porgn:__txid568076] NOT "Fusarium sp. FIESC_5 CS3069"[porgn:__txid1318460] NOT "Fusarium pseudograminearum CS3487"[porgn:__txid1318458] NOT "Fusarium pseudograminearum CS3427"[porgn:__txid1318457] NOT "Fusarium pseudograminearum CS3220"[porgn:__txid1318456] NOT "Fonsecaea multimorphosa"[porgn:__txid979981] NOT "Cladophialophora immunda"[porgn:__txid569365] NOT "Aspergillus nidulans FGSC A4"[porgn:__txid227321] NOT "Candida viswanathii"[porgn:__txid5486] NOT "Zygosaccharomyces bailii"[porgn:__txid4954] NOT "Metarhizium anisopliae"[porgn:__txid5530] NOT "Aspergillus flavus"[porgn:__txid5059] NOT "Talaromyces atroroseus"[porgn:__txid1441469] NOT "[Candida] auris"[porgn:__txid498019] NOT "Zygosaccharomyces rouxii"[porgn:__txid4956] NOT "[Candida] boidinii"[porgn:__txid5477] NOT "Komagataella phaffii"[porgn:__txid460519] NOT "Aspergillus fumigatus"[porgn:__txid746128] NOT "Candida albicans SC5314"[porgn:__txid237561] NOT "Yarrowia lipolytica"[porgn:__txid4952]) AND "Kluyveromyces lactis"[porgn:__txid28985] 

... I'm pretty sure that it's NCBI Reference Sequence: XM_452728.1, titled "Kluyveromyces lactis uncharacterized protein (KLLA0_C11825g), partial mRNA" BECAUSE that entry has this note "cerevisiae YDR082W STN1 Protein involved in telomere length regulation functions in telomere metabolism during late S phase." S. ceriviasa yeast was easier to find cause it's got stn1p in the title :) https://www.ncbi.nlm.nih.gov/nuccore/NM_001180390.1
                     
#### Attempt to Replicate Lyer 2005 Figure 4
I used the R msa library to align the three yeast organism STN1 proteins that were obtained above. I'm not too happy with this alignment :( The authors state that "The protein
sequences were aligned in ClustalW using default values." HOWEVER, the R msa package that I used was set to use ClustalW default values and it didn't replicate Figure 4 from Lyer 2005. The protein sequences seem to match (at least by eye). I'm not really sure where I went wrong here. I've only played with sequence alignment a little bit, so I'm probably doing something wrong. Anyway, here's the code:

```r
setwd("/media/david/Linux/Introns_Exons_and_Promoters/Yeast_STN1_ALT/Proteins")

library(msa)
All_Yeast_STN1_Protein <- "All_Yeast_STN1_Protein.fasta"
All_Yeast_STN1_Protein <- readDNAStringSet(All_Yeast_STN1_Protein)
All_Yeast_STN1_Protein_alignment <- msa(All_Yeast_STN1_Protein)
All_Yeast_STN1_Protein_alignment
msaPrettyPrint(All_Yeast_STN1_Protein_alignment, output="pdf", showNames="none",
showLogo="none", askForOverwrite=FALSE, verbose=FALSE)
```

![Aligning_Yeast_STN1_Proteins](/Assets/Aligning_Yeast_STN1_Proteins.jpg "Aligning_Yeast_STN1_Proteins")


#### Aligning Human, Yeast, and Frog STN1
Cohen 2002 reported that Xenopus laevis form extrachromosomal circular telomeric DNA. This is commonly associated with ALT! I briefly looked around for reptile ALT and this is the closest thing I could find. I'm not convinced that the activity reported by Cohen 2002 was actually ALT. This alignment isn't really related to anything. I just made it mostly for fun ... well, it's kinda connected to ALT and I know a herpetologist that might enjoy seeing frogs getting included in this repo ;) 

Cohen 2002 Formation of extrachromosomal circles from telomeric DNA in Xenopus laevis 
```r
Klactis_Xenopus_Human_STN1_Proteins <- "Klactis_Xenopus_Human_STN1_Proteins.fasta"
Klactis_Xenopus_Human_STN1_Proteins_AA <- readAAStringSet(Klactis_Xenopus_Human_STN1_Proteins)
Klactis_Xenopus_Human_STN1_Proteins_AA_alignment <- msa(Klactis_Xenopus_Human_STN1_Proteins_AA)
Klactis_Xenopus_Human_STN1_Proteins_AA_alignment
msaPrettyPrint(Klactis_Xenopus_Human_STN1_Proteins_AA_alignment, output="pdf", showNames="none",
showLogo="none", askForOverwrite=FALSE, verbose=FALSE)
```
![Klactis_Xenopus_Human_STN1_Proteins_AA_alignment](/Assets/Klactis_Xenopus_Human_STN1_Proteins_AA_alignment.jpg "Klactis_Xenopus_Human_STN1_Proteins_AA_alignment")

# Citations
* Cheng 2012 Caenorhabditis elegans POT-2 telomere protein represses a mode of alternative lengthening of telomeres with normal telomere lengths
* Cong 1999 The human telomerase catalytic subunit hTERT: organization of the gene and characterization of the promoter
* Lovejoy 2012 PLoS Genet Loss of ATRX, genome instability, altered DNA damage response hallmarks of ALT pathway
* Clynes 2013 Curr Opin Genet Dev ATRX and the replication of structured DNA
* Gonzalo 2006 DNA methyltransferases control telomere length and telomere recombination in mammalian cells.pdf
* Atkinson 2005 Lack of Telomerase Gene Expression in Alternative Lengthening of Telomere Cells Is Associated with Chromatin Remodeling of the hTR and hTERT Gene Promoters
* Kumakura 2005 Reversible Conversion of Immortal Human Cells from Telomerase-Positive to Telomerase-Negative Cells
* Cong 1999 The human telomerase catalytic subunit hTERT: organization of the gene and characterization of the promoter
* Dessain 2000 Methylation of the Human Telomerase Gene CpG Island
* Cheng 2012 Caenorhabditis elegans POT-2 telomere protein represses a mode of alternative lengthening of telomeres with normal telomere lengths
* Boeck 2016 The time resolved transcriptome of C. elegans
* Iyer 2005 A Mutation in the STN1 Gene Triggers an Alternative Lengthening of Telomere-Like Runaway Recombinational Telomere Elongation and Rapid Deletion in Yeast
* Cohen 2002 Formation of extrachromosomal circles from telomeric DNA in Xenopus laevis 
* Hu 2017 Imetelstat, a Telomerase Inhibitor, Is Capable of Depleting Myelofibrosis Hematopoietic Stem Cells and Progenitor Cells
* Weissman 2001 Telomere Shortening Accompanies Increased Cell Cycle Activity during Serial Transplantation of Hematopoietic Stem Cells
* Cesare 2010 Alternative lengthening of telomeres: models, mechanisms and implications
* Shay 2012 Cancer and Telomeres −− An ALTernative to Telomerase
* Harley 2008 Telomerase and cancer therapeutics
* Sahin 2012 Antitelomerase therapy provokes ALT and mitochondrial adaptive mechanisms in cancer
* Gocha 2013 Human Sarcomas Are Mosaic for Telomerase-Dependent and Telomerase-Independent Telomere Maintenance Mechanisms
* Kalmbach 2014 Telomere Length Reprogramming in Embryos and Stem Cells
* Huang 2014 Telomere regulation in pluripotent stem cells
* Dagg 2017 Extensive Proliferation of Human Cancer Cells with Ever-Shorter Telomeres

