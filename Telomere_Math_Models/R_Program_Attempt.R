#!/usr/bin/Rscript

chromosome_lengths <- c(246, 202, 193, 184, 173, 160, 146, 125, 110, 103, 100, 93, 85, 81, 62, 67, 48, 52)
telomere_lengths <- c(5681, 4987, 5018, 4589, 4302, 4127, 3922, 3708, 4045, 3624, 3460, 3109, 3077, 3007, 2750, 2913, 2735, 2806)
telomere_lengths <- telomere_lengths + 2000

chromosome_X_length <- 14.93*164 + 3920.25
#chromosome_X_length
chromosome_Y_length <- 14.93*59 + 3920.25
#chromosome_Y_length
chromosome_1_length <- telomere_lengths[1]
chromosome_2_length <- telomere_lengths[1]
chromosome_9_length <- telomere_lengths[8]
chromosome_10_length <- telomere_lengths[8]
chromosome_11_length <-telomere_lengths[8]
chromosome_12_length <- telomere_lengths[8]

chrom_1 <- chromosome_1_length
chrom_2 <- chromosome_2_length
chrom_3 <- telomere_lengths[2]
chrom_4 <- telomere_lengths[3]
chrom_5 <- telomere_lengths[4]
chrom_6 <- telomere_lengths[5]
chrom_7 <- telomere_lengths[6]
chrom_8 <- telomere_lengths[7]
chrom_9 <- chromosome_9_length
chrom_10 <- chromosome_10_length
chrom_11 <- chromosome_11_length
chrom_12 <- chromosome_12_length
chrom_13 <- telomere_lengths[9]
chrom_14 <- telomere_lengths[10]
chrom_15 <- telomere_lengths[11]
chrom_16 <- telomere_lengths[12]
chrom_17 <- telomere_lengths[13]
chrom_18 <- telomere_lengths[14]
chrom_19 <- telomere_lengths[15]
chrom_20 <- telomere_lengths[16]
chrom_21 <- telomere_lengths[17]
chrom_22 <- telomere_lengths[18]
chrom_X <- chromosome_X_length
chrom_Y <- chromosome_Y_length

telomerase_chromosome_lengths <- c(chrom_1, chrom_2, chrom_3, chrom_4, chrom_5, chrom_6, chrom_7, chrom_8, chrom_9, chrom_10, chrom_11, chrom_12,chrom_13,chrom_14,chrom_15,chrom_16,chrom_17,chrom_18,chrom_19,chrom_20, chrom_21, chrom_22, chrom_X, chrom_Y)
telomerase_chromosome_names <- c((1:22), "X", "Y")

above_zero <- TRUE
division_number <- 0
telomerase_chromosome_lengths_play <- telomerase_chromosome_lengths
while(above_zero) {
  telomerase_chromosome_lengths_play <- telomerase_chromosome_lengths_play - 50
  print(telomerase_chromosome_lengths_play)
  print(division_number)
  plot(telomerase_chromosome_lengths_play, names.arg=telomerase_chromosome_names, las=2, main ="TEL+ Telomere Lengths for Chromosomes", xlab="Chromosome", ylab="Telomere Lenghts (bp)")
  division_number <- division_number + 1
  if(sum(telomerase_chromosome_lengths_play < 0) > 0) {
    above_zero <- FALSE
  }
}
plot(telomerase_chromosome_lengths_play)
print("poop")
telomerase_chromosome_lengths_play
barplot(telomerase_chromosome_lengths_play)
