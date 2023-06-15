library(cubatr)
# setwd("./tests")
sequences <- readDNAStringSet("sequences.fasta")
codonFreq <-  calculateCodonFreq(sequences)
aaFreq <- calculateAminoAcidFreq(sequences)


optimalCodons <- findOptimalCodonsBinomial("high.fasta", "low.fasta")$significantPositive
calculateFOP(codonFreq,optimalCodons)
