#' calculateCodonFreq
#'
#' Create codonFreq to store codon frequency.
#' @param sequences DNAStringSet of \code{sequences}
#' @return codonFreq
#' @author Hebin Zhang
#' @export
calculateCodonFreq <- function(sequences) {
  # speed up possibility
  # if(class(sequences) == "data.frame"){
  #   return(sequences)
  # }

  # Create codonFreq to store codon frequency
  codonFreq <- trinucleotideFrequency(step = 3, with.labels = TRUE, sequences)
  codonFreq <- as.data.frame(codonFreq, stringsAsFactors = FALSE)
  ID <- names(sequences)
  codonFreq <- cbind(ID,codonFreq)
  return(codonFreq)
}



#' calculateAminoAcidFreq
#'
#' Create aaFreq to store codon frequency.
#' @param sequences DNAStringSet of \code{sequences}
#' @param codonTable select codon table
#' @return aaFreq
#' @author Hebin Zhang
#' @export
calculateAminoAcidFreq <- function(sequences, codonTable = "standard") {
  # There is aaFreq for amino acid frequency
  aminoAcids <- c("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His",
                  "Ile","Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp",
                  "Tyr", "Val", "Stop")
  selectedCodonTable = getCodonTable(codonTable)
  ID <- names(sequences)
  codonFreq <-  calculateCodonFreq(sequences)
  # Create an empty aaFreq to store data
  aaFreq <- data.frame(matrix(nrow = length(ID), ncol = length(aminoAcids),
                               dimnames = list(NULL, aminoAcids)))
  aaFreq <- cbind(ID, aaFreq)
  aaFreq[is.na(aaFreq)] <- 0

  # The codon frequency corresponding to that amino acid is found and then sum.
  for (col in colnames(codonFreq[, -1])) {
    targetAA <- selectedCodonTable$aminoAcid[selectedCodonTable$codon == col]
    targetColumn <- aaFreq[[targetAA]]
    aaFreq[, targetAA] <- codonFreq[[col]] + targetColumn
  }

  return(aaFreq)
}




