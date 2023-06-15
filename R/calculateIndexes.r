#' calculatecalculateRSCU
#'
#' calculate sequences RSCU.
#' @param sequences DNAStringSet of \code{sequences}
#' @param codonTable select codon table \code{codonTable}
#' @return rscuSheet
#' @author Hebin Zhang
#' @export
calculateRSCU <- function(sequences, codonTable = "standard") {
  # RSCU：https://pubmed.ncbi.nlm.nih.gov/3526280/
  selectedCodonTable <-  getCodonTable(codonTable)

  # Create rscuSheet to store RSCU
  rscuSheet <-  calculateCodonFreq(sequences)
  for (col in colnames(codonFreq[, -1])) {
    codonCount <- as.numeric(codonFreq[[col]])
    targetAA <-
      selectedCodonTable$aminoAcid[selectedCodonTable$codon == col]
    aaCount <- as.numeric(aaFreq[[targetAA]])
    aaMultiplicity <- sum(selectedCodonTable$aminoAcid == targetAA)
    rscu <-
      sum(selectedCodonTable$aminoAcid == targetAA) * codonCount / aaCount
    rscuSheet[, col] <- ifelse(aaCount == 0, 0, rscu)
  }
  stopCodons <-
    selectedCodonTable$codon[selectedCodonTable$aminoAcid == "Stop"]

  rscuSheet <-
    subset(rscuSheet, select = !names(rscuSheet) %in% stopCodons)

  return(rscuSheet)
}



#' calculateRelAdapt
#'
#' calculate relative adaptability.
#' @param highExpressedGene highExpressedGene DNAStringSet of \code{highExpressedGene}
#' @param codonTable select codon table \code{codonTable}
#' @return relAdapt
#' @author Hebin Zhang
#' @export
calculateRelAdapt <-
  function(highExpressedGene, codonTable = "standard") {
    highExpressedGene <- readDNAStringSet(highExpressedGene)
    highExprcodonFreq <-
      colSums(calculateCodonFreq(highExpressedGene)[, -1])

    aminoAcids <- c("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly",
                    "His","Ile","Leu", "Lys", "Met", "Phe", "Pro", "Ser",
                    "Thr", "Trp", "Tyr", "Val", "Stop")
    selectedCodonTable <-  getCodonTable(codonTable)

    relAdapt <- vector("numeric")
    for (aa in aminoAcids) {
      synCodons <-
        selectedCodonTable$codon[selectedCodonTable$aminoAcid == aa]
      maxSynCodonFreq <- max(highExprcodonFreq[synCodons])

      for (col in synCodons) {
        relAdapt[col] <-
          as.numeric(highExprcodonFreq[[col]]) / maxSynCodonFreq
      }
    }

    return(relAdapt)
  }

#' calculateCAI
#'
#' calculate cai(codon adaptation index).
#' @param sequences DNAStringSet of \code{sequences}
#' @param highExpressedGene highExpressedGene DNAStringSet of \code{highExpressedGene}
#' @param codonTable select codon table \code{codonTable}
#' @return cai
#' @author Hebin Zhang
#' @export
# cai2: xia Hua
calculateCAI <-
  function(sequences,highExpressedGene,codonTable = "standard") {
    if (class(sequences) == "DNAStringSet"){
      codonFreq <-  calculateCodonFreq(sequences)
    }else{
      codonFreq <- sequences
    }

    if (is.character(highExpressedGene)){
      RelAdapt <- calculateRelAdapt(highExpressedGene, codonTable)
    }else{
      RelAdapt <- highExpressedGene
    }

    codonFreq <- calculateCodonFreq(sequences)

    # creat cai dataFrame
    ID <- names(sequences)
    cai <- data.frame(matrix(
      nrow = length(ID),
      ncol = 2,
      dimnames = list(NULL, c("CAI", "CAI2"))
    ))
    cai <- cbind(ID, cai)
    cai[is.na(cai)] <- 0

    for (col in colnames(codonFreq[,-1])) {
      cai[, "CAI2"] <- cai[, "CAI2"] + codonFreq[[col]] * RelAdapt[col]
      cai[, "CAI"] <- cai[, "CAI"] + codonFreq[[col]] * log(RelAdapt[col])
    }

    cai[, "CAI2"] <- cai[, "CAI2"] / rowSums(codonFreq[, -1])
    cai[, "CAI"] <- exp(cai[, "CAI"] / rowSums(codonFreq[, -1]))
    return(cai)
  }

#' findOptimalCodonsBinomial
#'
#' Finding Optimal codons using binomial regression.
#' @param highFile Enter the fasta file of the highly expressed gene \code{highFile}
#' @param lowFile Enter the fasta file of the lowly expressed gene \code{lowFile}
#' @return binomialResults
#' @author Hebin Zhang
#' @export
findOptimalCodonsBinomial <- function(highFile, lowFile) {
  high <- readDNAStringSet(highFile)
  low <- readDNAStringSet(lowFile)
  codonFreqHigh <- calculateCodonFreq(high)
  codonFreqLow <- calculateCodonFreq(low)

  # Remove the ID column and transpose
  transposedHigh <- t(codonFreqHigh[, -1])
  transposedLow <- t(codonFreqLow[, -1])

  # Create the response variable (label): high expression as 1, low expression as 0
  expression <- factor(c(rep(1, ncol(transposedHigh)), rep(0, ncol(transposedLow))))

  # Combine the data
  combinedData <- data.frame(t(cbind(transposedHigh, transposedLow)))
  combinedData$Expression <- expression

  # Fit a binomial regression model
  model <- glm(Expression ~ ., data = combinedData, family = binomial)

  # Get the model summary information
  modelSummary <- summary(model)

  # Extract coefficients, standard errors, and p-values
  coefficients <- coef(modelSummary)
  standardErrors <- coef(modelSummary)[, "Std. Error"]
  pValues <- coef(modelSummary)[, "Pr(>|z|)"]

  # Create empty vectors to store significantly positive and negative parameters
  significantPositive <- c()
  significantNegative <- c()

  # Categorize based on p-values and coefficient signs
  for (i in 2:(length(coefficients)/4)) {
    if (pValues[i] < 0.05) {
      if (coefficients[i] > 0) {
        significantPositive <- append(significantPositive, rownames(coefficients)[i])
      } else if (coefficients[i] < 0) {
        significantNegative <- append(significantNegative, rownames(coefficients)[i])
      }
    }
  }

  return(list(
    coefficients = coefficients,
    standardErrors = standardErrors,
    pValues = pValues,
    significantPositive = significantPositive,
    significantNegative = significantNegative
  ))
}

#' calculateFOP
#'
#' calculate FOP(frequency of optimal codons)
#' @param codonFreq  codon frequency \code{codonFreq}
#' @param optimalCodons Enter optimal Codons \code{lowFile}
#' @return fop
#' @author Hebin Zhang
#' @export
calculateFOP <- function(codonFreq,optimalCodons){

  fop <- data.frame(matrix(
    nrow = length(rownames(codonFreq)),
    ncol = 1,
    dimnames = list(NULL, c("FOP"))
  ))
  fop <- cbind(ID=codonFreq[,1], fop)
  fop[is.na(fop)] <- 0
  fop[, "FOP"] <- rowSums(codonFreq[, optimalCodons])/rowSums(codonFreq[, -1])
  return(fop)
}
