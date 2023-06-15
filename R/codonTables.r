#' getCodonTable
#'
#' return a codon table
#' @param specie select a specie \code{specie}
#' @return codonTable
#' @author Hebin Zhang
#' @export
getCodonTable <- function(specie) {
  codonTable <- switch(specie,
                        standard = data.frame(
                          codon = c("TTT", "TTC", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG",
                                    "ATT", "ATC", "ATA", "ATG", "GTT", "GTC", "GTA", "GTG",
                                    "TCT", "TCC", "TCA", "TCG", "CCT", "CCC", "CCA", "CCG",
                                    "ACT", "ACC", "ACA", "ACG", "GCT", "GCC", "GCA", "GCG",
                                    "TAT", "TAC", "TAA", "TAG", "CAT", "CAC", "CAA", "CAG",
                                    "AAT", "AAC", "AAA", "AAG", "GAT", "GAC", "GAA", "GAG",
                                    "TGT", "TGC", "TGA", "TGG", "CGT", "CGC", "CGA", "CGG",
                                    "AGT", "AGC", "AGA", "AGG", "GGT", "GGC", "GGA", "GGG"),
                          aminoAcid = c("Phe", "Phe", "Leu", "Leu", "Leu", "Leu", "Leu", "Leu",
                                        "Ile", "Ile", "Ile", "Met", "Val", "Val", "Val", "Val",
                                        "Ser", "Ser", "Ser", "Ser", "Pro", "Pro", "Pro", "Pro",
                                        "Thr", "Thr", "Thr", "Thr", "Ala", "Ala", "Ala", "Ala",
                                        "Tyr", "Tyr", "Stop", "Stop", "His", "His", "Gln", "Gln",
                                        "Asn", "Asn", "Lys", "Lys", "Asp", "Asp", "Glu", "Glu",
                                        "Cys", "Cys", "Stop", "Trp", "Arg", "Arg", "Arg", "Arg",
                                        "Ser", "Ser", "Arg", "Arg", "Gly", "Gly", "Gly", "Gly")
                        ),
                        custom = list(
                          # codon_table <- list(...)
                        ),
                        stop("Invalid specie type. Please choose 'standard' or others.")
  )

  return(codonTable)
}
