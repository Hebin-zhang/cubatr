library(cubatr)
library(ggplot2)

# setwd("./tests")
sequences <- readDNAStringSet("sequences.fasta")
codonFreq <-  calculateCodonFreq(sequences)
# invertebrate yeast vertebrate moldProtozoan
aaFreq <- calculateAminoAcidFreq(sequences, codonTable = "yeast")




# ----RSCU----
RSCU <- calculateRSCU(sequences)
rscuMatrix <- as.matrix(RSCU[,-1])

rscu_long <- reshape2::melt(RSCU, id.vars = "ID")

ggplot(rscu_long, aes(x = variable, y = ID, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 6))

# ----CAI----
calculateRelAdapt("high.fasta")

cai <- calculateCAI(sequences, "high.fasta")
ggplot(cai, aes(x = ID)) +
  geom_bar(aes(y = CAI), position = "dodge", stat = "identity", width = 0.7) +
  geom_bar(aes(y = CAI2), position = "dodge", stat = "identity", width = 0.7,
           alpha = 0.5) +
  coord_flip()

# ----FOP----

high <- readDNAStringSet("high.fasta")
low <- readDNAStringSet("low.fasta")
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
modelSummary

prediction <- predict(model, newdata = codonFreq, type = "response")
prediction

# -------
optimalCodons <- findOptimalCodonsBinomial("high.fasta", "low.fasta")$significantPositive
calculateFOP(codonFreq,optimalCodons)
