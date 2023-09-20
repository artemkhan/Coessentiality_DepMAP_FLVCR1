# Load packages
library(tidyverse)

# Read in the file
df <- read.csv("CRISPR_(DepMap_Public_23Q2+Score,_Chronos).csv") 
# Select the column with FLVCR1
flvcr1 <- df["FLVCR1"]
Sys.time()
# Initialize the table
acc <- c()

Sys.time()

# Iterate through each gene and compute pearson correlation
for (i in 2:ncol(df)) {
  
  # Select the data for gene of interest
  current_gene <- df[,i]
  # Compute Pearson correlation and p-value for the association between the gene of interest and FLVCR1
  parameters <- cor.test(unlist(flvcr1), unlist(current_gene) )
  pval <- parameters$p.value
  pc <- parameters$estimate
  name <- colnames(df)[i]
  to_acc <- c(name, pc, pval)
  # Add to the accumulator variable
  acc <- rbind(acc, to_acc)
  print(i) 
  }

Sys.time()

acc <- as_tibble(acc)
# Compute correlation
acc$cor <- as.numeric(acc$cor)
# Compute rank
acc$rank <- rank(acc$cor)
plot(acc$rank, acc$cor)
# Remove FLVCR1 (check that it should be 1)
acc <- acc[acc$V1 != "FLVCR1",]

# Write the output
write.csv(acc, "20230627_FLVCR1_Coessentiality.csv")
