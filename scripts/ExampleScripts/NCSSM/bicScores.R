#################################################
## Author: Hannah Gahagan and Nikhil Milind    ##
## Date: January 12 2016                       ##
#################################################

#################################################
## ENVIRONMENTAL SETUP                         ##
#################################################

# Clear environmental variables
rm(list=ls())
setwd("/Users/gotwals/desktop/BTBR/")
# Load QTL library to do genome scans
# Import required libraries.
library(qtl)
library(ggplot2)
library(qtlnet)
# Load data.
load(file = "BTBR.clean.data.Rdata")

source("important_func.R")



#################################################
## FUNCTIONS                                   ##
#################################################

# Get gene expression data for a gene with a given name from a tissue
# gene.name: The common MGI gene symbol
# data.set: The gene expression dataset, such as adipose.rz
gene.exp <- function(gene.name, data.set) {
  
  #if (gene.data.exists(gene.name)) {
    
    return(data.set[,annot$a_gene_id[which(annot$gene_symbol==gene.name)]])
  #}
}

# Get clinical data for a certain parameter
# clin.name: The name of the clinical trait in the dataset
# data.set: Either phenotypes or phenotypes.rz
clinical <- function(clin.name, data.set) {
  
  if (clin.name %in% names(data.set)) {
    
    return(data.set[,clin.name])
  }
}

# Get the genotype of the SNPs at a certain position on a chromosome
# chr: The chromosome on which to look for the genotype data
# pos: The position at which to look for genotype data
genotype <- function(chr, pos) {
  
  if (chr > 0 && chr < 21 && pos >= 0 && pos <= 100) {
    
    return(f2g$geno[[chr]]$data[,find.marker(f2g, chr=chr, pos=pos)])
  }
}

# BIC Score model analysis
# X: gene expression data for a given gene
# Y: quantitative measurement of clinical phenotype
# Q: genotype at a given marker
triple.fit <- function(X, Y, Q) {
  
  # Remove any NA values from the data
  indx <- sort(unique(c(which(is.na(X)), which(is.na(Y)), which(is.na(Q)))))
  X <- X[-indx]
  Y <- Y[-indx]
  Q <- Q[-indx]
  print(paste("Removed ", length(indx), " rows with NA values from data.", sep=""))
  
  # Calculate BIC scores for models
  bic.independent <- BIC(lm(X~Q)) + BIC(lm(Y~Q)) # X<-Q->Y
  bic.reactive <- BIC(lm(X~Y)) + BIC(lm(Y~Q)) # Q->Y->X
  bic.causal <- BIC(lm(X~Q)) + BIC(lm(Y~X)) # Q->X->Y
  bic.complex <- BIC(lm(X~Q)) + BIC(lm(Y~Q+X))
  
  # Print out the scores from each model
  print("BIC Scores of each model")
  scores <- c(bic.independent, bic.reactive, bic.causal, bic.complex)
  names(scores) <- c("independent", "reactive", "causal", "complex")
  print(scores)
  
  # Make lowest BIC score 0 and linearize all other scores accordingly to calculate Delta values
  deltas <- scores - min(scores)
  
  # Estimate the strength of evidence for each model
  strengths <- exp(-0.5 * deltas) / sum(exp(-0.5 * deltas))
  
  # Print out the probabilities of each model being the likely explanation for the data
  print("Probability of each model explaining the data")
  print(strengths * 100)
  
  # Print out how many more times likely the best model is
  print("The factor by which the best model is better than the rest")
  print(max(strengths) / strengths)
  
  # References
  # [1] Burnham, K. P., and D. R. Anderson. 2002. Model selection and multimodel inference: a practical information-theoretic approach. 
  #     Second edition. Springer, New York, USA.
  # [2] Anderson, D. R. 2008. Model based inference in the life sciences: a primer on evidence.
  #     Springer, New York, USA.
}

