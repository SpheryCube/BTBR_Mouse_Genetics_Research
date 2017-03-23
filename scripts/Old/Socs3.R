rm(list=ls())
directory <- "/home/daniel14/Il6_Mouse_Research/data/"
setwd(directory)
getwd()

# Load QTL library to do genome scans.
library(qtl)
library(ggplot2)
# Load data.
load(file = "BTBR.clean.data.Rdata")




####
# fit causal models to a triplet with BIC scoring
# X is a transcript used here as first argument to make "apply" easy. Gene expression
# Y is a clincal trait
# Q is a genotype (factor or numeric)

triple.fit <- function(X,Y,Q){
  #remove any rows with missing values
  indx <- sort(unique(
    c(which(is.na(X)),which(is.na(Y)),which(is.na(Q)))
  ))
  X <- X[-indx]
  Y <- Y[-indx]
  Q <- Q[-indx]
  
  # fit models and compute scores
  b1 <- BIC(lm(X~Q)) + BIC(lm(Y~Q))	#independent X<-Q->Y
  b2 <- BIC(lm(X~Y)) + BIC(lm(Y~Q))	#reactive	 Q->Y->X
  b3 <- BIC(lm(X~Q)) + BIC(lm(Y~X))	#causal		 Q->X->Y
  b4 <- BIC(lm(X~Q)) + BIC(lm(Y~Q+X)) #complex
  scores <- c(b1,b2,b3,b4)
  names(scores) <- c("independent","reactive","causal","complex")
  scores
}


####
# some useful plotting functions
# see documentation for "pairs" function 
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, use="complete.obs")
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor*0.5, col=c("gray60", "black")[(abs(r)>0.65)+1])
}
#
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5))
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}

###################################################################
#############     End of useful plotting function    ##############
###################################################################


###################### Establish variables and perliminary scatterplots ######################3
# Looking into Socs3
grep("Socs3", annot$gene1, value = TRUE)
annot[grep("Socs3", annot$gene1), 1]

Socs3_adipose <- adipose.rz[, annot[grep("Socs3", annot$gene1), 1]]
Socs3_gastroc <- gastroc.rz[, annot[grep("Socs3", annot$gene1), 1]]

Insrr_gastroc <- gastroc.rz[, annot[grep("Insrr", annot$gene1), 1]]
Insrr_adipose <- adipose.rz[, annot[grep("Insrr", annot$gene1), 1]]

Il6_adipose <- adipose.rz[, annot[grep("Il6$", annot$gene1), 1]]
Il6_gastroc <- gastroc.rz[, annot[grep("Il6$", annot$gene1), 1]]
Il6ra_adipose <- adipose.rz[, annot[grep("Il6ra", annot$gene1), 1]]
Il6st_adipose <- adipose.rz[, annot[grep("Il6st", annot$gene1), 1]]


phenotypes.rz$IL.6[phenotypes.rz$IL.6 > 0 & is.numeric(phenotypes.rz$IL.6)] <- 1  # Fix IL-6 data. Replace all 0 with n/a.
phenotypes.rz$IL.6[phenotypes.rz$IL.6 <= 0 & is.numeric(phenotypes.rz$IL.6)] <- -1  # Fix IL-6 data. Replace all 0 with n/a.


f2g$pheno <- cbind(f2g$pheno[,c("MouseNum","Sex","pgm")], phenotypes.rz[,c("GLU.10wk","INS.10wk", "IL.6")], Insrr_gastroc, Insrr_adipose, Il6_adipose, Il6ra_adipose, Il6st_adipose, Il6_gastroc, Socs3_adipose, Socs3_gastroc)
names(f2g$pheno)
par(mfrow = c(1,1))
pairs(f2g$pheno[,4:length(f2g$pheno)], upper.panel=panel.cor,diag.panel=panel.hist)


################################# QTLs ###########################################

# Pull out sex as a numeric variable so that it can be used as a covariate
# in genome scans.
sex <- as.numeric(f2g$pheno$Sex)

# Calculate genotype probabilities before running genome scans.
f2g <- calc.genoprob(f2g, step = 1)

# Run genome scans with sex as a covariate. This will allow the average
# phenotype values to differ between the two sexes.
scan1 <- scanone(f2g,  pheno.col = c("IL.6", "GLU.10wk", "INS.10wk", "Insrr_gastroc", "Insrr_adipose", "Il6_adipose", "Il6ra_adipose", "Il6st_adipose", "Il6_gastroc", "Socs3_adipose", "Socs3_gastroc"), method = "hk", addcovar = sex)

# Identify LOD significance thresholds for gene expression traits.
perm1 <-scanone(f2g, pheno.col = 7, addcovar = sex, method = "hk", n.perm = 100, perm.Xsp = TRUE)


summary(perm1)

length(scan1)

# View genome scan summaries in different formats.
summary(scan1)
summary(scan1, format = "tabByChr")
summary(scan1, format = "tabByChr", perms=perm1, alpha = 0.05)






