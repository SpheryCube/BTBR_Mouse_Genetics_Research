# Test 3

# Hypothesis: Decreased levels of autophagy gene expression in 
# pancreatic islet will result in severe diabetic traits. 
# Model: 
# Decreased Ccng2 genes -> increased insulin -> higher blood glucose

rm(list=ls())
directory <- "/home/daniel14/Il6_Mouse_Research/data/"
setwd(directory)
getwd()

# Load QTL library to do genome scans.
library(qtl)
library(ggplot2)
# Load data.
load(file = "BTBR.clean.data.Rdata")


phenotypes.rz$IL.6[phenotypes.rz$IL.6 > 0 & is.numeric(phenotypes.rz$IL.6)] <- 1  # Fix IL-6 data. Replace all 0 with n/a.
phenotypes.rz$IL.6[phenotypes.rz$IL.6 <= 0 & is.numeric(phenotypes.rz$IL.6)] <- -1  # Fix IL-6 data. Replace all 0 with n/a.



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


# Find the parkin gene (gene symbol Insr) in the annotation data.
grep("Insrr", annot$gene1)
grep("Insr", annot$gene1, value = TRUE)
grep("Il6", annot$gene1, value = TRUE)

# Find the ID number for the Park2 gene in the annotation data.
# Use this ID number to pull out expression data from pancreatic islet.
annot[grep("Insrr", annot$gene1),]
annot[grep("Insrr", annot$gene1), 1]

Insrr_gastroc <- gastroc.rz[, annot[grep("Insrr", annot$gene1), 1]]
Insrr_adipose <- adipose.rz[, annot[grep("Insrr", annot$gene1), 1]]

Il6_adipose <- adipose.rz[, annot[grep("Il6$", annot$gene1), 1]]
Il6_gastroc <- gastroc.rz[, annot[grep("Il6$", annot$gene1), 1]]
Il6ra_adipose <- adipose.rz[, annot[grep("Il6ra", annot$gene1), 1]]
Il6st_adipose <- adipose.rz[, annot[grep("Il6st", annot$gene1), 1]]

#phenotypes.rz$IL.6[phenotypes.rz$IL.6 < 0 & is.numeric(phenotypes.rz$IL.6)] <- NA
f2g$pheno <- cbind(f2g$pheno[,c("MouseNum","Sex","pgm")], phenotypes.rz[,c("GLU.10wk","INS.10wk", "IL.6")], Il6_adipose, Il6ra_adipose)

# look at all pairwise scatterplots of clinical and expression traits
names(f2g$pheno)
par(mfrow = c(1,1))
pairs(f2g$pheno[,4:length(f2g$pheno)], upper.panel=panel.cor,diag.panel=panel.hist)




# Move the clinical and gene expression phenotypes in to the cross object.

#phenotypes.rz$IL.6[phenotypes.rz$IL.6 < 0 & is.numeric(phenotypes.rz$IL.6)] <- NA
f2g$pheno <- cbind(f2g$pheno[,c("MouseNum","Sex","pgm")], phenotypes.rz[,c("GLU.10wk","INS.10wk", "IL.6")], Insrr_gastroc, Insrr_adipose, Il6_adipose, Il6ra_adipose, Il6st_adipose, Il6_gastroc)

# look at all pairwise scatterplots of clinical and expression traits
names(f2g$pheno)
par(mfrow = c(1,1))
pairs(f2g$pheno[,4:length(f2g$pheno)], upper.panel=panel.cor,diag.panel=panel.hist)





# Pull out sex as a numeric variable so that it can be used as a covariate
# in genome scans.
sex <- as.numeric(f2g$pheno$Sex)


################################# QTLs ###########################################

# Calculate genotype probabilities before running genome scans.
f2g <- calc.genoprob(f2g, step = 1)

# Run genome scans with sex as a covariate. This will allow the average
# phenotype values to differ between the two sexes.
scan1 <- scanone(f2g,  pheno.col = c("Insrr_adipose", "Insrr_gastroc", "IL.6", "Il6_adipose", "Il6ra_adipose", "Il6st_adipose", "INS.10wk", "GLU.10wk", "Il6_gastroc"), method = "hk", addcovar = sex)

# Identify LOD significance thresholds for gene expression traits.
perm1 <-scanone(f2g, pheno.col = 7, addcovar = sex, method = "hk", n.perm = 100, perm.Xsp = TRUE)


summary(perm1)

length(scan1)

# View genome scan summaries in different formats.
summary(scan1)
summary(scan1, format = "tabByChr")
summary(scan1, format = "tabByChr", perms=perm1, alpha = 0.05)




#####################  Plotting all original QTLs   #############################  

par(mfrow=c(3,1))
plot(scan1, lodcolumn = 1, main = "Quantitative Trait Loci for Insrr Expression in Adipose Tissue", ylab = "Insrr Adipose Exp", ylim = c(0,7) )
add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed",lwd=2,col="orange")
add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed",lwd=2,col="purple")

plot(scan1, lodcolumn = 2, main = "Quantitative Trait Loci for Insrr Expression in Gastrocnemius Tissue", ylab = "Insrr Muscle Exp", ylim = c(0,7))
add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed",lwd=2,col="orange")
add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed",lwd=2,col="purple")

plot(scan1, lodcolumn = 3, main = "Quantitative Trait Loci for Clinical IL-6 Plasma Levels", ylab = "IL-6 Levels", ylim = c(0,7))
add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed",lwd=2,col="orange")
add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed",lwd=2,col="purple")

par(mfrow=c(3,1))
plot(scan1, lodcolumn = 4, main = "Quantitative Trait Loci for Il6 Expression in Adipose Tissue", ylab = "IL-6 Adipose Exp", ylim = c(0,7))
add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed",lwd=2,col="orange")
add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed",lwd=2,col="purple")

plot(scan1, lodcolumn = 5, main = "Quantitative Trait Loci for Il6ra Expression in Adipose Tissue", ylab = "Ilra Adipose Exp", ylim = c(0, 7))
add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed",lwd=2,col="orange")
add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed",lwd=2,col="purple")

plot(scan1, lodcolumn = 6, main = "Quantitative Trait Loci for Il6st Expression in Adipose Tissue", ylab = "Il6st Adipose Exp", ylim = c(0,7))
add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed",lwd=2,col="orange")
add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed",lwd=2,col="purple")

par(mfrow=c(3,1))
plot(scan1, lodcolumn = 7, main = "Quantitative Trait Loci for Insulin at 10 weeks", ylab = "Insulin Levels", ylim = c(0,7))
add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed",lwd=2,col="orange")
add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed",lwd=2,col="purple")

plot(scan1, lodcolumn = 8, main = "Quantitative Trait Loci for Glucose at 10 weeks", ylab = "Glucose Levels", ylim = c(0,7))
add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed",lwd=2,col="orange")
add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed",lwd=2,col="purple")

plot(scan1, lodcolumn = 9, main = "Quantitative Trait Loci for Il6 Expression in Gastrocnemius Tissue", ylab = "Il6 Muscle Exp", ylim = c(0, 7))
add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed",lwd=2,col="orange")
add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed",lwd=2,col="purple")


par(mfrow=c(5, 2))
plot(scan1, lodcolumn = 1, main = "Quantitative Trait Loci for Insrr Expression in Adipose Tissue", ylab = "Insrr Adipose Exp", ylim = c(0,7) )
add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed",lwd=2,col="orange")
add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed",lwd=2,col="purple")

plot(scan1, lodcolumn = 2, main = "Quantitative Trait Loci for Insrr Expression in Gastrocnemius Tissue", ylab = "Insrr Muscle Exp", ylim = c(0,7))
add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed",lwd=2,col="orange")
add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed",lwd=2,col="purple")

plot(scan1, lodcolumn = 3, main = "Quantitative Trait Loci for Clinical IL-6 Plasma Levels", ylab = "IL-6 Levels", ylim = c(0,7))
add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed",lwd=2,col="orange")
add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed",lwd=2,col="purple")

plot(scan1, lodcolumn = 4, main = "Quantitative Trait Loci for Il6 Expression in Adipose Tissue", ylab = "IL-6 Adipose Exp", ylim = c(0,7))
add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed",lwd=2,col="orange")
add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed",lwd=2,col="purple")

plot(scan1, lodcolumn = 5, main = "Quantitative Trait Loci for Il6ra Expression in Adipose Tissue", ylab = "Ilra Adipose Exp", ylim = c(0, 7))
add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed",lwd=2,col="orange")
add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed",lwd=2,col="purple")

plot(scan1, lodcolumn = 6, main = "Quantitative Trait Loci for Il6st Expression in Adipose Tissue", ylab = "Il6st Adipose Exp", ylim = c(0,7))
add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed",lwd=2,col="orange")
add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed",lwd=2,col="purple")

plot(scan1, lodcolumn = 7, main = "Quantitative Trait Loci for Insulin at 10 weeks", ylab = "Insulin Levels", ylim = c(0,7))
add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed",lwd=2,col="orange")
add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed",lwd=2,col="purple")

plot(scan1, lodcolumn = 8, main = "Quantitative Trait Loci for Glucose at 10 weeks", ylab = "Glucose Levels", ylim = c(0,7))
add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed",lwd=2,col="orange")
add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed",lwd=2,col="purple")

plot(scan1, lodcolumn = 9, main = "Quantitative Trait Loci for Il6 Expression in Gastrocnemius Tissue", ylab = "Il6 Muscle Exp", ylim = c(0, 7))
add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed",lwd=2,col="orange")
add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed",lwd=2,col="purple")





par(mfrow=c(3,1))
plot(scan1, lodcolumn = 7, main = "Insulin Scan for Clinical Insulin Plasma Levels at 10 Weeks", ylab = "Insulin Levels", ylim = c(0, 7))
add.threshold(scan1, perms=perm1, alpha=0.05, lty="dashed", lwd=2, col="orange")
add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed", lwd=2, col="purple")

plot(scan_Il6, main = "Insulin Scan with IL-6 levels as Covariate", ylab = "Insulin")
add.threshold(scan1, perms=perm1, alpha=0.05, lty="dashed", lwd=2, col="orange")
add.threshold(scan1,perms=perm1, alpha=0.10, lty="dashed", lwd=2, col="purple")

plot(scan_Il6_adipose, main = "Insulin Scan with Il6 adipose Expression as Covariate",ylab = "Insulin")
add.threshold(scan1,perms=perm1, alpha=0.05,lty="dashed", lwd=2, col="orange")
add.threshold(scan1,perms=perm1, alpha=0.10,lty="dashed", lwd=2, col="purple")


# Insrr_adipose and Insrr_gastroc expression share the same peak on chr2.
# Their confidence intervals are the same
# What on chr9 drives insrr expression?
# We can download a list of all genes within the largest confidence
# interval (between ci.low and ci.high) from 51.3 cM to 65.5cM.
# Locate the markers nearest 51.3 cM and 65.5 cM. 
#Note that these two confidences intervals are within the confidence interval for Il6.

find.marker(cross = f2g, chr = 9, pos = 51.3)
find.marker(cross = f2g, chr = 9, pos = 65.5)

#Mine: rs1348035 (95860805 bp) and rs6335220 (113901404 bp)
chr9_genes <- scan(file = "/home/daniel14/CompBioProjects/BTBRxB6/data/gene_9_CI_for_Ins_Exp.txt",  what = "character",   skip = 1)


#IL.6 and INS.10wk share a 0.05 significance lod peak in chromosome 5: 16.3 cm to 57.1 cm
find.marker(cross = f2g, chr = 5, pos = 16.3)
find.marker(cross = f2g, chr = 5, pos = 65.3)
#Mine: rs13478154 (27114401 bp) and rs13478458 (110340025 bp)





chr5_genes <- scan(file = "/home/daniel14/Il6_Mouse_Research/data/ch5_genes(large).txt",  what = "character",   skip = 1)



# The confidence interval for insulin (50.6 cm to 57.1 cm) was nested 
# inside an outer confidence interval for the Il6 gene expression in 
# adipose tissue (42.7 cm to 65.3 cm).

find.marker(cross = f2g, chr = 5, pos = 42.7)
find.marker(cross = f2g, chr = 5, pos = 65.3)
#Mine: rs13478311 (67637841 bp) and rs13478458 (110340025 bp)

chr5_genes_small <- scan(file = "/home/daniel14/Il6_Mouse_Research/data/ch5_genes(small).txt",  what = "character",   skip = 1)



################ Biggest Interval #######################33
find.marker(cross = f2g, chr = 5, pos = 15.4)
find.marker(cross = f2g, chr = 5, pos = 71.3)
# rs13478154 (27114401 ) to rs13478536 ( 132255439)

chr5_genes_big <- scan(file = "/home/daniel14/Il6_Mouse_Research/data/ch5_genes_big.txt",  what = "character",   skip = 1)

# To download gene symbols, click on BioMart, choose database
# Ensembl Genes, choose dataset Mus musculus genes. 
# Click Filters. Open up Region, select chromosome 5, then enter
# the two base pair positions for gene start and end.
# Click Attributes. Open up Gene, uncheck Ensembl Gene ID and Ensembl
# Transcript ID. Open up External. Check MGI symbol.
# Click Results button at upper left. Export a file as TSV.
# The file should be named mart_export.txt. Change the name
# to something descriptive, like chr2_genes.txt.   
# Read the chromosome 2 gene list into R.



##########################################################
### I. conditional genome scans
##########################################################

# Looking for a gene expression trait that "blocks" the chr 5 QTL peak

# scan insulin conditional on Il6 gene expression trait in adipose
scan_Il6 <- scanone(f2g, pheno.col="INS.10wk", addcovar=f2g$pheno$IL.6, method="hk")
scan_Il6_adipose <- scanone(f2g, pheno.col="INS.10wk",addcovar=f2g$pheno$Il6_adipose, method="hk")
scan_Il6_gastroc <- scanone(f2g, pheno.col="INS.10wk",addcovar=f2g$pheno$Il6_gastroc, method="hk")

#quartz()
par(mfrow=c(3,1))

plot(scan1, lodcolumn = 3, main = "Quantitative Trait Loci for Clinical IL-6 Plasma Levels", ylab = "IL-6 Levels", ylim = c(0, 7))
add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed",lwd=2,col="orange")
add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed",lwd=2,col="purple")

plot(scan1, lodcolumn = 4, main = "Quantitative Trait Loci for Il6 Expression in Adipose Tissue", ylab = "Il6 Adipose Expression", ylim = c(0, 7))
add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed",lwd=2,col="orange")
add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed",lwd=2,col="purple")

plot(scan1, lodcolumn = 7, main = "Quantitative Trait Loci for Clinical Insulin Plasma Levels at 10 Weeks", ylab = "Insulin Levels", ylim = c(0, 7))
add.threshold(scan1, perms=perm1, alpha=0.05, lty="dashed", lwd=2, col="orange")
add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed", lwd=2, col="purple")




summary(scan1)

par(mfrow=c(3,1))

plot(scan1, lodcolumn = 4, main = "Quantitative Trait Loci for Il6 adipose expression")
add.threshold(scan1, perms=perm1, alpha=0.05, lty="dashed", lwd=2, col="orange")
add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed", lwd=2, col="purple")

plot(scan_Il6, main = "Insulin Scan with IL-6 levels as Covariate", ylab = "Insulin")
add.threshold(scan1, perms=perm1, alpha=0.05, lty="dashed", lwd=2, col="orange")
add.threshold(scan1,perms=perm1, alpha=0.10, lty="dashed", lwd=2, col="purple")

plot(scan_Il6_adipose, main = "Insulin Scan with Il6 adipose Expression as Covariate",ylab = "Insulin")
add.threshold(scan1,perms=perm1, alpha=0.05,lty="dashed", lwd=2, col="orange")
add.threshold(scan1,perms=perm1, alpha=0.10,lty="dashed", lwd=2, col="purple")




par(mfrow=c(4,1))
plot(scan1, lodcolumn = 7, main = "Quantitative Trait Loci for Clinical Insulin Plasma Levels at 10 Weeks", ylab = "Insulin Levels", ylim = c(0, 7))
add.threshold(scan1, perms=perm1, alpha=0.05, lty="dashed", lwd=2, col="orange")
add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed", lwd=2, col="purple")

plot(scan_Il6, main = "Insulin Scan with IL-6 levels as Covariate", ylab = "Insulin Levels", ylim = c(0, 7))
add.threshold(scan1, perms=perm1, alpha=0.05, lty="dashed", lwd=2, col="orange")
add.threshold(scan1,perms=perm1, alpha=0.10, lty="dashed", lwd=2, col="purple")

plot(scan_Il6_adipose, main = "Insulin Scan with Il6 Adipose Expression as Covariate",ylab = "Insulin Levels", ylim = c(0, 7))
add.threshold(scan1,perms=perm1, alpha=0.05,lty="dashed", lwd=2, col="orange")
add.threshold(scan1,perms=perm1, alpha=0.10,lty="dashed", lwd=2, col="purple")

plot(scan_Il6_gastroc, main = "Insulin Scan with Il6 Gastrocnemius Expression as Covariate",ylab = "Insulin Levels", ylim = c(0, 7))
add.threshold(scan1,perms=perm1, alpha=0.05,lty="dashed", lwd=2, col="orange")
add.threshold(scan1,perms=perm1, alpha=0.10,lty="dashed", lwd=2, col="purple")

summary(scan_Il6)
summary(scan_Il6_adipose)

# We found that when insulin was scanned conditionally on Il6 expression, the QTL peak on chromosome 5 dropped 
# from 4.75 to 2.73, lowering it below the 0.10 significance threshold.



##########################################################################################################    
######## Now let's find some genes in the CI that we can use for conditional scanning             ########
##########################################################################################################        
# Loop through all genes in the chromosome 5 interval, adding each
# in as covariate, running the scan for insulin, IL.6 and Il6_adipose

# Add chr5 gene expression traits into cross object.
f2g$pheno <- cbind(f2g$pheno,  adipose.rz[,match(chr5_genes_big, annot$gene1, nomatch = 0)])
names(f2g$pheno)

# Column names for gene expression traits are the microarray probe IDs.
# Replace probe IDs with gene symbols.
names(f2g$pheno)[13:ncol(f2g$pheno)] <- annot$gene1[match(names(f2g$pheno)[13:ncol(f2g$pheno)],annot$a_gene_id)]
names(f2g$pheno)

# scan insulin, Il.6, and Il6_adipose conditional on chr5 gene expression traits.
# Scan only for chromosome 5.

scan_cond <- scanone(f2g,  pheno.col=c("INS.10wk", "Il6_adipose"), addcovar=f2g$pheno$Mfsd7a, method="hk") #Mfsda is hte 13th element

# note how we used "cbind" to concatenate the results. This is where we add in all the genes in our confidence interval as covariates.
# "cbind" calls the specialized function "cbind.scanone"
# Do rest of conditional scans. The scan above was just to get it started and to show how to do the first one manually.
# This part takes awhile.


for(i in 14:ncol(f2g$pheno)){ #Start at the 14th element
  scan_cond <- cbind(scan_cond, scanone(f2g, pheno.col=c("INS.10wk", "Il6_adipose"), addcovar=f2g$pheno[,i], method="hk"))
}
summary(scan_cond)

# Re-name the conditional scan columns with the gene symbol. Every 3 columns
# will reference a gene scanned conditionally on insulin, IL.6, and Il6_adipose in
# that order.

head(names(scan_cond))
dim(scan_cond)
for (i in 13:(ncol(f2g$pheno))) {
  names(scan_cond)[(2*(i-12)+1):(2*(i-12)+2)] <- names(f2g$pheno)[i]
}
head(names(scan_cond))
names(scan_cond)
summary(names(scan_cond))



par(mfrow=c(6,1), mar=c(3,4,1,4) + 0.1)
for( i in 1:(ncol(scan_cond))){
  plot(scan_cond, lodcolumn=i, main=names(scan_cond)[i+2], chr =5)
  add.threshold(scan1, perms=perm1, alpha=0.05, lty="dashed", lwd=2, col="orange")
  add.threshold(scan1, perms=perm1, alpha=0.10, lty="dashed", lwd=2, col="purple")
}




my_scan_func <- function(x){
  par(mfrow=c(2, 2))
  
  
  plot(scan1, chr =5, lodcolumn = 7, main = "QTL for Clinical Insulin Plasma Levels at 10 Weeks", ylab = "Insulin Levels", ylim = c(0, 7))
  add.threshold(scan1, perms=perm1, alpha=0.05, lty="dashed", lwd=2, col="orange")
  add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed", lwd=2, col="purple")
  
  plot(scan1, chr =5,  lodcolumn = 4, main = "Quantitative Trait Loci for Il6 Expression in Adipose Tissue", ylab = "Il6 Adipose Expression", ylim = c(0, 7))
  add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed",lwd=2,col="orange")
  add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed",lwd=2,col="purple")
  
  
  #Find lod column number of Chst1 insulin, Chst1, il6, and Chst1 il6 adipose exp in scan_cond
  
  plot(scan_cond, chr = 5, lodcolumn = 2*(x-1)+1, main = paste("Insulin Scan with",  names(scan_cond[2*x-1]), "Adipose Expression as Covariate"), ylab= "Insulin Levels", ylim = c(0, 7)) #insulin
  add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed", lwd=2, col="orange")
  add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed", lwd=2, col="purple")
  
  plot(scan_cond, chr = 5, lodcolumn = 2*(x-1)+2, main = paste("Il6 Adipose Exp Scan with",  names(scan_cond[2*x]), "Adipose Expression as Covariate"), ylab = "Il6 Adipose Exp", ylim = c(0, 7)) #il6 adipose
  add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed", lwd=2, col="orange")
  add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed", lwd=2, col="purple")
}




# Which genes drop the chromosome 5 peak below the
# significance threshold of 3.67 for insulin,, Il6 Adipose Expression, AND Clinical IL-6 levels?

summary(subset(scan_cond, chr = 5),  thresholds = c(9.1, 7.1, 6.3, 6.3,3.3))

x = 1
lod_threshold = 3.67
candidate_genes <- vector()
while (2+2*x < length(summary(subset(scan_cond, chr = 5))))
{
  scores = vector()
  for (lod_score in summary(subset(scan_cond, chr = 5, threshold = c(9.1, 7.1, 6.3, 6.3, 3.3)))[(2*x+1):(2*x+2)]) #We want to skip the first two So that's why x starts at 1.
  {
    scores <- c(scores, lod_score)
    print(scores)
  }
  if (length(scores[scores < lod_threshold]) == 2)
  {
    print(names(scan_cond[2*x]))
    cat("Lod column numbers in scan_cond object: ", 2*(x-1)+1, " through, ", 2*(x-1)+2)
    print("")
    print("-----------------------")
    candidate_genes <- c(candidate_genes, names(scan_cond[2*x]))
    print(x)
    my_scan_func(x)
  }
  x = x + 1
}
print(candidate_genes)

# "Scfd2"  "Pdgfra" "Srp72"  "Gtf3c2" "Thap6"  "Cmklr1" "Ccng2"  "Cxcl5" 


summary(subset(scan_cond, chr = 5))[,which(names(scan_cond) %in% c("Scfd2", "Pdgfra", "Srp72", "Gtf3c3", "Thap6", "Cmklr1", "Ccng2", "Cxcl5"))]




#######################################################################################
####################### More Scatterplot Matrices #####################################
#######################################################################################

#Define some variables


Scfd2_adipose <- adipose.rz[, annot[grep("Scfd2", annot$gene1), 1]]
Pdgfra_adipose <- adipose.rz[, annot[grep("Pdgfra", annot$gene1), 1]]
Srp72_adipose <- adipose.rz[, annot[grep("Srp72", annot$gene1), 1]]
#Gtf3c3_adipose <- adipose.rz[, annot[grep("Gtf3c3", annot$gene1), 1]]
Thap6_adipose <- adipose.rz[, annot[grep("Thap6", annot$gene1), 1]]
Cmklr1_adipose <- adipose.rz[, annot[grep("Cmklr1", annot$gene1), 1]]
Ccng2_adipose <- adipose.rz[, annot[grep("Ccng2", annot$gene1), 1]]
Cxcl5_adipose <- adipose.rz[, annot[grep("Cxcl5$", annot$gene1), 1]]


f2g$pheno <- cbind(f2g$pheno[,c("MouseNum","Sex","pgm")], phenotypes.rz[,c("INS.10wk", "IL.6")], Il6_adipose, Scfd2_adipose, Pdgfra_adipose, Srp72_adipose, Thap6_adipose, Cmklr1_adipose, Ccng2_adipose, Cxcl5_adipose)
# look at all pairwise scatterplots of clinical and expression traits
names(f2g$pheno)
par(mfrow=c(1,1))
pairs(f2g$pheno[,4:length(f2g$pheno)], upper.panel=panel.cor,diag.panel=panel.hist)
# Ccng2 has a 0.21 positive correlation with IL-6 and a 0.12 correlation with Il6 gene expression in adipose tissue.

#Check out Pdgfra and Cxcl5




scan2 <- scanone(f2g,  pheno.col = c("IL.6", "INS.10wk", "Il6_adipose", "Scfd2_adipose", "Pdgfra_adipose", "Srp72_adipose", "Thap6_adipose","Cmklr1_adipose","Ccng2_adipose","Cxcl5_adipose"), method = "hk", addcovar = sex)

perm2 <-scanone(f2g, pheno.col = 7, addcovar = sex, method = "hk", n.perm = 100, perm.Xsp = TRUE)
summary(perm2)
length(scan2)

# View genome scan summaries in different formats.
summary(scan2)
summary(scan2, format = "tabByChr")
summary(scan2, format = "tabByChr", perms=perm2, alpha = 0.05)


plot(scan2, chr = 5, lodcolumn = c(1, 2, 3), col = c("red", "green", "blue"))

#Shared peak on chromosome 6 between insulin (42.8 cm), Cmklr1 (52.8 cm) and Cxcl5 (46.0)
#Chromosome 5 peak. Notice Ccng2_adipose shares a peak with IL6 (51.9). It is strongly associated with IL.6 (0.21 correlation)


f2g$pheno <- cbind(f2g$pheno[,c("MouseNum","Sex","pgm")], phenotypes.rz[,c("INS.10wk", "IL.6")], Il6_adipose, Ccng2_adipose, Pdgfra_adipose, Cxcl5_adipose)
par(mfrow=c(1,1))
pairs(f2g$pheno[,4:length(f2g$pheno)], upper.panel=panel.cor,diag.panel=panel.hist)

####################################################################
#############           Effect Plots and BIC analysis ##############
####################################################################

########################################################################
##################            Ccng2                   ##################            
########################################################################



print("##################            Ccng2     ##########################")
f2g$pheno <- transform(f2g$pheno, Q5 = as.factor(f2g$geno[[5]]$data[,find.marker(f2g, 5, 51.9)])) #The 51.9 cM location was found from scan2 above.
levels(f2g$pheno$Q5) <- c("B", "H", "R")

print("##################            Ccng2 - Insulin     #########################")
triple.fit(X = f2g$pheno$Ccng2_adipose, Y = f2g$pheno$INS.10wk, Q = f2g$pheno$Q5)
#Inconclusive

print("##################            Ccng2 - IL6     #############################")
triple.fit(X = f2g$pheno$Ccng2_adipose, Y = f2g$pheno$IL.6, Q = f2g$pheno$Q5)
#Causal model.   Q -> X -> Y

print("##################  Ccng2 - Il6 adipose expression     ##########")
triple.fit(X = f2g$pheno$Ccng2_adipose, Y = f2g$pheno$Il6_adipose, Q = f2g$pheno$Q5)
#Inconclusive.

#Effect plots.
f2g$pheno <- transform(f2g$pheno, Ccng2_adipose)
par(mfrow=c(3,1))
qplot(Ccng2_adipose, INS.10wk, color=Q5, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q5), method="lm", se=FALSE)
qplot(Ccng2_adipose, IL.6, color=Q5, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q5), method="lm", se=FALSE)
qplot(Ccng2_adipose, Il6_adipose, color=Q5, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q5), method="lm", se=FALSE)


##########################################################
# We found that the causal model linking Q5, IL-6 and Ccng2 adipose expression to be significant,
# but now we have to check to make sure the conditions for drawing such a conclusion are passed.

# check 4 conditions for Ccng2 gene expression as a mediator of Q5 effect on IL6
# Normally this is done before doing BIC modeling, but we are going to do it in reverse.
# If the conditions are met, we have a case for a pathway.

#davidakenny.net/cm/mediate.htm
####  i) IL-6 is linked to Q5                    # X -> Y
anova(lm(f2g$pheno$IL.6 ~ Sex + Q5, data = f2g$pheno))
# significant

####  ii) Ccng2 gene expression is linked to Q5     # X -> M
anova(lm(Ccng2_adipose ~ Sex + Q5, data = f2g$pheno))
# significant

####  iii) IL-6 not linked after accounting for Q5   
anova(lm(f2g$pheno$IL.6 ~ Sex + Ccng2_adipose + Q5, data = f2g$pheno))
# not significant * .

####  iv) Ccng2 gene expression is still linked after 
# accounting for IL-6
anova(lm(Ccng2_adipose ~ Sex + f2g$pheno$IL.6 + Q5, data = f2g$pheno))
# significant ***

# all 4 conditions for a mediator are satisfied
# so we have a case for a causal pathway.

##########################################################



#  III. establish mediator using model selection with BIC scoring

# note that missing values will mess up BIC analysis
apply(is.na(f2g$pheno), 2, sum)
# easiest to use the triple.fit function that removes missing data
# and fits the three models

# compute BIC scores for the causal triplet models
#   using Ccng2_adipose
with(f2g$pheno,triple.fit(Ccng2_adipose, IL.6, Q5))
#  causal model has lowest score
#  suggests that Ccng2_adipose is a strong mediator of IL.6



##########################################################


# Continuing from BIC analysis on IL6 clinical levels. 
# Goal 1: Find relationship between Ccng2, IL6, and Insulin receptors in order to make the causal chain between Q5, Ccng2 and IL6 longer.
# Goal 2: Identify relationship between Ccng2 and Pparg:  https://www.ncbi.nlm.nih.gov/pubmed/20844002/

grep("Ppar", annot$gene1, value = TRUE)
grep("Pparg$", annot$gene1, value = TRUE)



Insr_adipose <- adipose.rz[, annot[grep("Insr$", annot$gene1), 1]]
Insr_islet <- islet.rz[, annot[grep("Insr$", annot$gene1), 1]]
Pparg_adipose <- adipose.rz[, annot[grep("Pparg$", annot$gene1), 1]]
Pparg_islet <- islet.rz[, annot[grep("Pparg$", annot$gene1), 1]]


f2g$pheno <- cbind(f2g$pheno[,c("MouseNum","Sex","pgm")], phenotypes.rz[,c("INS.10wk", "IL.6")], Il6_adipose, Ccng2_adipose, Insr_adipose, Insr_islet)

pairs(f2g$pheno[,4:length(f2g$pheno)], upper.panel=panel.cor,diag.panel=panel.hist)



scan3 <- scanone(f2g,  pheno.col = c("IL.6", "INS.10wk", "Il6_adipose", "Ccng2_adipose", "Insr_adipose", "Insr_islet", "Pparg_adipose","Pparg_islet"), method = "hk", addcovar = sex)
# Identify LOD significance thresholds for gene expression traits.
perm3 <-scanone(f2g, pheno.col = 7, addcovar = sex, method = "hk", n.perm = 100, perm.Xsp = TRUE)

summary(perm3)
length(scan3)

# View genome scan summaries in different formats.
summary(scan3)
summary(scan3, format = "tabByChr")
summary(scan3, format = "tabByChr", perms=perm3, alpha = 0.05)


# https://www.ncbi.nlm.nih.gov/pubmed/20844002/
# Ccng2 = Cyclin-G2
# Positive regulator of adipogenesis through ppar gamma  (PPARG) coactivation




