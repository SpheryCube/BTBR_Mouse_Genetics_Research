# Hypothesis: Decreased levels of autophagy gene expression in 
# pancreatic islet will result in severe diabetic traits. 
# Model: 
# decreased Atg genes -> increased insulin -> higher blood glucose

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

phenotypes.rz$IL.6[phenotypes.rz$IL.6 > 0 & is.numeric(phenotypes.rz$IL.6)] <- 1  # Fix IL-6 data. Replace all 0 with n/a.
phenotypes.rz$IL.6[phenotypes.rz$IL.6 <= 0 & is.numeric(phenotypes.rz$IL.6)] <- -1  # Fix IL-6 data. Replace all 0 with n/a.

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


#Plot genome scans for each phenotype.
# par(mfrow=c(3, 1))
# for(i in 1:(length(scan1)-2)){
#   plot(scan1, lodcolumn = i)
#   add.threshold(scan1,
#                 perms=perm1, alpha=0.05,
#                 lty="dashed",lwd=2,col="orange")
#   add.threshold(scan1,
#                 perms=perm1, alpha=0.10,
#                 lty="dashed",lwd=2,col="purple")
# }


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

scan_cond <- scanone(f2g,  pheno.col=c("INS.10wk", "IL.6", "Il6_adipose"), addcovar=f2g$pheno$Mfsd7a, method="hk")

# note how we used "cbind" to concatenate the results. This is where we add in all the genes in our confidence interval as covariates.
# "cbind" calls the specialized function "cbind.scanone"
# Do rest of conditional scans. The scan above was just to get it started and to show how to do the first one manually.
# This part takes awhile.


for(i in 14:ncol(f2g$pheno)){
  scan_cond <- cbind(scan_cond, scanone(f2g, pheno.col=c("INS.10wk", "Il6_adipose", "IL.6"), addcovar=f2g$pheno[,i], method="hk"))
}
summary(scan_cond)

# Re-name the conditional scan columns with the gene symbol. Every 3 columns
# will reference a gene scanned conditionally on insulin, IL.6, and Il6_adipose in
# that order.

head(names(scan_cond))
dim(scan_cond)
for (i in 13:(ncol(f2g$pheno))) {
  names(scan_cond)[(3*(i-11)):(3*(i-11)+2)] <- names(f2g$pheno)[i]
}
head(names(scan_cond))
names(scan_cond)
summary(names(scan_cond))
#Plot all the conditional scans
# graphics.off()
# par(mfrow=c(6,1), mar=c(3,4,1,4) + 0.1)
# for( i in 1:ncol(scan_cond)){
#   plot(scan_cond, lodcolumn=i, main=names(scan_cond)[i+2])
#   add.threshold(scan1, perms=perm1, alpha=0.05, lty="dashed", lwd=2, col="orange")
#   add.threshold(scan1, perms=perm1, alpha=0.10, lty="dashed", lwd=2, col="purple")
# }
# graphics.on()

# Several genes drop the chromosome 2 peaks for all 3 phenotypes
# (insulin, Atg5, Atg7) down below significance. Find these genes in 
# the scan object. Subtract 2 from each index number to account 
# for the first two columns in the scan object, which are 
# chromosome number and cM position.


par(mfrow=c(6,1), mar=c(3,4,1,4) + 0.1)
for( i in 1:(ncol(scan_cond))){
  plot(scan_cond, lodcolumn=i, main=names(scan_cond)[i+2], chr =5)
  add.threshold(scan1, perms=perm1, alpha=0.05, lty="dashed", lwd=2, col="orange")
  add.threshold(scan1, perms=perm1, alpha=0.10, lty="dashed", lwd=2, col="purple")
}


my_scan_func <- function(x){
  par(mfrow=c(2, 3))
  
  
  plot(scan1, chr =5, lodcolumn = 7, main = "QTL for Clinical Insulin Plasma Levels at 10 Weeks", ylab = "Insulin Levels", ylim = c(0, 7))
  add.threshold(scan1, perms=perm1, alpha=0.05, lty="dashed", lwd=2, col="orange")
  add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed", lwd=2, col="purple")
  
  plot(scan1, chr =5,  lodcolumn = 4, main = "Quantitative Trait Loci for Il6 Expression in Adipose Tissue", ylab = "Il6 Adipose Expression", ylim = c(0, 7))
  add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed",lwd=2,col="orange")
  add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed",lwd=2,col="purple")
  
  plot(scan1, chr =5, lodcolumn = 3, main = "Quantitative Trait Loci for Clinical IL-6 Plasma Levels", ylab = "IL-6 Levels", ylim= c(0, 7))
  add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed",lwd=2,col="orange")
  add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed",lwd=2,col="purple")
  
  #Find lod column number of Chst1 insulin, Chst1, il6, and Chst1 il6 adipose exp in scan_cond
  
  plot(scan_cond, chr = 5, lodcolumn = 3*(x-1)+1, main = paste("Insulin Scan with",  names(scan_cond[3*x]), " Adipose Expression as Covariate"), ylab= "Insulin Levels", ylim = c(0, 7)) #insulin
  add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed", lwd=2, col="orange")
  add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed", lwd=2, col="purple")
  
  plot(scan_cond, chr = 5, lodcolumn = 3*(x-1)+2, main = paste("Il6 Adipose Exp Scan with",  names(scan_cond[3*x+1]), "Adipose Expression as Covariate"), ylab = "Il6 Adipose Exp", ylim = c(0, 7)) #il6 adipose
  add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed", lwd=2, col="orange")
  add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed", lwd=2, col="purple")
  
  
  plot(scan_cond, chr = 5, lodcolumn = 3*(x-1)+3, main = paste("IL-6 Scan with",  names(scan_cond[3*x+2]), "Adipose Expression as Covariate"), ylim = c(0, 7), ylab = "IL-6 Levels") #il6
  add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed", lwd=2, col="orange")
  add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed", lwd=2, col="purple")
}




# Which genes drop the chromosome 5 peak below the
# significance threshold of 3.67 for insulin,, Il6 Adipose Expression, AND Clinical IL-6 levels?
summary(subset(scan_cond, chr = 5),  thresholds = c(9.1, 7.1, 6.3, 6.3,3.3))

x = 1
lod_threshold = 3.67
candidate_genes <- vector()
while (2+3*x < length(summary(subset(scan_cond, chr = 5))))
{
  scores = vector()
  for (lod_score in summary(subset(scan_cond, chr = 5, threshold = c(9.1, 7.1, 6.3, 6.3, 3.3)))[(3*x+1):(3*x+3)]) #We want to skip the first five So that's why x starts at 2.
  {
    scores <- c(scores, lod_score)
    print(scores)
  }
  if (length(scores[scores < lod_threshold]) == 3)
  {
    print(names(scan_cond[3*x+1]))
    cat("Lod column numbers in scan_cond object: ", 3*(x-1)+1, " through, ", 3*(x-1)+3)
    print("")
    print("-----------------------")
    candidate_genes <- c(candidate_genes, names(scan_cond[3*x]))
    my_scan_func(x)
  }
  x = x + 1
}
print(candidate_genes)

# "Scfd2"    "Pdgfra"   "Lnx1"     "Thap6"    "Hsd17b13" "Cxcl5"  Cmklr1

summary(subset(scan_cond, chr = 5))[,which(names(scan_cond) %in% c("Scfd2", "Pdgfra", "Lnx1", "Thap6", "Hsd17b13", "Cxcl5", "Cmklr1"))]



# Pkd2 (IL-6, il6)

# Lnx1 (Insulin, Il6 adipose)

# Maybe BC005561 on IL-6, insulin
# Grxcr1 (insulin, IL-6)
# Tmem175 (INsulin, IL-6)
# Tmprss11f (Insulin, IL-6)
# Uso1 (Insulin, IL-6)
# Grxcr1 (Insulin, IL-6)
# BC005561 (Insulin, IL-6)
# Afp (insulin, IL-6)

# Cxcl5
# Hsd17b13
# Pdgfra
# Thap6



#######################################################################################
####################### More Scatterplot Matrices #####################################
#######################################################################################

#Define some variables

grep("Ins", annot$gene1)

#Primary Genes
Pdgfra_adipose <- adipose.rz[, annot[grep("Pdgfra", annot$gene1), 1]] #39.55 
Cxcl5_adipose <- adipose.rz[, annot[grep("Cxcl5", annot$gene1), 1]] #44.78
Hsd17b13_adipose <- liver.rz[, annot[grep("Hsd17b13", annot$gene1), 1]] #50.46
Thap6_adipose <- liver.rz[, annot[grep("Thap6", annot$gene1), 1]] #45.62
Cmklr1_adipose <- adipose.rz[, annot[grep("Cmklr1", annot$gene1), 1]] #55.55


#Secondary Genes
Tmem175_adipose <- adipose.rz[, annot[grep("Tmem175", annot$gene1), 1]]
Tmprss11f_adipose <- adipose.rz[, annot[grep("Tmprss11f", annot$gene1), 1]]
Uso1_adipose <- adipose.rz[, annot[grep("Uso1", annot$gene1), 1]]
Grxcr1_adipose <- adipose.rz[, annot[grep("Grxcr1", annot$gene1), 1]]
BC005561_adipose <- adipose.rz[, annot[grep("BC005561", annot$gene1), 1]]
Afp_adipose <- adipose.rz[, annot[grep("Afp", annot$gene1), 1]]



Ins_islet <- islet.rz[, annot[grep("Ins$", annot$gene1), 1]]



f2g$pheno <- cbind(f2g$pheno[,c("MouseNum","Sex","pgm")], phenotypes.rz[,c("INS.10wk", "IL.6")], Il6_adipose, Cxcl5_adipose, Hsd17b13_adipose, Pdgfra_adipose, Thap6_adipose, Cmklr1_adipose)
# look at all pairwise scatterplots of clinical and expression traits
names(f2g$pheno)
par(mfrow=c(1,1))
pairs(f2g$pheno[,4:length(f2g$pheno)], upper.panel=panel.cor,diag.panel=panel.hist)

#Pdgfra IL.6
#Pdgfra INS
#Pdgfra Il6_adipose

#Scatterplots show that Pdgfra seems to be only revelant one.

#Add primary AND secondary genes now.
f2g$pheno <- cbind(f2g$pheno[,c("MouseNum","Sex","pgm")], phenotypes.rz[,c("INS.10wk", "IL.6")], Il6_adipose, Cxcl5_adipose, Hsd17b13_adipose, Pdgfra_adipose, Thap6_adipose, Cmklr1_adipose)
                                      
#Tmem175_adipose, Tmprss11f_adipose, Uso1_adipose, Grxcr1_adipose, BC005561_adipose, Afp_adipose)

# The Q should be set to the closest marker to whatever gene we are looking at.

#39.55  44.78  50.46  45.62 55.55
f2g$pheno <- transform(f2g$pheno, Q5 = as.factor(f2g$geno[[5]]$data[,find.marker(f2g, 5, 44.78)]))
levels(f2g$pheno$Q5) <- c("B", "H", "R")
names(f2g$pheno)

# B -> 
# H -> Heterozygous
# R -> BTBR



####################################################################
#############           Effect Plots and BIC analysis ##############
######################################################################


# Analysis with Primary Candidate Genes (Genes that blocked all three peaks)

#Plot Pdgfra_adipose expression against Insulin


f2g$pheno <- transform(f2g$pheno, Hsd17b13_adipose)
par(mfrow=c(1,1))
qplot(Hsd17b13_adipose, INS.10wk, color=Q5, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q5), method="lm", se=FALSE)

f2g$pheno <- transform(f2g$pheno, Thap6_adipose)
par(mfrow=c(1,1))
qplot(Thap6_adipose, INS.10wk, color=Q5, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q5), method="lm", se=FALSE)


#Pdgfra. Nothing interesting
# Cxcl5_adipose. B genotype
# Hsd17b13. B genotype
# Thap6. All four diverge.

#Insulin
print("BIC Analysis with Insulin and IL-6. Primary Candidate Genes")


##################            Pdgfra
print("##################            Pdgfra     ##########################")
f2g$pheno <- transform(f2g$pheno, Q5 = as.factor(f2g$geno[[5]]$data[,find.marker(f2g, 5, 39.55)]))
levels(f2g$pheno$Q5) <- c("B", "H", "R")
print("##################            Pdgfra - Insulin     ##########################")
triple.fit(X = f2g$pheno$Pdgfra_adipose, Y = f2g$pheno$INS.10wk, Q = f2g$pheno$Q5)
print("##################            Pdgfra - IL6     ##########################")
triple.fit(X = f2g$pheno$Pdgfra_adipose, Y = f2g$pheno$IL.6, Q = f2g$pheno$Q5)

f2g$pheno <- transform(f2g$pheno, Pdgfra_adipose)
par(mfrow=c(2,1))
qplot(Pdgfra_adipose, INS.10wk, color=Q5, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q5), method="lm", se=FALSE)
qplot(Pdgfra_adipose, f2g$pheno$IL.6, color=Q5, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q5), method="lm", se=FALSE)


##################            Cxcl5
print("##################            Cxcl5     ##########################")
f2g$pheno <- transform(f2g$pheno, Q5 = as.factor(f2g$geno[[5]]$data[,find.marker(f2g, 5, 44.78)]))
levels(f2g$pheno$Q5) <- c("B", "H", "R")
print("##################            Cxcl5 - Insulin     ##########################")
triple.fit(X = f2g$pheno$Cxcl5_adipose, Y = f2g$pheno$INS.10wk, Q = f2g$pheno$Q5) #Reactive 
print("##################            Cxcl5 - IL6     ##########################")
triple.fit(X = f2g$pheno$Cxcl5_adipose, Y = f2g$pheno$IL.6, Q = f2g$pheno$Q5) #Reactive 


f2g$pheno <- transform(f2g$pheno, Cxcl5_adipose)
par(mfrow=c(2,1))
qplot(Cxcl5_adipose, INS.10wk, color=Q5, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q5), method="lm", se=FALSE)
qplot(Cxcl5_adipose, f2g$pheno$IL.6, color=Q5, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q5), method="lm", se=FALSE)


##################            Hsd16b13     ##########################
print("##################            Hsd16b13     ##########################")
f2g$pheno <- transform(f2g$pheno, Q5 = as.factor(f2g$geno[[5]]$data[,find.marker(f2g, 5, 50.46)]))
levels(f2g$pheno$Q5) <- c("B", "H", "R")
print("##################            Hsd16b13 - Insulin     ##########################")
triple.fit(X = f2g$pheno$Hsd17b13_adipose, Y = f2g$pheno$INS.10wk, Q = f2g$pheno$Q5) 
print("##################            Hsd16b13 - IL6     ##########################")
triple.fit(X = f2g$pheno$Hsd17b13_adipose, Y = f2g$pheno$IL.6, Q = f2g$pheno$Q5) 


f2g$pheno <- transform(f2g$pheno, Hsd17b13_adipose)
par(mfrow=c(2,1))
qplot(Hsd17b13_adipose, INS.10wk, color=Q5, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q5), method="lm", se=FALSE)
qplot(Hsd17b13_adipose, f2g$pheno$IL.6, color=Q5, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q5), method="lm", se=FALSE)

##################            Thap6
print("##################            Thap6     ##########################")
f2g$pheno <- transform(f2g$pheno, Q5 = as.factor(f2g$geno[[5]]$data[,find.marker(f2g, 5, 45.62)]))
levels(f2g$pheno$Q5) <- c("B", "H", "R")
print("##################            Thap6 - Insulin     ##########################")
triple.fit(X = f2g$pheno$Thap6_adipose, Y = f2g$pheno$INS.10wk, Q = f2g$pheno$Q5)
print("##################            Thap6 - Insulin     ##########################")
triple.fit(X = f2g$pheno$Thap6_adipose, Y = f2g$pheno$IL.6, Q = f2g$pheno$Q5)


f2g$pheno <- transform(f2g$pheno, Thap6_adipose)
par(mfrow=c(2,1))
qplot(Thap6_adipose, INS.10wk, color=Q5, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q5), method="lm", se=FALSE)
qplot(Thap6_adipose, f2g$pheno$IL.6, color=Q5, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q5), method="lm", se=FALSE)


##################            Cmklr1          
print("##################            Cmklr1     ##########################")
f2g$pheno <- transform(f2g$pheno, Q5 = as.factor(f2g$geno[[5]]$data[,find.marker(f2g, 5, 55.55)]))
levels(f2g$pheno$Q5) <- c("B", "H", "R")

print("##################            Cmklr1 - Insulin     ##########################")
triple.fit(X = f2g$pheno$Cmklr1_adipose, Y = f2g$pheno$INS.10wk, Q = f2g$pheno$Q5)
print("##################            Cmklr1 - IL6     ##########################")
triple.fit(X = f2g$pheno$Cmklr1_adipose, Y = f2g$pheno$IL.6, Q = f2g$pheno$Q5)


f2g$pheno <- transform(f2g$pheno, Cmklr1_adipose)
par(mfrow=c(2,1))
qplot(Cmklr1_adipose, INS.10wk, color=Q5, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q5), method="lm", se=FALSE)
qplot(Cmklr1_adipose, f2g$pheno$IL.6, color=Q5, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q5), method="lm", se=FALSE)

#####################################################################################













##########################################################################################################################

# Analysis with Secondary Candidate Genes (Genes that blocked Insulin and IL-6 peaks, but not Il6 adipose gene expression)

# Tmem175 (INsulin, IL-6)
# Tmprss11f (Insulin, IL-6)
# Uso1 (Insulin, IL-6)
# Grxcr1 (Insulin, IL-6)
# BC005561 (Insulin, IL-6)
# Afp (insulin, IL-6)


f2g$pheno <- transform(f2g$pheno, Tmprss11f_adipose)
par(mfrow=c(1,1))
qplot(Tmprss11f_adipose, INS.10wk, color=Q5, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q5), method="lm", se=FALSE)

f2g$pheno <- transform(f2g$pheno, Uso1_adipose)
par(mfrow=c(1,1))
qplot(Uso1_adipose, INS.10wk, color=Q5, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q5), method="lm", se=FALSE)

f2g$pheno <- transform(f2g$pheno, Grxcr1_adipose)
par(mfrow=c(1,1))
qplot(Grxcr1_adipose, INS.10wk, color=Q5, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q5), method="lm", se=FALSE)

f2g$pheno <- transform(f2g$pheno, BC005561_adipose)
par(mfrow=c(1,1))
qplot(BC005561_adipose, INS.10wk, color=Q5, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q5), method="lm", se=FALSE)

f2g$pheno <- transform(f2g$pheno, Afp_adipose)
par(mfrow=c(1,1))
qplot(Afp_adipose, INS.10wk, color=Q5, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q5), method="lm", se=FALSE)


print("BIC analyses with secondary candidate genes")


##################            Tmem175     ##########################
triple.fit(X = f2g$pheno$Tmem175_adipose, Y = f2g$pheno$INS.10wk, Q = f2g$pheno$Q5)
triple.fit(X = f2g$pheno$Tmem175_adipose, Y = f2g$pheno$IL.6, Q = f2g$pheno$Q5)

f2g$pheno <- transform(f2g$pheno, Tmem175_adipose)
par(mfrow=c(2,1))
qplot(Tmem175_adipose, INS.10wk, color=Q5, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q5), method="lm", se=FALSE)
qplot(Tmem175_adipose, f2g$pheno$IL.6, color=Q5, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q5), method="lm", se=FALSE)

##################            Tmprss11f     ##########################
triple.fit(X = f2g$pheno$Tmprss11f_adipose, Y = f2g$pheno$INS.10wk, Q = f2g$pheno$Q5) #Reactive 6.595
triple.fit(X = f2g$pheno$Tmprss11f_adipose, Y = f2g$pheno$IL.6, Q = f2g$pheno$Q5)

f2g$pheno <- transform(f2g$pheno, Tmprss11f_adipose)
par(mfrow=c(2,1))
qplot(Tmprss11f_adipose, INS.10wk, color=Q5, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q5), method="lm", se=FALSE)
qplot(Tmprss11f_adipose, f2g$pheno$IL.6, color=Q5, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q5), method="lm", se=FALSE)

##################            Uso1     ##########################
triple.fit(X = f2g$pheno$Uso1_adipose, Y = f2g$pheno$INS.10wk, Q = f2g$pheno$Q5)
triple.fit(X = f2g$pheno$Uso1_adipose, Y = f2g$pheno$IL.6, Q = f2g$pheno$Q5)

f2g$pheno <- transform(f2g$pheno, Uso1_adipose)
par(mfrow=c(2,1))
qplot(Uso1_adipose, INS.10wk, color=Q5, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q5), method="lm", se=FALSE)
qplot(Uso1_adipose, f2g$pheno$IL.6, color=Q5, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q5), method="lm", se=FALSE)

##################            Cmklr1     ##########################
triple.fit(X = f2g$pheno$Grxcr1_adipose, Y = f2g$pheno$INS.10wk, Q = f2g$pheno$Q5) #Reactive 5.555
triple.fit(X = f2g$pheno$Grxcr1_adipose, Y = f2g$pheno$IL.6, Q = f2g$pheno$Q5)

f2g$pheno <- transform(f2g$pheno, Grxcr1_adipose)
par(mfrow=c(2,1))
qplot(Grxcr1_adipose, INS.10wk, color=Q5, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q5), method="lm", se=FALSE)
qplot(Grxcr1_adipose, f2g$pheno$IL.6, color=Q5, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q5), method="lm", se=FALSE)

##################            BC005561     ##########################
triple.fit(X = f2g$pheno$BC005561_adipose, Y = f2g$pheno$INS.10wk, Q = f2g$pheno$Q5)
triple.fit(X = f2g$pheno$BC005561_adipose, Y = f2g$pheno$IL.6, Q = f2g$pheno$Q5)

f2g$pheno <- transform(f2g$pheno, BC005561_adipose)
par(mfrow=c(2,1))
qplot(BC005561_adipose, INS.10wk, color=Q5, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q5), method="lm", se=FALSE)
qplot(BC005561_adipose, f2g$pheno$IL.6, color=Q5, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q5), method="lm", se=FALSE)

##################            Afp_adipose     ##########################
triple.fit(X = f2g$pheno$Afp_adipose, Y = f2g$pheno$INS.10wk, Q = f2g$pheno$Q5)
triple.fit(X = f2g$pheno$Afp_adipose, Y = f2g$pheno$IL.6, Q = f2g$pheno$Q5)

f2g$pheno <- transform(f2g$pheno, Afp_adipose)
par(mfrow=c(2,1))
qplot(Afp_adipose, INS.10wk, color=Q5, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q5), method="lm", se=FALSE)
qplot(Afp_adipose, f2g$pheno$IL.6, color=Q5, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q5), method="lm", se=FALSE)






##########################################################
triple.fit(X = f2g$pheno$Cxcl5_adipose, Y = f2g$pheno$INS.10wk, Q = f2g$pheno$Q5) #Reactive. 7.656
triple.fit(X = f2g$pheno$Cxcl5_adipose, Y = f2g$pheno$IL.6, Q = f2g$pheno$Q5) #Reactive. 7.656
##########################################################
# check 4 conditions for Cxcl5 gene expression 
# as a mediator of Q5 effect on insulin
# Normally this is done before doing BIC modeling, but we are going to do it in reverse.
#If the conditions are met, we have a case for a pathway.

#davidakenny.net/cm/mediate.htm
####  i) Insulin is linked to Q5                    # X -> Y
anova(lm(INS.10wk ~ Sex + Q5, data = f2g$pheno))
# significant

####  ii) Clxcl5 gene expression is linked to Q5     # X -> M
anova(lm(Cxcl5_adipose ~ Sex + Q5, data = f2g$pheno))
# significant

####  iii) Insulin not linked after accounting for Q5   
anova(lm(INS.10wk ~ Sex + Cxcl5_adipose + Q5, data = f2g$pheno))
# not significant * .

####  iv) Cxc gene expression is still linked after 
# accounting for insulin
anova(lm(Cxcl5_adipose ~ Sex + INS.10wk + Q5, data = f2g$pheno))
# significant ***

# all 4 conditions for a mediator are satisfied
#######################################################################3


##########################################################
####  II. Mediation analysis
#  use linear models to compute mediation on the peak chr 8 marker
#  this is just another way to get at the same question

###
# get the genotypes at peak marker for insulin
#
# fill in the missing genotypes
f2g <- fill.geno(f2g, method="argmax")
#
#create a three level factor Q2 and add to pheno data
f2g$pheno <- transform(f2g$pheno,Q2 = as.factor(f2g$geno[[2]]$data[,find.marker(f2g, 2, 75.2)]))
levels(f2g$pheno$Q2) <- c("B","H","R")
names(f2g$pheno)
# plot insulin against Pdrg1 expression
pdrg1_islet <- islet.rz[, annot[grep("Pdrg1", annot$gene1), 1]]
f2g$pheno <- transform(f2g$pheno, pdrg1_islet)
par(mfrow=c(1,1))
qplot(pdrg1_islet, INS.10wk, color=Q2, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q2), method="lm", se=FALSE)
# There's reasonably good correlation between genotype classes.



#####
# check 4 conditions for Pdrg1 gene expression 
# as a mediator of Q2 effect on insulin

#davidakenny.net/cm/mediate.htm

####  i) Insulin is linked to Q2                    # X -> Y
anova(lm(INS.10wk ~ Sex + Q2, data = f2g$pheno))
# significant

####  ii) Pdrg1 gene expression is linked to Q2     # X -> M
anova(lm(pdrg1_islet ~ Sex + Q2, data = f2g$pheno))
# significant

####  iii) Insulin not linked after accounting for Q2   
anova(lm(INS.10wk ~ Sex + pdrg1_islet + Q2, data = f2g$pheno))
# not significant * .

####  iv) Pdrg1 gene expression is still linked after 
# accounting for insulin
anova(lm(pdrg1_islet ~ Sex + INS.10wk + Q2, data = f2g$pheno))
# significant ***

# all 4 conditions for a mediator are satisfied

###############################################
#  III. establish mediator using model selection with BIC scoring

# note that missing values will mess up BIC analysis
apply(is.na(f2g$pheno), 2, sum)
#
# easiest to use the triple.fit function that removes missing data
# and fits the three models

# compute BIC scores for the causal triplet models
#   using pdrg1_islet
with(f2g$pheno,triple.fit(pdrg1_islet, INS.10wk, Q2))
#  causal model has lowest score
#  suggests that pdrg1_islet is a strong mediator of insulin

# compute BIC scores for the causal triplet models
#   using atg5_islet and atg7_islet
with(f2g$pheno,triple.fit(pdrg1_islet, atg5_islet, Q2))
#  causal model has lowest score
#  suggests that pdrg1_islet is a strong mediator of Atg5 expression in islet

#  causal model has lowest score but differs by only 5 from complex model
#  therefore inconclusive



###############################################





























########## random stuff. code that was taken out.



#########################################################################################


#par(mfrow=c(3, 3))


# plot(scan1, chr =5, lodcolumn = 7, main = "QTL for Clinical Insulin Plasma Levels at 10 Weeks", ylab = "Insulin Levels", ylim = c(0, 7))
# add.threshold(scan1, perms=perm1, alpha=0.05, lty="dashed", lwd=2, col="orange")
# add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed", lwd=2, col="purple")
# 
# plot(scan1, chr =5,  lodcolumn = 4, main = "Quantitative Trait Loci for Il6 Expression in Adipose Tissue", ylab = "Il6 Adipose Expression", ylim = c(0, 7))
# add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed",lwd=2,col="orange")
# add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed",lwd=2,col="purple")
# 
# plot(scan1, chr =5, lodcolumn = 3, main = "Quantitative Trait Loci for Clinical IL-6 Plasma Levels", ylab = "IL-6 Levels", ylim= c(0, 7))
# add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed",lwd=2,col="orange")
# add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed",lwd=2,col="purple")
# 
# #Pdgfra
# 
# plot(scan_cond, chr = 5, lodcolumn = 106, main = paste("Insulin Scan with Pdgfra Adipose Expression as Covariate"), ylab= "Insulin Levels", ylim = c(0, 7)) #insulin
# add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed", lwd=2, col="orange")
# add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed", lwd=2, col="purple")
# 
# plot(scan_cond, chr = 5, lodcolumn = 107, main = paste("Il6 Adipose Exp Scan with Pdgfra Adipose Expression as Covariate"), ylab = "Il6 Adipose Exp", ylim = c(0, 7)) #il6 adipose
# add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed", lwd=2, col="orange")
# add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed", lwd=2, col="purple")
# 
# plot(scan_cond, chr = 5, lodcolumn = 108, main = paste("IL-6 Scan with Pdgfra Adipose Expression as Covariate"), ylim = c(0, 7), ylab = "IL-6 Levels") #il6
# add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed", lwd=2, col="orange")
# add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed", lwd=2, col="purple")
# 
# #Thap6
# plot(scan_cond, chr = 5, lodcolumn = 307, main = paste("Insulin Scan with Thap6 Adipose Expression as Covariate"), ylab= "Insulin Levels", ylim = c(0, 7)) #insulin
# add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed", lwd=2, col="orange")
# add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed", lwd=2, col="purple")
# 
# plot(scan_cond, chr = 5, lodcolumn = 308, main = paste("Il6 Adipose Exp Scan with Thap6 Adipose Expression as Covariate"), ylab = "Il6 Adipose Exp", ylim = c(0, 7)) #il6 adipose
# add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed", lwd=2, col="orange")
# add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed", lwd=2, col="purple")
# 
# plot(scan_cond, chr = 5, lodcolumn = 309, main = paste("IL-6 Scan with Thap6 Adipose Expression as Covariate"), ylim = c(0, 7), ylab = "IL-6 Levels") #il6
# add.threshold(scan1, perms=perm1, alpha=0.05,lty="dashed", lwd=2, col="orange")
# add.threshold(scan1, perms=perm1, alpha=0.10,lty="dashed", lwd=2, col="purple")




# #Some practice with 
# lm(formula = f2g$pheno$INS.10wk ~ sex)
# summary(lm(formula = f2g$pheno$INS.10wk ~ sex))
# BIC(lm(formula = f2g$pheno$INS.10wk ~ sex))
# 
# 
# lm(formula = f2g$pheno$INS.10wk ~ sex + f2g$pheno$Chst1)
# summary(lm(formula = f2g$pheno$INS.10wk ~ sex + f2g$pheno$Chst1))
# BIC(lm(formula = f2g$pheno$INS.10wk ~ sex + f2g$pheno$Chst1))
# #This model is better. (It has a lower BIC score). A difference of at least 5 is needed to differentiate between models.
# 
# 
# #We can automate this process

#If the complex model is smallest, we are at the dead end cuz we have to assume that is the best model.

# f2g$pheno <- cbind(f2g$pheno[,c("MouseNum","Sex","pgm")], phenotypes.rz[,c("INS.10wk", "GLU.10wk", "IL.6")], Chst1_adipose, Chst1_gastroc,  Il6_adipose, Il6ra_adipose)
# par(mfrow=c(1,1))
# pairs(f2g$pheno[,4:length(f2g$pheno)], upper.panel=panel.cor,diag.panel=panel.hist)


