###########################################################
#  This is a script to demo the use of Attie BTBR eQTL data
#  This script looks at lipid phenotypes
#  updated January 7, 2014 - GAC
###########################################################
# 
#set working directory
#setwd("/Users/garyc/GARY/Projects/180_BTBR_eQTL/BTBR_data")


# Clear environmental variables
rm(list=ls())
directory <- "/home/daniel14/CompBioProjects/BTBRxB6/data"
setwd(directory)
getwd()
# Load data.
load(file = "BTBR.clean.data.Rdata")
#Fix Il.6 data. Replaces every data point with less than 0 as an NA
phenotypes.rz$IL.6[phenotypes.rz$IL.6 < 0 & is.numeric(phenotypes.rz$IL.6)] <- NA

#load qtl library
library(qtl)

#load graphics package  (http://had.co.nz/ggplot2/)
library(ggplot2)

####
# some useful plotting functions fron Gary
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
####end useful plotting functions


#plot(phenotypes.rz$IL.6, phenotypes.rz$insulin.10wk)
#str(phenotypes.rz$Il.6)
phenotypes.rz$IL.6[phenotypes.rz$IL.6 < 0 & is.numeric(phenotypes.rz$IL.6)] <- NA #Fix Il.6 data. Replaces every data point with less than 0 as an NA
#plot(phenotypes.rz$IL.6, phenotypes.rz$insulin.10wk)


quartz()     #For windows. Use X11() for ubuntu.
pairs(f2g$pheno[,c("INS.10wk","IL.6","Urine.Volume")], upper.panel=panel.cor,diag.panel=panel.hist) #Using shortname


####################################################################################
##############  QTL analysis on IL-6 levels ########################################
####################################################################################

f2g$pheno <- f2g$pheno[,c("MouseNum", "Sex", "pgm")]
my.phenos <- phenotypes.rz[,c("IL.6", "INS.10wk", "GLU.10wk")]

my.phenos
head(my.phenos)
f2g$pheno <- f2g$pheno[]

f2g$pheno <- cbind(f2g$pheno, my.phenos) #Binds my.phenos with whats already f2g$pheno

names(f2g$pheno) #Let's see if we are right.

#We have markers for the genotypes at certain places across the genenome for each chromosome

library(qtl)  #Make sure the qtl package is already installed.

#Quantiative trait mapping procedure
f2g <- calc.genoprob(f2g, step = 1) #Every 1 centimorgan calculate the gene probability
#We aren't setting f2g equal to something else really, rather we are more or less modifying it.
#Now we can get on to actually doing QTL mapping
sex <- as.numeric(f2g$pheno$Sex)
my.first.scan <- scanone(f2g, pheno.col = 4:6, addcovar = sex, method = "hk") #We start at 4 since mousenum sex and pgm aren't phenotypes
#But we want to seperate by sex (covariate)
#Method = haley knot

#We want to plot it and use a for loop

for (i in 1:3)  #This gives us some nice graphs. (3 graphs, to be exact. One for each phenotype)
{
  plot(my.first.scan, lodcolumn = i)
}
#End of first set of QTL analyses
####################################################################################
####################################################################################

#QTL analysis on insr levels in gastronic tissue and IR-6 levels


#Run QTL analysis on the ratio between insulin levels and IL-6 levels

#Make scatterplot matrices for IL6, insulin levels, and IR expression in gastronomous and adipose tissue
X11()

