###########################################################
#  This is a script to demo the use of Attie BTBR eQTL data
#  This script looks at lipid phenotypes
#  updated January 7, 2014 - GAC
###########################################################
# 
#set working directory
#setwd("/Users/garyc/GARY/Projects/180_BTBR_eQTL/BTBR_data")

#load qtl library
library(qtl)

#load graphics package  (http://had.co.nz/ggplot2/)
library(ggplot2)

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
####end useful plotting function

##########################################################
#read cross data into R environment
load("BTBR.clean.data.Rdata")
ls()

names(phenotypes.rz)
#we are going to look at insulin and interleukin-6 levels.

#set up the cross object with out phenotypes
summary(f2g)
names(f2g$pheno)
f2g$pheno <- cbind(f2g$pheno[,c("MouseNum","Sex","pgm")],
                   phenotypes.rz[,c("INS.10wk","IL.6","Urine.Volume")])
names(f2g$pheno)

#look at histrograms of data by sex
X11()  #quartz() for windows.
qplot(INS.10wk, facets=Sex~., data=f2g$pheno)
X11()
qplot(IL.6, facets=Sex~., data=f2g$pheno)
X11()
qplot(Urine.Volume, facets=Sex~., data=f2g$pheno)



# note higher CHOL and HDL in males, also LDL but less so

#look at raw data cholesterol traits
#log scaling is helpful
X11()
pairs(log(phenotypes[,c("INS.10wk","IL.6","Urine.Volume")]),
      upper.panel=panel.cor,diag.panel=panel.hist)

#normally LDL is computed trait
quartz()
qplot(CHOL-HDL, LDL, data=phenotypes)
quartz()
qplot(LDL+HDL, CHOL, data=phenotypes)		
#LDL appears to be directly measured here
#ask Mark about this

#look at transformed cholesterol traits
quartz()
pairs(f2g$pheno[,c("LDL","HDL","CHOL")],
      upper.panel=panel.cor,diag.panel=panel.hist)
#all three are positively correlated
#HDL and CHOL are tightly correlated 