load("BTBR.clean.data.Rdata"))
names(f2g$pheno)
f2g$pheno <- f2g$pheno[,c("MouseNum", "Sex", "pgm")]
my.phenos <- phenotypes.rz[,c("TG.homogenate", "INS.4wk", "GLU.4wk")]

my.phenos
head(my.phenos)
f2g$pheno <- f2g$pheno[]

f2g$pheno <- cbind(f2g$pheno, my.phenos) #Binds my.phenos with whats already f2g$pheno

names(f2g$pheno) #Let's see if we are right.

  #We have markers for the genotypes at certain places across the genenome for each chromosome
  #However, we

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

#We need to permute the genotypes and the phenotypes (mix up the data)
#We need to do this to find a LOD threshold to see what's significant on our actualy data.

my.first.scan <- scanone(f2g, pheno.col = 4:6, addcovar = sex, method = "hk", n.perm = 100, perm.Xsp = TRUE) 
#Don't get scared when X chromosome permutations go over 100
for (i in 1:3)  #This gives us some nice graphs. (3 graphs, to be exact. One for each phenotype)
{
  plot(my.first.scan, lodcolumn = i)
}


