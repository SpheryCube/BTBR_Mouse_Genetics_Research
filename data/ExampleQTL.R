## Code illustrating the QDG and QTLnet routines of the R/qtlnet package 

############################################################################
## Load package as R library and load environmental variables
############################################################################

rm(list=ls())
directory <- "/home/daniel14/CompBioProjects/BTBRxB6/data"
load("BTBR.clean.data.Rdata")
#Fix Il.6 data. Replaces every data point with less than 0 as an NA
phenotypes.rz$IL.6[phenotypes.rz$IL.6 < 0 & is.numeric(phenotypes.rz$IL.6)] <- NA

library(qtlnet)


############################################################################
## Edit f2g object
############################################################################

set.seed(12345)

# Map <- sim.map(len = rep(100, 5), n.mar = 11, eq.spacing = TRUE, 
#                include.x = FALSE)
# f2g <- sim.f2g(map = Map, n.ind = 500, type = "f2")
# f2ges <- vector(mode = "list", length = 5)
# add.effects <- c(1, 1, 1, 1, 1)
# for (i in 1:5) {
#   map <- sim.map(len = rep(100, i), n.mar = 11, eq.spacing = TRUE, 
#                  include.x = FALSE)
#   f2ges[[i]] <- sim.f2g(map = map, n.ind = 500, type = "f2", 
#                             model = c(i, 50, add.effects[i], 0))
#   f2g$geno[[i]] <- f2ges[[i]]$geno[[i]]
# }

beta <- 1   
# f2g$pheno[, 1] <- f2ges[[1]]$pheno
# f2g$pheno[, 2] <- f2ges[[2]]$pheno + beta * f2g$pheno[, 1]
# f2g$pheno[, 3] <- f2ges[[3]]$pheno + beta * f2g$pheno[, 2]
# f2g$pheno[, 4] <- f2ges[[4]]$pheno + beta * f2g$pheno[, 2] 
# f2g$pheno[, 5] <- f2ges[[5]]$pheno + beta * f2g$pheno[, 3] + 
#   beta * f2g$pheno[,4]

f2g$pheno <- cbind(f2g$pheno[,c("MouseNum","Sex","pgm")], ppil.expr,
                   phenotypes.rz[,c("CHOL","vWF","LDL","liver.TG")])

names(f2g$pheno) <- paste("y", 1:5, sep = "")


############################################################################
## Determine LOD score permutation threshold
############################################################################

f2g <- calc.genoprob(f2g, step = 1)
set.seed(12345)
perm.test <- scanone(f2g, n.perm = 1000, method = "hk")
summary(perm.test)


############################################################################
## QDG routines ############################################################ 
############################################################################

############################################################################
## QDG routines - perform QTL mapping and create QTL objects
############################################################################

Scan <- scanone(f2g, pheno.col = 1:5, method = "hk")

f2g <- sim.geno(f2g, n.draws = 1)
marker.nms <- allqtls <- vector(mode = "list", length = 5)
names(marker.nms) <- names(allqtls) <- paste("y", 1:5, sep = "")
for (i in 1:5) {
  aux <- summary(Scan[, c(1, 2, i + 2)], thr = 3.04)
  marker.nms[[i]] <- find.marker(f2g, chr = aux[, 1], pos = aux[, 2])
  allqtls[[i]] <- makeqtl(f2g, chr = aux[, 1], pos = aux[, 2])
}


############################################################################
## QDG routines - fit qdg algorithm
############################################################################

out1 <- qdg(f2g = f2g,
            phenotype.names = paste("y", 1:5, sep = ""),
            marker.names = marker.nms,
            QTL = allqtls,
            alpha = 0.005,
            n.qdg.random.starts = 10,
            addcov = NULL, 
            intcov = NULL,
            skel.method = "pcskel")


############################################################################
## QDG routines - output and plot
############################################################################

out1$UDG
out1$DG
out1$Solutions

gr1 <- graph.qdg(out1, include.qtl = FALSE)
plot(gr1)
gr2 <- graph.qdg(out1, include.qtl = TRUE)
plot(gr2)


############################################################################
## Checking conditional and unconditional QTL mapping results ############## 
############################################################################

#postscript("uncond.cond.QTL.eps", width = 14, height = 8)
par(mfrow = c(2, 5), cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
uncond.nms <- paste("Y", 1:5, sep = "")
for (i in 1:5) {
  plot(Scan, lodcolumn = i, main = uncond.nms[i], ylab = "lod")
}
plot(Scan, lodcolumn = 1, main = uncond.nms[1], ylab = "lod")
cond.nms <- c("Y1", "Y2 | Y1", "Y3 | Y2", "Y4 | Y2", "Y5 | Y3, Y4")
pheno.parents <- list(NULL, 1, 2, 2, c(3, 4))
for (i in 2:5) {
  CondScan <- scanone(f2g, pheno.col = i, method = "hk",
                      addcov = f2g$pheno[, pheno.parents[[i]]])
  plot(CondScan, main = cond.nms[i], ylab = "lod")
}
#dev.off()

############################################################################
## QTLnet routines ######################################################### 
############################################################################

############################################################################
## QTLnet routines - basic functionality - fit QTLNet algorithm
############################################################################

out2 <- mcmc.qtlnet(f2g = f2g,
                    pheno.col = 1:5, 
                    threshold = 3.04, 
                    addcov = NULL, 
                    intcov = NULL,
                    nSamples = 1000, 
                    thinning = 3, 
                    max.parents = 4, 
                    M0 = NULL,
                    burnin = 0.2, 
                    method = "hk", 
                    random.seed = 987654321, 
                    init.edges = 0,
                    saved.scores = NULL, 
                    rev.method = "nbhd",
                    verbose = TRUE)


############################################################################
## QTLnet routines - basic functionality - inspect results
############################################################################

summary(out2)
print(out2)
loci.qtlnet(out2)
plot(out2)
par(mfrow = c(1, 1))
plotbic.qtlnet(out2, smooth = FALSE)


############################################################################
## QTLnet routines - "parallel functionality"
############################################################################

######################################################
## STEP 1: Preparation. Fast. Needed in steps 2 and 3.

pheno.col <- 1:5
max.parents <- 4
size.qtlnet(pheno.col, max.parents)
parents <- parents.qtlnet(pheno.col, max.parents)
groups <- group.qtlnet(parents = parents, group.size = 10)

save(f2g, pheno.col, max.parents, parents, groups,
     file = "Step1.RData", compress = TRUE)

parents
groups

pa <- summary(parents)
g <- nrow(groups)
N <- rep(NA, g)
for (i in 1:g) {
  N[i] <- sum(pa[seq(groups[i, 1], groups[i, 2]), 2])
}
N

######################################################
## STEP 2: Pre-compute BIC scores for selected parents.

load("Step1.RData")
for (i in seq(nrow(groups))) {
  bic <- bic.qtlnet(f2g, 
                    pheno.col, 
                    threshold = 3.04,
                    max.parents = max.parents,
                    parents = parents[seq(groups[i,1], groups[i,2])])
  
  save(bic, file = paste("bic", i, ".RData", sep = ""), compress = TRUE)
}

## Read in saved BIC scores and combine into one object.

load("Step1.RData")
bic.group <- list()
for (i in seq(nrow(groups))) {
  load(paste("bic", i, ".RData", sep = ""))
  bic.group[[i]] <- bic
  cat("group =", i, "\n")
}
saved.scores <- bic.join(f2g, pheno.col, bic.group, max.parents = 4)
saved.scores


######################################
## STEP 3: Sample Markov chain (MCMC).

set.seed(54321)
n.runs <- 3
for (i in seq(n.runs)) {
  cat("run =", i, "\n")
  ## Run MCMC with randomized initial network.
  mcmc <- mcmc.qtlnet(f2g, 
                      pheno.col, 
                      threshold = 3.04, 
                      thinning = 1,
                      max.parents = max.parents, 
                      saved.scores = saved.scores, 
                      init.edges = NULL)
  
  save(mcmc, file = paste("mcmc", i, ".RData", sep = ""), 
       compress = TRUE)
}


###############################################
## STEP 4: Combine results for post-processing.

n.runs <- 3
outs.qtlnet <- list()
for (i in seq(n.runs)) {
  load(paste("mcmc", i, ".RData", sep = ""))
  outs.qtlnet[[i]] <- mcmc
}
out3 <- c.qtlnet(outs.qtlnet)


############################
## Inspect combined results. 

summary(out3)
plot(out3)
plotbic.qtlnet(out3, smooth = FALSE)