//Remove data from a script

> hist(phenotypes$IL.6, breaks=100)
> which(phenotypes$IL.6 > 400)
[1] 126
> hist(phenotypes$IL.6[c(1:125, 127:518)], breaks=100)
> which(phenotypes$IL.6 = 0)
Error: unexpected '=' in "which(phenotypes$IL.6 ="
> which(phenotypes$IL.6 < 1)
> which(phenotypes$IL.6 < 1)
> s = which(phenotypes$IL.6 < 1)
> newphenotypesIL <- which(phenotypes$IL.6 < 0.1)
> hist(newphenotypesIL)