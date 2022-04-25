# load package and data
library(ragt2ridges)
library(Biobase)
data(hpvP53)

# reformat and zero center data
Y <- centerVAR1data(longitudinal2array(t(exprs(hpvP53rna))))

# fit the model
VAR1hat <- ridgeVAR1(Y=Y, lambdaA=100, lambdaP=1)

# support determination 
zerosA <- sparsifyVAR1(A=VAR1hat$A, SigmaE=symm(solve(VAR1hat$P)), 
                       threshold="top", top=50)$zeros
VAR1hat$A[zerosA] <- 0
VAR1hat$P <- sparsify(VAR1hat$P, threshold="top", top=10)$sparseParCor

# plot time-series chain graph
graphVAR1(VAR1hat$A, VAR1hat$P, nNames=featureNames(hpvP53rna))

# motif detection
motifStatsVAR1(VAR1hat$A)

