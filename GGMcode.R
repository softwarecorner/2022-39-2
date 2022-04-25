# needed package
library("rags2ridges")

# load, extract and scale data for AD Class 2
data("ADdata")
ADclass2 <- scale(t(ADmetabolites[,sampleInfo$ApoEClass=="Class 2"]))

# precision matrix estimation with given penalty value**
P <- ridgeP(covML(ADclass2), lambda = .15)

# extract Network**
P0 <- sparsify(P, threshold="localFDR", FDRcut=.999)

# visualize Network with node-coloring**
PcorP  <- pruneMatrix(P0$sparseParCor)
Colors <- rownames(PcorP)
Colors[grep("Amine",     rownames(PcorP))] <- "lightblue"
Colors[grep("Org.Acid",  rownames(PcorP))] <- "orange"
Colors[grep("Lip",       rownames(PcorP))] <- "yellow"
Colors[grep("Ox.Stress", rownames(PcorP))] <- "purple"

# plot network
Ugraph(PcorP, type="fancy", lay="layout_with_fr", 
       Vcolor=Colors, Vsize=7, Vcex=.3)

# find and visualize the communities for the extracted network
Commy <- Communities(PcorP, Vcolor=Colors, Vsize=7,  Vcex=.3)

