# BiocManager::install("mixOmicsTeam/mixOmics@MultiAssayExperiment")
library(MultiAssayExperiment)
data("miniACC")

## matching cells only
miniACC <- MatchedAssayExperiment(miniACC)
## Unique cell/sample names required
miniACC <- make_unique_feature_names(miniACC)

## look at available assays amd sample data
names(colData(miniACC))
names(experiments(miniACC))


## ----------------------------- Unsupervised ----------------------------- ##
ncomp <- 2
keep_features <- list(RNASeq2GeneNorm = c(20, 10),
                      RPPAArray = c(5, 15),
                      miRNASeqGene = c(12, 20))

mmsplsda <- MultiModalSparsePLS(data = miniACC, formula = RNASeq2GeneNorm ~ RPPAArray + miRNASeqGene, 
                                ncomp = ncomp, design = 'full', keep_features = keep_features)

plotIndiv(mmspls, pch = 16, group = miniACC$protein.cluster, legend = TRUE, legend.title = 'Protein Cluster')


## ------------------------------ Supervised ------------------------------ ##
mmsplsda <- MultiModalSparsePLSDA(data = miniACC, formula = protein.cluster ~ RNASeq2GeneNorm + RPPAArray + miRNASeqGene, 
                                ncomp = ncomp, design = 'full', keep_features = keep_features)

plotIndiv(mmsplsda, pch = 16, group = miniACC$protein.cluster, legend = TRUE, legend.title = 'Protein Cluster')
circosPlot(mmsplsda, cutoff = 0.6)
