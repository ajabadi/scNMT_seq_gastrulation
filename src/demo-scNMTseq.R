library(MultiAssayExperiment)

gastrulation_mae <- readRDS(here('output/scnmtseq_gastrulation_mae_826-cells_orderedFeatures.rds'))  # MAE experiment
keep_cells <- (gastrulation_mae$lineage10x_2 %in% c('Mesoderm', 'Endoderm', 'Ectoderm')) & (gastrulation_mae$stage == 'E7.5')
gastrulation_E7.5 <- gastrulation_mae[, keep_cells,1:13]

gastrulation_E7.5 <- make_unique_feature_names(gastrulation_E7.5)

gastrulation_E7.5_subset <- MultiAssayExperiment(experiments = lapply(experiments(gastrulation_E7.5), function(x) subset_pn(x, n = 300)), colData = colData(gastrulation_E7.5))

## matching cells only
 gastrulation_E7.5_subset <- MatchedAssayExperiment( gastrulation_E7.5_subset)
## Unique cell/sample names required

## look at available assays amd sample data
names(colData( gastrulation_E7.5_subset))
names(experiments( gastrulation_E7.5_subset))


## ----------------------------- Unsupervised ----------------------------- ##
ncomp <- 2
keep_features <- list(rna = c(20, 10),
                      acc_promoter = c(5, 15),
                      met_promoter = c(12, 20))

mmsplsda <- MultiModalSparsePLS(data =  gastrulation_E7.5_subset, formula = rna ~ acc_promoter + met_promoter, 
                                ncomp = ncomp, design = 'null', keep_features = keep_features)

plotIndiv(mmsplsda, pch = 16, group =  gastrulation_E7.5_subset$stage_lineage, legend = TRUE, legend.title = 'Protein Cluster')


## ------------------------------ Supervised ------------------------------ ##
mmsplsda <- MultiModalSparsePLSDA(data =  gastrulation_E7.5_subset, formula = protein.cluster ~ RNASeq2GeneNorm + RPPAArray + miRNASeqGene, 
                                  ncomp = ncomp, design = 'full', keep_features = keep_features)

plotIndiv(mmsplsda, pch = 16, group =  gastrulation_E7.5_subset$protein.cluster, legend = TRUE, legend.title = 'Protein Cluster')
circosPlot(mmsplsda, cutoff = 0.6)
