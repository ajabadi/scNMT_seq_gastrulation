library(MultiAssayExperiment)
library(here)
library(ggplot2)
library(knitr)

mmsplsda <- multimodal_analysis_wrapper(mae = gastru.mae, study_assays = NULL, ncomp = 3, scale = FALSE, design = 'null', lineages = NULL, stages = NULL, DA = 'stage_lineage', keepX = NULL, save = FALSE)
saveRDS(mmspls, file = 'savedata/MultiModalSparsePLSDA-All.rds')
