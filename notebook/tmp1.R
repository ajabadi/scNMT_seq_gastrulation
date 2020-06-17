library(MultiAssayExperiment)
library(here)
library(ggplot2)
library(knitr)


mmsplsda_E7.5 <- multimodal_analysis_wrapper(mae = gastru.mae, study_assays = NULL, ncomp = 3, scale = FALSE, design = 'null', lineages = NULL, stages = "E7.5", DA = 'stage_lineage', keepX = NULL, save = FALSE)
saveRDS(mmspls, file = 'savedata/MultiModalSparsePLSDA-E7.5.rds')