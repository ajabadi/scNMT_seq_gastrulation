library(MultiAssayExperiment)
library(here)
library(ggplot2)
library(knitr)

mmspls_E7.5 <- multimodal_analysis_wrapper(mae = gastru.mae, study_assays = NULL, ncomp = 3, scale = FALSE, design = 'null', lineages = NULL, stages = "E7.5", DA = NULL, keepX = NULL, save = FALSE)
saveRDS(mmspls, file = 'savedata/MultiModalSparsePLS-E7.5.rds')