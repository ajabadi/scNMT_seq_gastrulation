params <-
  list(on_hpc = c(user = TRUE), on_mac = c(user = FALSE))

params <-
list(on_hpc = c(user = FALSE), on_mac = c(user = TRUE))

## ---- warning=FALSE, message=FALSE-------------------------------------------------------------------------------------------------
setup_libPaths()

library(MultiAssayExperiment)
library(MOFA)
library(mixOmics)
library(reticulate)

# Using a specific python binary

which_py <- ifelse(params$on_hpc, "/usr/local/easybuild/software/Python/2.7.13-GCC-6.2.0-bare/bin/python", "/usr/local/bin/python")
use_python(which_py, required = TRUE)


## ---- echo=FALSE-------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(eval = TRUE, cache = FALSE)


## ----------------------------------------------------------------------------------------------------------------------------------
library(here)
gastru.mae <- readRDS(here('output/scnmtseq_gastrulation_mae_826-cells_orderedFeatures.rds'))  # MAE experiment
gastru.mae


## ----------------------------------------------------------------------------------------------------------------------------------
study_assays <- c("rna", "met_promoter", "acc_promoter", "met_p300", "acc_p300")


## ----------------------------------------------------------------------------------------------------------------------------------
study_mae <- gastru.mae[,,study_assays]


## ---- eval=params$on_mac, echo=params$on_mac---------------------------------------------------------------------------------------
## for subsetting the P x N data on local computer
subset_pn <- function(mat, n=50, p=100) {
  n <- min(n, dim(mat)[2])
  p <- min(p, dim(mat)[1])
  mat[seq_len(p), seq_len(n)]
}

subset_mae <- function(mae, samples=50, max_feat=150) {
  for (i in seq_along(experiments(mae))) {
    feats <- min(dim( mae[[i]])[1], max_feat)
    mae[[i]] <-  mae[[i]][seq_len(feats),]
  }
  mae[,seq_len(samples),]
}

if (params$on_mac) {
  gastru.mae <- subset_mae(gastru.mae)
}



## ----------------------------------------------------------------------------------------------------------------------------------
runtimes <- list()
MOFAobjects <- list()


## ----------------------------------------------------------------------------------------------------------------------------------
MOFAobject <- createMOFAobject(study_mae)
MOFAobject


## ----------------------------------------------------------------------------------------------------------------------------------
plotDataOverview(MOFAobject)


## ----------------------------------------------------------------------------------------------------------------------------------
DataOptions <- getDefaultDataOptions()
DataOptions$scaleViews <- TRUE
DataOptions


## ----------------------------------------------------------------------------------------------------------------------------------
ModelOptions <- getDefaultModelOptions(MOFAobject)
ModelOptions$numFactors <- 2
ModelOptions


## ----------------------------------------------------------------------------------------------------------------------------------
TrainOptions <- getDefaultTrainOptions()

# Automatically drop factors that explain less than 2% of variance in all omics
TrainOptions$DropFactorThreshold <- 0

TrainOptions$seed <- 2017

TrainOptions


## ----------------------------------------------------------------------------------------------------------------------------------
MOFAobject <- prepareMOFA(
  MOFAobject, 
  DataOptions = DataOptions,
  ModelOptions = ModelOptions,
  TrainOptions = TrainOptions
)


## ----------------------------------------------------------------------------------------------------------------------------------
# MOFAobject <- regressCovariates(
#   object = MOFAobject,
#   views = c("Drugs","Methylation","mRNA"),
#   covariates = MOFAobject@InputData$Gender
# )


## ----MOFA-run, include=FALSE, echo=TRUE--------------------------------------------------------------------------------------------
runtimes$gastru$MOFA <-  system.time({
  MOFAobject <- runMOFA(MOFAobject)
})



## ----------------------------------------------------------------------------------------------------------------------------------
if (params$on_hpc) {
  saveRDS2(MOFAobject, file = here("output/MOFAobject.rds"), suffix = "auto", log = paste0(study_assays, collapse = ", "))
}


## ----------------------------------------------------------------------------------------------------------------------------------
X <- lapply(as.list(experiments(study_mae)), t)


## ----------------------------------------------------------------------------------------------------------------------------------
create_keepX <- function(mae, keepX) {
  rep(list(keepX), length(mae)) %>% set_names(names(mae))
}


## ----------------------------------------------------------------------------------------------------------------------------------
design <- 0.5*(1-diag(length(study_mae)))
keepX <- c(5,10,20,50)
runtimes$gastru$tuno_diablo <- system.time({
  tune_diablo_res <-
    tune.block.splsda(
      X = X,
      Y = study_mae$lineage10x_2,
      ncomp = 2,
      test.keepX = create_keepX(study_mae, keepX = keepX),
      design = design,
      cpus = parallel::detectCores()
    )
})



## ----------------------------------------------------------------------------------------------------------------------------------
if (params$on_hpc) {
  saveRDS2(tune_diablo_res, file = here("output/tune_diablo_res.rds"), suffix = "auto", log = paste0(study_assays, collapse = ", "))
}


## ----------------------------------------------------------------------------------------------------------------------------------
table(runtimes$gastru)


## ----------------------------------------------------------------------------------------------------------------------------------
sessionInfo()


## ---- eval=FALSE-------------------------------------------------------------------------------------------------------------------
## knitr::purl("mofa.Rmd")

