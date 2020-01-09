# met_dist_mats <- readRDS("~/Projects/analysis/R/multiOmics/scNMT_seq_gatrulation/spartan/output/met_dist_mats.rds")
# acc_dist_mats <- readRDS("~/Projects/analysis/R/multiOmics/scNMT_seq_gatrulation/spartan/output/acc_dist_mats.rds")

rm(list = ls())
PATH <- "spartan/output"
# PATH="output"
rds <- list.files(path = PATH, pattern = "*.rds", full.names = TRUE)
## basenames
rds_bn <- gsub(x = basename(rds), pattern = ".rds", "")

mapply(FUN = function(x,y) assign(x = x, value = readRDS(y), envir = .GlobalEnv), rds_bn, rds)
