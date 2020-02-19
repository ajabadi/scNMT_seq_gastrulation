###### inputs
source("src/utils/hpc_utils.R")

rmd_file <- "notebook/single-omic/single-omic-visualisations.Rmd"

setup_hpc_for_file <- function(rmd_file) {
  setup_libPaths()
  .libPaths()
  pkgs <- get_required_pkgs(rmd_file)
  pkgs
  pkgs <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
  installer(pkgs = pkgs, hpc_libPaths = "~/R_libs", update_all = TRUE, reinstall = NULL)
}

setup_hpc_for_file(rmd_file = rmd_file)
