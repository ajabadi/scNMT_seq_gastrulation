## ---- function to set up libPaths()
setup_libPaths <- function(username = "alabadi", libpath = "~/R_libs") {
  if (Sys.info()["user"] != username) {
    .libPaths(libpath)
  }
  invisible(NULL)
}
list(nipals_ncomp = 20L, rerun = TRUE)
## ----  main code
setup_libPaths()

library(here)
# list.files(here("notebook"))
rmd_file <- here("notebook/single-omic/nipals.Rmd")
## source utilities
sapply(list.files(here('src', 'utils'), pattern = "*.R", full.names = TRUE), source)
## check/install/load packages
pkgs <- get_required_pkgs(rmd_file)
installer(pkgs = pkgs, hpc_libPaths = "~/R_libs", update_all = TRUE, reinstall = NULL)
load_pkgs(pkgs)
## render
rmarkdown::render(rmd_file, output_format = "html_document", params = list(nipals_ncomp = 20L, subset = c(0, 1000), rerun = TRUE, saveRDS = TRUE, 
                                                                           repeat_runs = 1L, props_NA = c(0.1, 0.25, 0.4)))