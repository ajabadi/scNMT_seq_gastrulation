## ---- function to set up libPaths()
setup_libPaths <- function(username = "alabadi", libpath = "~/R_libs") {
  if (Sys.info()["user"] != username) {
    .libPaths(libpath)
  }
  invisible(NULL)
}

## ----  main code
setup_libPaths()

library(here)
# list.files(here("notebook"))
rmd_file <- here("notebook/scnmtseq-gastrulation.Rmd")
## source utilities
sapply(list.files(here('src/utils'), pattern = "*.R", full.names = TRUE), source)
## check/install/load packages
pkgs <- get_required_pkgs(rmd_file)
installer(pkgs = pkgs, hpc_libPaths = "~/R_libs", update_all = TRUE, reinstall = NULL)
load_pkgs(pkgs)
## render
rmarkdown::render(rmd_file, output_format = "html_document", params = list())