# # opts_in <-           " --recursive --verbose --progress --compress --include='*.slurm' --include='notebook/*' --include='src/*' --include='output' --exclude='*' "
# # opts_out <- " --update --recursive --verbose --progress --compress exclude='data' exclude='output' "
# 
# ## ---------- rsync options
# ## mac -> hpc
# opts_out <-         " --recursive  --exclude='data' --exclude='ouput' "
# opts_out
# ## mac <- hpc
# opts_in <- " --recursive --exclude='data' exclude='output' "
# 
# ## ---------- directories
# prj_dir <- "Projects/analysis/R/multiOmics/scNMT_seq_gatrulation" ## no trailing slash
# mac_dir <- file.path("/Users/alabadi", prj_dir)
# hpc_dir <- file.path("/data/cephfs/punim0613/AL", prj_dir) ## do mkdir --parents THIS on hpc
# mac_rsync_in <- "hpc/rsync" ## rsync back into this from hpc to avoid inadvertant deletions - relative to mac_dir - no slashes
# 
# ## ---------- mac -> hpc
# rsync_out <- function(del=FALSE, exc=c("*", ".*", "./*"), inc=c("hpc", "notebook", "*.out")) {
#   
#   options(useFancyQuotes=FALSE)
#   excs <- incs <- ""
#   if (!is.null(exc))
#     excs <- paste0(" --exclude=", sQuote(exc), collapse = " ")
#   if (!is.null(inc))
#     incs <- paste0(" --include=", sQuote(inc), collapse = " ")
#   
#   delit <- ifelse(del, " --delete ", " ")
#   cmd <- sprintf("rsync -rpzv %s %s %s %s/ ajabadi@spartan.hpc.unimelb.edu.au:%s/", delit, incs, excs, mac_dir, hpc_dir)
#   system(cmd)
# }
# 
# ## ---------- mac <- hpc
# rsync_in <- function(del=FALSE, exc=c("output", "data"), inc=NULL) {
#   options(useFancyQuotes=FALSE)
#   mac_hpc_dir <- file.path(mac_dir, mac_rsync_in)
#   excs <- incs <- ""
#   if (!is.null(exc))
#     excs <- paste0(" --exclude=", sQuote(exc), collapse = " ")
#   if (!is.null(exc))
#     incs <- paste0(" --include=", sQuote(inc), collapse = " ")
#   
#   delit <- ifelse(del, " --delete ", " ")
#   cmd <- sprintf("rsync -rpzv %s %s %s ajabadi@spartan.hpc.unimelb.edu.au:%s/ %s/ ", delit, incs, excs, hpc_dir, mac_hpc_dir)
#   system(cmd)
# }
# 
