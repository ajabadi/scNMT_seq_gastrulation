
check_files_exist <- function(io_list){
  ## check that a list of directories and files exist
  
  
  ## flatten the list
  io_flat <- as.list(unlist(io_list, recursive=TRUE))
  ## find non-exisiting ones
  non_exist <- io_flat  %>% .[!unlist(lapply(., file.exists))]
  if (length(non_exist))
    ## if non-NULL, throw condition
    stop(paste0('could not find the following file(s): \n',
                paste0(non_exist, collapse = " \n")))
}

symlink <- function(src, pattern, repl) {
  ## sym link src to a dest resulting from pattern 2 repl replacement
  ## this is to avoid copying all datasets
  dest <- sub(pattern,repl, src)
  if (!file.exists(dest)) {
    ## create the parent first
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    system(sprintf("ln -s %s %s", src, dest))
  }

}

cp <- function(src, pattern, repl) {
  ## copy src to a dest resulting from pattern 2 repl replacement
  ## this is to avoid copying all datasets
  dest <- sub(pattern,repl, src)
  if (!file.exists(dest)) {
    ## create the parent first
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    file.copy(src, dest)
  }
  
}