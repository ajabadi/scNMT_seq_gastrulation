## ----------- Set up .libPaths on HPC ----------- ##

#' Set up .libPaths on HPC
#'
#' @param username local user name
#' @param libpath R libpath created on home of HPC
setup_libPaths <- function(username = "alabadi", libpath = "~/R_libs") {
    if (Sys.info()["user"] != username) {
        .libPaths(libpath)
    }
    invisible(NULL)
}

## ----------- Find required packages for a script ----------- ##
get_required_pkgs <- function(file_name) {
    setup_libPaths()
    if (!requireNamespace('qdapRegex', quietly = TRUE)) install.packages('qdapRegex')
    library(qdapRegex)
    ## read the script
    the.script <- readLines(file_name)
    ## look for the libraries used
    req.pkgs <- rm_between(the.script, c('library(', 'require('), c(')', ')'), extract=TRUE)
    ## get a vector of required packages not already installed
    req.pkgs <- unique(req.pkgs)
    req.pkgs <- unlist(req.pkgs[!is.na(req.pkgs)])
    return(req.pkgs)
}

## ----------- Install Required Package on HPC ----------- ##
installer <- function(pkgs=c( "data.table", "magrittr", "mixOmics", "scater", "purrr", "uwot"),
                      update_all=FALSE, 
                      hpc_libPaths = "~/R_libs",
                      reinstall = NULL) {
    setup_libPaths()
    
    ## to update old spartan packages in my R_libs
    hpc_pkgs <- installed.packages()[,c("Package") ]
    missing_pkgs <- setdiff(pkgs, hpc_pkgs)
    ## if a package on HPC is not updated
    missing_pkgs <- unique(c(missing_pkgs, reinstall))
    
    options("install.packages.compile.from.source" = "no")
    
    if (!requireNamespace("BiocManager"))
        install.packages('BiocManager', quiet = TRUE)
    
    for (i in missing_pkgs) {
        cat(sprintf("\nTrying to install: %s \n", i))
        msg <- "was successfully installed"
        tryCatch({BiocManager::install(i, update = update_all, ask = FALSE, site_repository = BiocManager::repositories()[1])}, warning = function(e) {
            msg <<- "failed to install"
        }
        )
        cat(sprintf("\npackage %s %s \n", i, msg))
    }
    
    failed <- setdiff(missing_pkgs,  installed.packages()[,c("Package") ])
    succeeded <- setdiff(missing_pkgs, failed)
    installed <- setdiff(pkgs, missing_pkgs)

    if (length(failed) == 0 )
        failed <- 'NONE'
    if (length(succeeded) == 0 )
        succeeded <- 'NONE'
    if (length(installed) == 0 )
        installed <- 'NONE'
    
    message(sprintf("\nunable to install:\n"))
    cat(paste0(failed, collapse = ", "))
    message(sprintf("\nsuccessfully installed:\n"))
    cat(paste0(succeeded, collapse = ", "))
    message(sprintf("\nalready installed:\n"))
    cat(paste0(installed, collapse = ", "))
    cat("\n")
    
}

## ----------- quietly load packages ----------- ##
load_pkgs <- function(pkgs) {
    setup_libPaths()
    invisible(sapply(pkgs, function(x) suppressMessages(library(x, character.only = TRUE))))
}

## ----------- single quote ----------- ##

sQuote2 <- function(foo=""){
    sprintf("'%s'", foo)
}
## ----------- rsync from R ----------- ##

prj_dir <- "Projects/analysis/R/multiOmics/scNMT_seq_gatrulation" ## no trailing slash
mac_dir <- file.path("/Users/alabadi", prj_dir)
hpc_dir <- file.path("/data/cephfs/punim0613/AL", prj_dir) ## do mkdir --parents THIS on hpc
mac_rsync_in <- "hpc/rsync" ## rsync back into this from hpc to avoid inadvertant deletions - relative to mac_dir - no slashes
## included files when  mac -> hpc
rsync_out_include = file.path(mac_dir, "hpc/include-out.txt")
rsync_out_include = file.path(mac_dir, "hpc/include-in.txt")
## excluded files when  mac -> hpc
rsync_out_exclude = file.path(mac_dir, "hpc/exclude-out.txt")
rsync_in_exclude = file.path(mac_dir, "hpc/exclude-in.txt")

#' rsync local -> HPC from within R
#' 
#' Having opened the project and logged into HPC, you can rsync from within R
#'
#' @param del pattern/files to delete
#' @param ignore pattern/files to ignore
rsync_out_exc <- function(del=FALSE, ignore=NULL, inc = NULL) {
    
    ignore <- ifelse(is.null(ignore), "", paste0(" --exclude=", sQuote2(ignore)))
    delit <- ifelse(del, " --delete ", " ")
    incs <- ifelse(is.null(inc), "", paste0(" --include=", sQuote2(inc), collapse = " "))
    cmd <- sprintf("rsync --recursive --verbose --update --progress --compress %s %s --exclude-from=%s %s %s/ ajabadi@spartan.hpc.unimelb.edu.au:%s/", delit, incs , sQuote2(rsync_out_exclude), ignore, mac_dir, hpc_dir)
    system(cmd)
}

#' rsync local <- HPC from within R
#' 
#' Having opened the project and logged into HPC, you can rsync from within R
#'
#' @param del pattern/files to delete
#' @param ignore pattern/files to ignore
rsync_in_exc <- function(del=FALSE, ignore=NULL, inc = NULL) {
    
    ignore <- ifelse(is.null(ignore), "", paste0(" --exclude=", sQuote2(ignore)))
    delit <- ifelse(del, " --delete ", " ")
    incs <- ifelse(is.null(inc), "", paste0(" --include=", sQuote2(inc), collapse = " "))
    cmd <- sprintf("rsync --recursive --verbose --progress --compress %s %s --exclude-from=%s %s ajabadi@spartan.hpc.unimelb.edu.au:%s/ %s/", delit, incs , sQuote2(rsync_in_exclude), ignore, hpc_dir, mac_dir)
    system(cmd)
}

#' rsync slurm.out files local <- HPC 
#' 
#' Having opened the project and logged into HPC, you can rsync from within R
#'
#' @param del whether to delete the current slurm files
rsync_in_slurmouts_only <- function(del=FALSE) {
    
    delit <- ifelse(del, " --delete ", " ")
    cmd <- sprintf("rsync -r -v%s  --include='*.out' --exclude='*' ajabadi@spartan.hpc.unimelb.edu.au:%s/ %s/spartan/", delit, hpc_dir, mac_dir)
    # cmd <- sprintf("rsync -r -v%s--files-from='./hpc/include.txt' ajabadi@spartan.hpc.unimelb.edu.au:%s %s/spartan", delit, mac_dir, hpc_dir)
    system(cmd)
}

rsync_out_inc <- function(del=FALSE, ignore=NULL) {
    
    ignore <- ifelse(is.null(ignore), "", paste0(" --exclude=", sQuote(ignore)))
    delit <- ifelse(del, " --delete ", " ")
    cmd <- sprintf("rsync --recursive --verbose --update --progress --compress %s --files-from=%s %s %s/ ajabadi@spartan.hpc.unimelb.edu.au:%s/", delit, sQuote2(rsync_out_include), ignore, mac_dir, hpc_dir)
    system(cmd)
}




## create spartan folder first

rsync_in_inc_exc <- function(del=FALSE, inc=NULL, ignore=c("*.gz", "*.rds", "notebook", "__*")) {
    
    ignore <- ifelse(is.null(ignore), "", paste0(" --exclude=", sQuote2(ignore)))
    include <- ifelse(is.null(inc), "", paste0(" --include=", sQuote2(inc)))
    delit <- ifelse(del, " --delete ", " ")
    cmd <- sprintf("rsync --recursive --verbose --progress --compress %s %s %s ajabadi@spartan.hpc.unimelb.edu.au:%s/ %s/hpc/rsync/", delit, include, ignore, hpc_dir, mac_dir)
    system(cmd)
}


rsync_in_inc <- function(del=FALSE, ignore=c("*.gz")) {
    
    ignore <- ifelse(is.null(ignore), "", paste0(" --exclude=", sQuote2(ignore)))
    delit <- ifelse(del, " --delete ", " ")
    cmd <- sprintf("rsync -r -v%s --files-from='./hpc/include-from-hpc.txt'%s ajabadi@spartan.hpc.unimelb.edu.au:%s/ %s/hpc/rsync/", delit, ignore, hpc_dir, mac_dir)
    # cmd <- sprintf("rsync -r -v%s--files-from='./hpc/include.txt' ajabadi@spartan.hpc.unimelb.edu.au:%s %s/spartan", delit, mac_dir, hpc_dir)
    system(cmd)
}


rsync_rm_slurmouts_on_hpc <- function() {
    file.remove(list.files("spartan/", pattern = "*.out", full.names = TRUE))
    delit <- " --delete "
    cmd <- sprintf("rsync -r -v%s  --include='*.out' --exclude='*' %s/spartan/ ajabadi@spartan.hpc.unimelb.edu.au:%s/", delit,  mac_dir, hpc_dir)
    # cmd <- sprintf("rsync -r -v%s--files-from='./hpc/include.txt' ajabadi@spartan.hpc.unimelb.edu.au:%s %s/spartan", delit, mac_dir, hpc_dir)
    system(cmd)
}

## ----------- saveRDS with backup ----------- ##

saveRDS2 <- function(object, file, force=FALSE, suffix=NULL, log=NULL, ...){
    
    parent_base_ext <- function(filePath){
        
        file_dir <- dirname(filePath)
        file_name <- tools::file_path_sans_ext(basename(filePath))
        file_ext <- tools::file_ext(filePath)
        return(list(parent=file_dir, base=file_name, ext=file_ext))
    }
    
    suffixer <- function(file, suffix){
        pl <- parent_base_ext(file) ## path list
        pl$sfx <- suffix
        with(pl, sprintf('%s/%s%s.%s', parent, base,sfx,ext))
    }
    
    file_sfx <- file
    ## if file aleardy exists and no forcing asked  ---------------------------------------
    if(file.exists(file) & !force){
        ## esnure suffix is valid
        if(!is(try(suffix), 'character')){
            stop(sprintf('You sure you want to overwrite %s?. Use force=TRUE to overwrite.', sQuote2(file)), call. = FALSE)
        } else {
            if(suffix=='auto'){ ## if auto suffix asked, number files
                file_name <- tools::file_path_sans_ext(basename(file))
                file_dir <- dirname(file)
                nf <- length(list.files(path = file_dir, pattern = file_name))
                file_sfx <- suffixer(file=file, suffix = sprintf('(%i)', nf))
            } else{
                file_sfx <- suffixer(file=file, suffix = suffix)
            }
        }
    }
    ## if logging asked  ---------------------------------------
    if(try(is(log,'character'))){
        ## check that __log dir exists in parent
        logfile <- sprintf('%s/__log/%s.log', dirname(file), basename(file))
        ## if logfile does not exist
        if(!file.exists(logfile)){
            ## check directory exists first
            if(!dir.exists(dirname(logfile))){
                dir.create(dirname(logfile))
            }
            file.create(logfile)
            cat('initialise log on', date(),'\n\n', file = logfile)
        }
        cat(sprintf('%s: %s \n\n', basename(file_sfx), log), file = logfile, append = TRUE)
    }
    
    saveRDS(object=object, file=file_sfx,...)
}

## ----------- paste parameters used in Rmd yaml for logging ----------- ##
paste_params <- function(params) {
    paste0(paste0(names(params), ": "), paste0(params), collapse = " ; ")
}