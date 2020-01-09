params <-
list(on_hpc = TRUE)

## ----setup, include=FALSE------------------------------------------------
## knitr::opts_chunk$set(echo = TRUE, eval=FALSE)


## ---- eval=FALSE---------------------------------------------------------
## ## link data folder
## # system("ln -s /Users/alabadi/Projects/analysis/R/multiOmics/Ricard/scNMT_gastrulation/data /Users/alabadi/Projects/analysis/R/multiOmics/scNMT_seq_gatrulation")


## ---- eval=TRUE----------------------------------------------------------
## whether it is being run on mac and not on HPC
# on_hpc <- params$on_hpc
on_hpc <- TRUE
on_mac <- !on_hpc


## ---- eval=TRUE----------------------------------------------------------
#########
## I/O ##
#########
io$data <- io$out <- io <- list()
## Rproject/slurm file home - must be FULL PATH so the extracted R code can run as well
io$home <- ifelse(Sys.info()["user"] == "alabadi",
       "/Users/alabadi/Projects/analysis/R/multiOmics/scNMT_seq_gatrulation",
       "/data/cephfs/punim0613/AL/Projects/analysis/R/multiOmics/scNMT_seq_gatrulation")
## outputs
io$outdir <- file.path(io$home, "output")


## ---- eval=TRUE----------------------------------------------------------
## ----------- source utils ----------- 
## list all R files in utils folder recursively
utls <- list.files(file.path(io$home, "src/utils"), pattern = ".R", full.names = TRUE, recursive = TRUE)
## source all
invisible(sapply(utls, source))


## ---- eval=on_hpc--------------------------------------------------------
installer(pkgs = c("data.table", "magrittr", "mixOmics", "scater", "purrr", "uwot"))


## ---- eval=TRUE----------------------------------------------------------
if (on_hpc) {
  .libPaths("~/R_libs")
}

suppressMessages(library(data.table))
suppressMessages(library(scater))
suppressMessages(library(purrr))
suppressMessages(library(magrittr))
suppressMessages(library(mixOmics))


## ---- eval=TRUE----------------------------------------------------------
#############
## Options ##
#############
opts <- list()
rtimes <- list() ## record run times
# Multiple testing correction options
# opts$threshold_fdr  <- 0.10

opts$met.annos <- c(
  Promoters = "prom_2000_2000",
  Genebody = "genebody",
  CGI = "CGI",
  p300 = "ESC_p300",
  CTCF = "ESC_CTCF",  ##CTCF binding sites.  The primary role of CTCF transcription factors is thought to be in regulating the 3D structure of chromatin.
  DHS = "ESC_DHS" ## In these specific regions of the genome, chromatin has lost its condensed structure, exposing the DNA and making it accessible.
  
)

opts$acc.annos <- c(
 Promoters = "prom_2000_2000",
  Genebody = "genebody",
  CGI = "CGI",
  p300 = "ESC_p300",
  CTCF = "ESC_CTCF",  ##CTCF binding sites.  The primary role of CTCF transcription factors is thought to be in regulating the 3D structure of chromatin.
  DHS = "ESC_DHS" ## In these specific regions of the genome, chromatin has lost its condensed structure, exposing the DNA and making it accessible.
  
)


## ---- eval=TRUE----------------------------------------------------------
## files/directories from data-alias that we'll need
io$data$base   <- file.path(io$home, "data/data-alias/gastrulation") ## data dir
io$data$sample.metadata <- file.path(io$data$base,"sample_metadata.txt")
# io$annos.dir  <- file.path(io$base_dir, "features/filt")
io$data$rna.file   <- file.path(io$data$base, "rna/parsed/SingleCellExperiment.rds")
io$data$met.dir   <- file.path(io$data$base, "met/parsed")
io$data$acc.dir   <- file.path(io$data$base, "acc/parsed")


## ---- eval=TRUE----------------------------------------------------------
## input files
io$data$mets <- lapply(opts$met.annos, function(x) sprintf("%s/%s.tsv.gz", io$data$met.dir, x))
io$data$accs <- lapply(opts$acc.annos, function(x) sprintf("%s/%s.tsv.gz", io$data$acc.dir, x))


## ---- eval=on_mac--------------------------------------------------------
## ## check that all files/dirs are valid
## check_files_exist(io)


## ---- eval=on_mac--------------------------------------------------------
## ## symlink only the needed datasets from full alias dir to ./data
## ## so only them are rsynced to hpc
## ##   ------- methylation
## invisible(lapply(io$data$mets, function(x) {
##   cp(x,
##           pattern = "/data/data-alias/",
##           repl = "/data/")
## }))
## 
## ##   ------- acc
## invisible(lapply(io$data$accs, function(x) {
##   cp(x,
##           pattern = "/data/data-alias/",
##           repl = "/data/")
## }))
## 
## ##   ------- rna
## invisible(lapply(io$data$rna.file, function(x) {
##   cp(x,
##           pattern = "/data/data-alias/",
##           repl = "/data/")
## }))
## 
## ##   ------- metadata
## invisible(lapply(io$data$sample.metadata, function(x) {
##   cp(x,
##           pattern = "/data/data-alias/",
##           repl = "/data/")
## }))


## ---- eval=TRUE----------------------------------------------------------
## now replace all /data-alias/ with /data/ in data directories
io$data <- rapply(io$data, function(x) sub(pattern = "/data/data-alias/", replacement = "/data/", x = x), how = "list")


## ---- eval=TRUE----------------------------------------------------------
## check that all files/dirs are valid
check_files_exist(io)


## ---- eval=TRUE----------------------------------------------------------
sample_metadata <- fread(io$data$sample.metadata)
# sample_metadata %>% head() %>% View()
sample_metadata <- sample_metadata %>% 
  .[pass_metQC==TRUE & pass_rnaQC==TRUE  & pass_accQC==TRUE ] %>% 
  .[,c("sample","id_met","id_acc","id_rna","stage","lineage10x_2")] %>%
  .[,stage_lineage:=as.factor(paste(stage,lineage10x_2,sep="_"))]

length(unique(sample_metadata$sample)) ## sample size


## ------------------------------------------------------------------------
## # Define which stage and lineages to look at
## # opts$stages <- c("E5.5","E6.5","E7.5")
## # opts$lineages <- c(
## #   "Nascent mesoderm",
## #   "Primitive Streak",
## #   "Mixed mesoderm",
## #   "Epiblast",
## #   "Rostral neurectoderm",
## #   "Gut",
## #   "Embryonic endoderm",
## #   "Paraxial mesoderm",
## #   "Notochord",
## #   "Pharyngeal mesoderm",
## #   "Mature mesoderm",
## #   "Intermediate mesoderm",
## #   "Somitic mesoderm",
## #   "Caudal mesoderm",
## #   "Rostral neuroectoderm",
## #   "Caudal neurectoderm",
## #   "Caudal Mesoderm",
## #   "Def. endoderm"
## # )
## # opts$stage_lineage <- as.vector(outer(opts$stages,opts$lineages, paste, sep="_"))
## # opts$stage_lineage <- opts$stage_lineage[!opts$stage_lineage %in% c("E7.5_Epiblast","E7.5_Primitive Streak")]
## 
## # # Define which cells to use
## # tmp <- fread(io$data$sample.metadata) %>%
## #   .[,stage_lineage:=paste(stage,lineage10x,sep="_")] %>%
## #   # .[assay=="scNMT"] %>%
## #   .[stage_lineage%in%opts$stage_lineage]
## # opts$met_cells <- tmp %>% .[pass_metQC==T,id_met]
## # opts$rna_cells <- tmp %>% .[pass_rnaQC==T,id_rna]
## # opts$acc_cells <- tmp %>% .[pass_accQC==T,id_acc]


## ---- eval=TRUE----------------------------------------------------------
# Load scater object
sce <- readRDS(io$data$rna.file)

# rna_passQC <- sample_metadata[pass_rnaQC==TRUE]$id_rna
# Filter cells
sce <- sce[,colnames(sce) %in% sample_metadata$id_rna]


## ---- eval=TRUE----------------------------------------------------------
## rename id_rna for sce
id_map <- sample_metadata[,c("sample", "id_rna")] %>% data.frame(row.names = 2)
colnames(sce) <- id_map[colnames(sce),]


## ------------------------------------------------------------------------
## ## mean gene library sizes
## hist(rowMeans(logcounts(sce))[rowMeans(logcounts(sce))>0], breaks = 100)
## ## mean gene library sizes > exp(0.3)
## # hist(rowMeans(logcounts(sce))[rowMeans(counts(sce))>0.3], breaks = 100)


## ------------------------------------------------------------------------
## ## Mean vs variance plot
## foo <- data.frame(sd=apply(exprs(sce),1,sd), mean=apply(exprs(sce),1,mean))
## ggplot(foo, aes(x=mean, y=sd)) +
##   geom_point() +
##   stat_smooth() +
##   scale_color_manual(values=c("black","red")) +
##   xlab('Mean') + ylab('Standard deviation')
## 


## ---- eval=TRUE----------------------------------------------------------
io$out$pca_scnmt <- file.path(io$outdir, "pca_scnmt.rds")


## ---- eval=on_hpc--------------------------------------------------------
pca_scnmt <- pca(t(logcounts(sce)), ncomp=10)
saveRDS(pca_scnmt, file = io$out$pca_scnmt)


## ------------------------------------------------------------------------
## pca_scnmt <- readRDS(io$out$pca_scnmt)
## # visualise rna
## plot(pca_scnmt)
## 
## plotIndiv(pca_scnmt, group = sce$stage, pch=16, legend = TRUE)
## ## comp 2,3
## plotIndiv(pca_scnmt, comp = c(2,3), group = sce$stage, pch=16, legend = TRUE)
## 
## plotIndiv(pca_scnmt, comp = c(1,3), group = sce$stage, pch=16, legend = TRUE)


## ---- eval=TRUE----------------------------------------------------------
io$out$umap_scnmt <- file.path(io$outdir, "umap_scnmt.rds")


## ---- eval=on_hpc--------------------------------------------------------
library(uwot)
Ncomps = 3
umap_scnmt <- uwot::umap(t(logcounts(sce)), n_components = Ncomps) %>%
  as.data.frame() %>%
  set_colnames(paste0("UMAP_", seq_len(Ncomps)))
saveRDS(umap_scnmt, file = io$out$umap_scnmt )


## ------------------------------------------------------------------------
## umap_scnmt <- readRDS(io$out$umap_scnmt)
## library(ggplot2)
## ggplot(umap_scnmt) + geom_point(aes(UMAP_1, UMAP_2, col = sce$stage))
## ggplot(as.data.frame(umap_scnmt)) + geom_point(aes(UMAP_1, UMAP_3, col = sce$stage))


## ------------------------------------------------------------------------
## ## function to read data and filter
## # io <- list()
## # io$data$met.dir <- "../data/gastrulation/met/parsed"
## # io$data$acc.dir <- "../data/gastrulation/acc/parsed"


## ---- eval=TRUE----------------------------------------------------------
io$out$met_dt_list <- file.path(io$outdir, "met_dt_list.rds")


## ---- eval=on_hpc--------------------------------------------------------
rtimes$met_dt_list <- system.time({
  
  met_dt_list <- mclapply(io$data$mets, function(tsv_name) {
    
    calc_site_stats(filePath = tsv_name,
                    sample_name = "id_met", ## name of sample to output
                    keep_samples = sample_metadata$id_met, ## pass QC samples
                    min_N = 3, ## min number of calls at site
                    min_cov = 400, ## min number of cells detecting
                    alpha = 0.1) ## for lower bound of CI
  }, mc.cores = parallel::detectCores())
  
})

saveRDS(met_dt_list, file = io$out$met_dt_list)


## ---- eval=TRUE----------------------------------------------------------
io$out$met_dist_mats <- file.path(io$outdir, "met_dist_mats.rds")


## ---- eval=FALSE---------------------------------------------------------
## # met_dist_mats <- list()
## #
## # for (dt_name in names(met_dt_list)) {
## #   met_dist_mats[[dt_name]] <- wt_euc_dist(met_dt = met_dt_list[[dt_name]],
## #                                           sample_name = "id_met",
## #                                           n_hvr = 2000)
## #   gc()
## # }
## #
## # saveRDS(met_dist_mats, file = io$out$met_dist_mats)


## ---- eval=on_hpc--------------------------------------------------------
n_hvr <- 1000
for (dt_name in names(met_dt_list)) {
  dist_mat <- wt_euc_dist(met_dt = met_dt_list[[dt_name]], 
                                          sample_name = "id_met", 
                                          n_hvr = n_hvr)
  saveRDS(dist_mat, file = sprintf("%s/met_dist_mat_%s_%s.rds", io$outdir, dt_name, n_hvr))
  rm(dist_mat)
  gc()
}


## ---- eval=TRUE----------------------------------------------------------
io$out$acc_dt_list <- file.path(io$outdir, "acc_dt_list.rds")


## ---- eval=on_hpc--------------------------------------------------------
rtimes$acc_dt_list <- system.time({
  
  acc_dt_list <- lapply(io$data$accs, function(tsv_name) {
    
    calc_site_stats(filePath = tsv_name,
                    sample_name = "id_acc", ## name of sample to output
                    keep_samples = sample_metadata$id_acc, ## pass QC samples
                    min_N = 3, ## min number of calls at site
                    min_cov = 500, ## min number of cells detecting
                    alpha = 0.1) ## for lower bound of CI
  })
  
})

saveRDS(acc_dt_list, file = io$out$acc_dt_list)


## ---- eval=TRUE----------------------------------------------------------
io$out$acc_dist_mats <- file.path(io$outdir, "acc_dist_mats.rds")


## ---- eval=FALSE---------------------------------------------------------
## acc_dist_mats <- list()
## 
## for (dt_name in names(acc_dt_list)) {
##   acc_dist_mats[[dt_name]] <- wt_euc_dist(met_dt = acc_dt_list[[dt_name]],
##                                           sample_name = "id_acc",
##                                           n_hvr = 2000)
##   gc()
## }
## 
## saveRDS(acc_dist_mats, file = io$out$acc_dist_mats)


## ---- eval=on_hpc--------------------------------------------------------
n_hvr <- 1000
for (dt_name in names(acc_dt_list)) {
  dist_mat <- wt_euc_dist(met_dt = acc_dt_list[[dt_name]], 
                                          sample_name = "id_acc", 
                                          n_hvr = n_hvr)
  saveRDS(dist_mat, file = sprintf("%s/acc_dist_mat_%s_%s.rds", io$outdir, dt_name, n_hvr))
  rm(dist_mat)
  gc()
}


## ------------------------------------------------------------------------
## met_dist_mats <- readRDS(io$out$met_dist_mats)


## ------------------------------------------------------------------------
## cmdscale(d = met_dist_mats[[1]]) %>% heatmap()


## ------------------------------------------------------------------------
## rtimes$dcast_met <- system.time({
##   dcast_met <- dcast_dt_list(met_dt_list, row = "id", col = "id_met", value.var = "rate")
## }) ## 8 secs
## 


## ------------------------------------------------------------------------
## lapply(dcast_met, dim)


## ------------------------------------------------------------------------
## ## quick check datasets to ensure features are properly ranked
## mapply(x=dcast_met, y=met_dt_list, FUN = function(x,y) {
##   ## is the first rowname in matrices the same as the max lci in data.table?
##   identical(x %>% rownames() %>% .[1],
##             y %>% .[which.max(.$lci)] %>% .$id)
## })


## ------------------------------------------------------------------------
## acc_dist_accs <- readRDS(io$out$acc_dist_mats)


## ------------------------------------------------------------------------
## cmdscale(d = acc_dist_mats[[1]]) %>% heatmap()


## ------------------------------------------------------------------------
## rtimes$dcast_acc <- system.time({
##   dcast_acc <- dcast_dt_list(acc_dt_list, row = "id", col = "id_acc", value.var = "rate")
## })
## 


## ------------------------------------------------------------------------
## lapply(dcast_acc, dim)


## ------------------------------------------------------------------------
## ## quick check datasets to ensure features are properly ranked
## mapply(x=dcast_acc, y=acc_dt_list, FUN = function(x,y) {
##   ## is the first rowname in matrices the same as the max lci in data.table?
##   identical(x %>% rownames() %>% .[1],
##             y %>% .[which.max(.$lci)] %>% .$id)
## })


## ------------------------------------------------------------------------
## library(MultiAssayExperiment)
## 
## exper.list <- ExperimentList(rna = logcounts(sce),
## 
##                            met_genebody = dcast_met$Genebody,
##                            met_promoter = dcast_met$Promoters,
##                            met_cgi = dcast_met$CGI,
##                            met_p300 = dcast_met$p300,
##                            met_CTCF = dcast_met$CTCF,
##                            met_DHS = dcast_met$DHS,
## 
##                            acc_genebody = dcast_acc$Genebody,
##                            acc_promoter = dcast_acc$Promoters,
##                            acc_cgi = dcast_acc$CGI,
##                            acc_p300 = dcast_acc$p300,
##                            acc_CTCF = dcast_acc$CTCF,
##                            acc_DHS = dcast_acc$DHS
##                            )


## ------------------------------------------------------------------------
## coldata <- sample_metadata %>% data.frame(row.names = TRUE) %>% .[,c("stage", "lineage10x_2",
## "stage_lineage")] %>% DataFrame()
## scnmtseq_gastrulation_mae <- MultiAssayExperiment(experiments = exper.list, colData = coldata)


## ------------------------------------------------------------------------
## saveRDS(scnmtseq_gastrulation_mae, file= "../output/scnmtseq_gastrulation_mae_826-cells_orderedFeatures.rds")


## ------------------------------------------------------------------------
## # scnmtseq_gastrulation_mae <- readRDS("../output/scnmtseq_gastrulation_mae.rds")


## ------------------------------------------------------------------------
## # sce <- readRDS("../data/gastrulation/rna/parsed/SingleCellExperiment.rds")
## # rowData(sce)


## ------------------------------------------------------------------------
## library(scran)
## fit <- trendVar(sce, use.spikes=FALSE)
## decomp.var <- decomposeVar(sce, fit = fit)
## library(ggplot2)
## ggplot(as.data.frame(decomp.var)) + geom_point(aes(x=mean, y=total), alpha=0.6) + theme_bw()
## ggplot(as.data.frame(decomp.var)) + geom_point(aes(x=mean, y=bio), alpha=0.6) + theme_bw()


## ------------------------------------------------------------------------
## ## add to rowData
## rowData(sce) <- cbind(rowData(sce), decomp.var[,c("bio", "total")])


## ------------------------------------------------------------------------
## library(MultiAssayExperiment)
## 
## exper.list <- ExperimentList(rna = sce,
##                            experiments(scnmtseq_gastrulation_mae)[-1]
##                            )
## 
## MultiAssayExperiment(experiments = exper.list, colData = colData(scnmtseq_gastrulation_mae))
## 
## coldata <- sample_metadata %>% data.frame(row.names = TRUE) %>% .[,c("stage", "lineage10x_2",
## "stage_lineage")] %>% DataFrame()
## scnmtseq_gastrulation_mae <- MultiAssayExperiment(experiments = exper.list, colData = coldata)
## scnmtseq_gastrulation_mae <- MatchedAssayExperiment(mae_scnmt_seq)
## 
## coldata <- sample_metadata %>% data.frame(row.names = TRUE) %>% .[,c("stage", "lineage10x_2",
## "stage_lineage")] %>% DataFrame()
## scnmtseq_gastrulation_mae <- MultiAssayExperiment(experiments = exper.list, colData = coldata)
## scnmtseq_gastrulation_mae <- MatchedAssayExperiment(mae_scnmt_seq)


## ---- eval=FALSE---------------------------------------------------------
## knitr::purl("scnmtseq-gastrulation.Rmd",documentation = 1L, output = "../src/scnmtseq-gastrulation.R")

