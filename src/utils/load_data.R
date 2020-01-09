################
## Define I/O ##
################

# io <- list()
# io$met.dir <- "data/gastrulation/met/parsed"
# io$acc.dir <- "data/gastrulation/acc/parsed"
# io$rna.file <- "data/gastrulation/rna/parsed/SingleCellExperiment.rds"

load_data <- function(rna.id=opts$rna_cells, met.id=opts$met_cells, met.anno=opts$met.annos[1], acc.id=opts$acc.annos[1], acc.anno=met.anno, min.cpg=5, min.gpc=5, io=io) {

  ## function to change ML estimates to MAP estimates for a dt
  ml2map <- function(dt) {
    dt[,rate:=rate/100]
    dt[,rate:=(rate*N+1)/(N+2)]
    dt
  }
  ##########
  ## Load ##
  ##########
  
  # Load DNA methylation data
  met_dt <- fread(cmd=sprintf("zcat < %s/%s.tsv.gz",io$met.dir,met.anno), stringsAsFactors=F, quote="", showProgress=F) %>%
    setnames(c("id_met","id","anno","Nmet","N","rate")) %>% ml2map()
  ## change rate to [0-1] scale
  
  # Load DNA accessibility data
  acc_dt <- fread(cmd=sprintf("zcat < %s/%s.tsv.gz",io$acc.dir,acc.anno), stringsAsFactors=F, quote="", showProgress=F) %>%
    setnames(c("id_acc","id","anno","Nmet","N","rate")) %>% ml2map()
  
  # Load RNA data
  sce <- readRDS(io$rna.file)
  sce <- sce[,colnames(sce)%in%rna.id]
  sce <- sce[rowMeans(logcounts(sce))>=0.3,]
  # rna_dt <- exprs(sce) %>% t %>% as.data.table(keep.rownames = "id_rna") %>% 
  #   melt(id.vars = "id_rna", value.name = "expr", variable.name = "ens_id") %>%
  #   merge(rowData(sce) %>% as.data.frame(row.names = rownames(sce)) %>% tibble::rownames_to_column("ens_id") %>% .[,c("symbol","ens_id")] %>% setnames("symbol","gene"))
  # 
  ############
  ## Filter ##
  ############
  
  # Select ID
  met_dt <- met_dt[id_met%in%met.id] %>% setnames("rate","value")
  acc_dt <- acc_dt[id_acc%in%acc.id] %>% setnames("rate","value")
  # rna_dt <- rna_dt[ens_id%in%rna.id] %>% setnames("expr","value")
  
  # Filter by coverage
  met_dt <- met_dt[N>=min.cpg]
  acc_dt <- acc_dt[N>=min.gpc]
  
  return(list("met"=met_dt, "acc"=acc_dt, "rna"=rna_dt))
}