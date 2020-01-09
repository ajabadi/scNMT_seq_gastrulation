## ----------- return sce with sorted top HVGs ----------- 
get_hvgs_sce <- function(sce = NULL,
                         n_genes = NULL, ## No. of genes, NULL for all genes ordered by variation
                         use_bio=TRUE ## choose based on biological variance or total variance?
){ 
  ## -- decompose variance
  require(scran)
  fit <- trendVar(x = sce, use.spikes=FALSE)
  sce <- decomposeVar(x = sce, fit = fit, use.spikes=FALSE)
  ## -- order genes based on variance
  gene_order <- ifelse(use_bio, order(-sce$bio), order(-sce$total))
  
  ## -- subset genes (if asked)
  if(is.null(n_genes)) {
    n_genes <- dim(sce)[1]
  }
  
  keep_genes <- gene_order[seq_len(n_genes)]
  
  ## -- subset sce
  sce <- sce[keep_genes,]
  return(sce)
}