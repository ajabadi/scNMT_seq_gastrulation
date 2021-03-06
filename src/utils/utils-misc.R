## ----------- set_if_null ----------- 
### from Altools
set_if_null <- function(null_arg=NULL, default_value) {
  return(
    if (is.null(null_arg))
      default_value
    else
      null_arg
  )
}
## ----------- suffix_colnames ----------- 
## add suffix to column names for integration
suffix_colnames <- function(x, indices=1:2){
  for (j in indices) {
    suff <- names(x)[j]
    colnames(x[[j]]) <- paste0(colnames(x[[j]]), '_' ,suff)
  }
  x
}
## ----------- all_identical ----------- 
## check if elements in a list are identical
all_identical <- function(lst) {
  for (j in seq_along(lst[-1])) {
    if (!identical(lst[[j]], lst[[j+1]]))
      stop(sprintf("not identical elements: %s and %s",j , j+1 ), call. = FALSE)
  }
  TRUE
}
## ----------- named_list ----------- 
## create a named list
named_list <- function(char) {
  out <- as.list(char)
  names(out) <- char
  out
}

## ----------- impute_by_mean ----------- 
## impute NAs as mean of columns
impute_by_means <- function(x) {
  apply(x, 2, function(y){
    y[is.na(y)] <- mean(y, na.rm=TRUE)
    y
  })
}
## ----------- ggplot color hue ----------- 
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}

## ----------- histogram of values using ggplot ----------- 
hist_values <- function(values, fill_col = 1, ...) {
  
  df <- data.frame(values = values)
  ggplot(df) + geom_histogram(aes(values), bins = 50, fill = gg_color_hue(30)[fill_col %% 30 + 1]) + theme_bw() +
    labs(...)
}

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

## ----------- create design matrix for MAE ----------- 
create_design <- function(mae, off_diag = 0.5, rna_only = FALSE) {
  

  if (rna_only) {
    design = matrix(0, ncol = length(mae), nrow = length(mae),
                    dimnames = list(names(mae), names(mae)))
    design[1,] <- design[,1] <- off_diag
  } else {
    design = matrix(off_diag, ncol = length(mae), nrow = length(mae),
                    dimnames = list(names(mae), names(mae)))
    
  }
  diag(design) =  0
  
  design 
}

## ----------- create repeated keepX for MAE ----------- 
create_keepX <- function(mae, keepX) {
  rep(list(keepX), length(mae)) %>% set_names(names(mae))
}

## ----------- subset MAE ----------- 
subset_mae <- function(mae, n=NULL, p=NULL, SEED=42) {
  for (i in seq_along(experiments(mae))) {
    feats <- min(dim( mae[[i]])[1], p)
    mae[[i]] <-  mae[[i]][seq_len(feats),]
  }
  n_cells <- dim(mae[[1]])[2]
  set.seed(SEED)
  mae[,sample(x = seq_len(n_cells), size = n),]
}

## for subsetting the P x N data on local computer
subset_pn <- function(mat, n=NULL, p=NULL) {
  `%||%` <- purrr::`%||%`
  p <- p %||% dim(mat)[2]
  n <- n %||% dim(mat)[1]
  
  p <- min(p, dim(mat)[2])
  n <- min(n, dim(mat)[1])
  mat[seq_len(n), seq_len(p)]
}