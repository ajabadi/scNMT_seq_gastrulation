
## ------------------------------------------------------------------------ ##
calc_na_prop <- function(mat) {
  mat <- as.matrix(mat)
  NA_prop <- sum(is.na(mat)) / prod(dim(mat))
  round(NA_prop, 2)
}
## ------------------------------------------------------------------------ ##
rand_na <- function(mat=matrix(1:100, ncol = 20), prop=0.3){
  mat <- as.matrix(mat)
  
  if (any(is.na(mat))) {
    mat_na <- calc_na_prop(mat = mat)
    message(sprintf("Matrix already contains %s%% NAs\n", 100*mat_na))
    if (prop > 0 & mat_na > prop)
      warning(sprintf("matrix already contains %s%% NAs, asked for: %s%%. Setting some NAs to mean.", 100*mat_na , 100*prop))
  }

  if (prop == 0) {
    out <- mat
  } else {
    repeat_rand <- TRUE
    seed_no <- 42
    i <- 1
    while (repeat_rand) {
      
      set.seed(seed_no)
      ## calculate the total number of NAs
      total_na <- floor(prop*prod(dim(mat)))
      vec <- as.vector(mat)
      already_na <- is.na(vec)
      total_na <- total_na - sum(already_na)
      vec[sample(seq_along(vec)[!already_na], size = total_na, replace = FALSE)] <- NA
      out <- matrix(vec, ncol = ncol(mat), dimnames = dimnames(mat))
      ## make sure no column/row has more than
      repeat_rand <- (max(colSums(is.na(out))) > 0.9*dim(out)[1]) | (max(rowSums(is.na(out))) > 0.9*dim(out)[2])
      seed_no <- seed_no  +  i
      i <- i  +  1
      if (i > 2000) {
        warning("Failed to create suitable matrix with disperse NAs", call. = FALSE)
        break
      }
      
    }
  }
  
  
  out
}
## ------------------------------------------------------------------------ ##
subset_please <- function(dataset, weights = NULL, np_subset = c(0, 0), pheno = NULL) {
  dataset <- as.matrix(dataset)
  set.seed(21)
  nzv_cols <- colVars(dataset, na.rm = TRUE) < 1e-3
  dataset <- dataset[,!nzv_cols]
  np_data <- dim(dataset)
  np_subset <- mapply(x = np_data, y = np_subset, FUN = function(x, y){
    if (y == 0 | y > x) 
      y <- x
    y
  })
  rows <- sample(x = np_data[1], size = np_subset[1], replace = FALSE)

  cols <- order(colVars(dataset[rows,], na.rm = TRUE), decreasing = TRUE) <= np_subset[2]

  nzv_cols <- colVars(dataset[rows, cols], na.rm = TRUE) < 1e-3
  dataset <- dataset[rows,cols][,!nzv_cols]
  
  out <- list(dataset = dataset, pheno = pheno[rows])
  if (!is.null(weights)) {
    out$weights <- weights[rows, cols]
  }
  return(out)
}
## ------------------------------------------------------------------------ ##

## ------------------------------------------------------------------------ ##
## run mixOmics:: and niapls::nipals on dataset with different NA proportions and records run times
benchmark_nipals <-
  function(dataset,
           weights = NULL,
           props = c(0.1, 0.25, 0.4),
           ncomp = 2,
           center = TRUE,
           scale = TRUE,
           gramschmidt = FALSE,
           reconst = FALSE,
           repeat_runs = 1,
           subset = c(0, 0),
           pheno=NULL
  ){
    
    max.iter <- 500
    tol <- 1e-9
    
    benchmark_nipals_helper <- function(dataset_na, repeat_runs, weights=NULL) {
      dataset_na <- dataset_na[,colVars(dataset_na, na.rm = TRUE) > 0]
      nipals_nipals <- nipals::nipals(x = dataset_na, ncomp = ncomp, center = FALSE, scale = FALSE, maxiter = max.iter, tol = tol, fitted = reconst, gramschmidt = FALSE, verbose = FALSE, force.na = TRUE)
      ## runtimes
      mb.res_nipals <- microbenchmark(
        "nipals::nipals" = nipals::nipals(x = dataset_na, ncomp = ncomp, center = FALSE, scale = FALSE, maxiter = max.iter, tol = tol, fitted = reconst, gramschmidt = FALSE, verbose = TRUE),
        times = repeat_runs, unit = "s")
      
      mb.class <- class(mb.res_nipals)
      
      if (is.null(weights)) {
        mixOmics_nipals <- mixOmics::nipals(X = dataset_na, ncomp = ncomp, reconst = reconst, max.iter = max.iter, tol = tol)
        
        mb.res_mixOmics <- microbenchmark(
          "mixOmics::nipals" = mixOmics::nipals(dataset_na, ncomp = ncomp, reconst = reconst, max.iter = max.iter, tol = tol),
          times = repeat_runs, unit = "s")
        mb.res <- Reduce(f = rbind, x = list(as.data.frame(mb.res_mixOmics), as.data.frame(mb.res_nipals)))
        class(mb.res) <- mb.class
      } else {
        nipals_empca <- nipals::empca(x = dataset, w = weights, ncomp = ncomp, center = center, scale = scale, maxiter = max.iter, tol = tol, fitted = reconst, gramschmidt = TRUE)
        
        mb.res_empca <- microbenchmark(
          "nipals::empca" = nipals::empca(x = dataset, w = weights, ncomp = ncomp, center = center, scale = scale, maxiter = max.iter, tol = tol, fitted = reconst, gramschmidt = gramschmidt),
          times = repeat_runs, unit = "s")
        
        mb.res <- Reduce(f = rbind, x = list(as.data.frame(mb.res_empca), as.data.frame(mb.res_nipals)))
        class(mb.res) <- mb.class
      }
      
      out <- list(
        dataset = dataset_na,
        microbenchmark = mb.res,
        nipals_nipals = nipals_nipals
      )
      if (is.null(weights)) {
        out$mixOmics_nipals <- mixOmics_nipals
      } else {
        out$nipals_empca <- nipals_empca
      }
      
      ## output
      return(out)
    }
    
    ## ensure subset n, p) is valid
    subset_res <- subset_please(dataset = dataset, weights = weights, np_subset = subset, pheno = pheno)
    
    dataset <- subset_res$dataset
    weights <- subset_res$weights
    pheno <- subset_res$pheno
    
    nipals_benchmark_results <- list()
    
    for (prop_i in props) {
      dataset_na <- rand_na(dataset, prop = prop_i)
      dataset_na <- base::scale(x = dataset_na, scale= scale, center = center)
      nipals_benchmark_results[[paste0("prop_NA_", prop_i)]] <- tryCatch(benchmark_nipals_helper(dataset_na = dataset_na, repeat_runs = repeat_runs, weights = weights), error = function(e) e)
    }
    if (is.null(weights)) {
      nipals_benchmark_results[["pca"]] <- mixOmics::pca(dataset, ncomp = ncomp, center = center, scale = scale, max.iter = max.iter, tol = tol)
    }
    
    nipals_benchmark_results[["pheno"]] <- pheno
    return(nipals_benchmark_results)
  }
## ------------------------------------------------------------------------ ##
plot_helper <- function(df, comps = c(1,2), title = NULL) {
  axes <- paste0("PC", comps)
  p <- ggplot(df)  + labs(x = axes[1], y = axes[2], title = title)
  
  if ("pheno" %in% colnames(df)) {
    p <- p  +  geom_point(aes_string(x = axes[1], y = axes[2], col = "pheno"))  + 
      guides(col = guide_legend(title = "Phenotype"))
    if (is.numeric(df$pheno)){
      p <- p  +  scale_color_gradient(low = "darkblue", high = "orange")
    }
  } else {
    p <- p  +  geom_point(aes_string(x = axes[1], y = axes[2]))
  }
  p
}
## ------------------------------------------------------------------------ ##
show_benchmark_res <- function(bm_res, which=1, comps=c(1,2), col=TRUE, empca = FALSE, subtitle = NULL) {
  if (is(bm_res[[which]], "error")) {
    cat("This run has encountered error due to abundance of NAs")
  } 
  else 
  {
    library(microbenchmark)
    show_microbench <- function(x) microbenchmark:::print.microbenchmark(x, unit = "s", signif = 2)
    if (is.null(bm_res$pheno)) {
      if (isTRUE(col))
        cat("No $pheno output for coloring")
      col <- FALSE
    }
    res <- list()
    res_i <- bm_res[[which]]
    dataset <- res_i[["dataset"]]
    
    res$benchmark <- show_microbench(res_i[["microbenchmark"]])
    
    if (empca) {
      NA_prop <- sum(is.na(dataset))/prod(dim(dataset))
      empca_res <- res_i[["nipals_empca"]]
      empca_df <- data.frame(empca_res$scores)
      if (isTRUE(col)) {
        empca_df$pheno <- bm_res$pheno
      }
     
      title <- sprintf("nipals::empca - prop_NA: %s", round(NA_prop, 2))
      res$empca <- plot_helper(empca_df, comps=comps, title = title)
      
      nipals_res <- res_i[["nipals_nipals"]]
      nipals_df <- data.frame(nipals_res$scores)
      if (isTRUE(col)) {
        nipals_df$pheno <- bm_res$pheno
      }
      title <- sprintf("nipals::nipals - prop_NA: %s", round(NA_prop, 2))
      res$nipals <- plot_helper(nipals_df, comps=comps, title = title)
    } else {
      mixo_res <- res_i[["mixOmics_nipals"]]
      pca_df <- bm_res[["pca"]]$variates$X %>% data.frame()
      if (isTRUE(col)) {
        pca_df$pheno <- bm_res$pheno
      }
      res$pca <-  plot_helper(pca_df, comps=comps, title = "PCA")
      
      mixo_df <- mixo_res$t %>% data.frame() %>% set_colnames(paste0("PC", seq_len(dim(mixo_res$p)[2]))) %>% set_rownames(rownames(dataset))
      
      if (isTRUE(col)) {
        mixo_df$pheno <- bm_res$pheno
      }
      title <- "mixOmics::nipals"
      if (!grepl("0$", x = title)) title <- sprintf("%s - %s", title, names(bm_res)[which])
      res$mixOmics <- plot_helper(mixo_df, comps=comps,  title =  title)
      
      nipals_res <- res_i[["nipals_nipals"]]
      nipals_df <- data.frame(nipals_res$scores)
      if (isTRUE(col)) {
        nipals_df$pheno <- bm_res$pheno
      }
      title <- "nipals::nipals"
      if (!grepl("0$", x = title)) title <- sprintf("%s - %s", title, names(bm_res)[which])
      if (grepl("0$", x = title)) title <- NULL
      res$nipals <- plot_helper(nipals_df, comps=comps, title = title)
    }
    
    if( !is.null(subtitle)) {
      res[-1] <- lapply(res[-1], function(gg) {
        gg  +  labs(subtitle = subtitle)
      })
    }
    
    res
  }
  
  
}

benchmark_nipals_across_datasets <- function(datasets, props_NA, phenos, ...) {
  mapply(x = datasets, y = phenos, FUN = function(x, y,  ...) {
    bm_res <- benchmark_nipals(dataset = x, props = props_NA, pheno = y, ...)
  })
}
## ------------------------------------------------------------------------ ##
get_datasets_runtimes <- function(propNA, bm_res) {
  theList <- lapply(bm_res, function(x){
    summary(x[[propNA]][["microbenchmark"]])[,c("expr","mean")]
  })
  rbindListWithNames(theList)
}

get_all_runtimes <- function(propNA_vec, bm_res) {
  theList <- lapply(named_list(propNA_vec), function(x) {
    get_datasets_runtimes(x, bm_res = bm_res)
  })
  rbindListWithNames(theList, new_col = "prop_NA")
}
## ------------------------------------------------------------------------ ##
ggplot_all_runtimes <- function(df, ...) {
  ggplot(df, aes(x = factor(dataset), y = mean))  +  geom_point(aes(col = expr), size=3)  + facet_grid(.~prop_NA)  +   theme_bw()  + 
    guides(col = guide_legend(title = "function"), shape = guide_legend(title = "NA proportion"))  + 
    labs(...)
}
