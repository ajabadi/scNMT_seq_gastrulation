## GSEA using MOFA
## https://raw.githack.com/bioFAM/MOFA2/master/MOFA2/vignettes/GSEA.html

# BiocManager::install("MOFAdata")
# BiocManager::install("MOFA2")
# remotes::install_github("bioFAM/MOFA2/MOFA2", build_opts = c("--no-resave-data --no-build-vignettes"))

library(data.table)
library(purrr)
library(ggplot2)
library(cowplot)
library(MOFAdata)
library(MOFA2)

# C2: curated gene sets from online pathway databases, publications in PubMed, and knowledge of domain experts.
data("MSigDB_v6.0_C2_mouse") 

# C5: extracted from the Gene Ontology data.base
data("MSigDB_v6.0_C5_mouse") 

run_enrichment2 <- function (W, feature.sets, comps = 2, set.statistic = c("mean.diff", 
                                                                                            "rank.sum"), statistical.test = c("parametric", "cor.adj.parametric", 
                                                                                                                              "permutation"), sign = c("all", "positive", "negative"), 
                             min.size = 10, nperm = 1000, cores = 1, p.adj.method = "BH", 
                             alpha = 0.1) 
{
  # if (!is(object, "MOFA")) 
  #   stop("'object' has to be an instance of MOFA")
  # if (!(is(feature.sets, "matrix") & all(feature.sets %in% 
  #                                        c(0, 1)))) 
  #   stop("feature.sets has to be a list or a binary matrix.")
  # view <- .check_and_get_views(object, view)
  # factors <- .check_and_get_factors(object, factors)
  sign <- match.arg(sign)
  set.statistic <- match.arg(set.statistic)
  statistical.test <- match.arg(statistical.test)
  data <- get_data(object, views = view, as.data.frame = FALSE)[[1]]
  if (is(data, "list")) 
    data <- Reduce(cbind, data)
  data <- t(data)
  W <- get_weights(object, views = view, factors = comps, 
                   scale = TRUE)[[1]]
  Z <- get_factors(object, factors = comps)
  if (is(Z, "list")) 
    Z <- Reduce(rbind, Z)
  stopifnot(rownames(data) == rownames(Z))
  idx <- apply(data, 2, function(x) var(x, na.rm = T)) == 0
  if (sum(idx) >= 1) {
    warning(sprintf("%d fetures were removed because they had no variance in the data.\n", 
                    sum(idx)))
    data <- data[, !idx]
    W <- W[!idx, ]
  }
  features <- intersect(colnames(data), colnames(feature.sets))
  if (length(features) == 0) 
    stop("Feature names in feature.sets do not match feature names in model.")
  message(sprintf("Intersecting features names in the model and the gene set annotation results in a total of %d features.", 
                  length(features)))
  data <- data[, features]
  W <- W[features, ]
  feature.sets <- feature.sets[, features]
  feature.sets <- feature.sets[rowSums(feature.sets) >= min.size, 
  ]
  if (sign == "positive") {
    W[W < 0] <- 0
  }
  else if (sign == "negative") {
    W[W > 0] <- 0
    W <- abs(W)
  }
  message("\nRunning feature set Enrichment Analysis with the following options...")
  message(sprintf("View: %s", view))
  message(sprintf("Number of feature sets: %d", nrow(feature.sets)))
  message(sprintf("Set statistic: %s", set.statistic))
  message(sprintf("Statistical test: %s", statistical.test))
  if (sign %in% c("positive", "negative")) 
    message(sprintf("Subsetting weights with %s sign", sign))
  if (statistical.test == "permutation") {
    message(sprintf("Cores: %d", cores))
    message(sprintf("Number of permutations: %d", nperm))
  }
  message("\n")
  if (nperm < 100) 
    warning("A large number of permutations (at least 1000) is required for the permutation approach!\n")
  if (statistical.test == "permutation") {
    doParallel::registerDoParallel(cores = cores)
    `%dopar%` <- foreach::`%dopar%`
    null_dist_tmp <- lapply(seq_len(nperm), function(i) {
      print(sprintf("Running permutation %d/%d...", i, 
                    nperm))
      perm <- sample(ncol(data))
      W_null <- W[perm, ]
      rownames(W_null) <- rownames(W)
      colnames(W_null) <- colnames(W)
      data_null <- data[, perm]
      rownames(data_null) <- rownames(data)
      s.background <- .pcgse(data = data_null, prcomp.output = list(rotation = W_null, 
                                                                    x = Z), pc.indexes = seq_along(comps), feature.sets = feature.sets, 
                             set.statistic = set.statistic, set.test = "parametric")$statistic
      return(abs(s.background))
    })
    null_dist <- do.call("rbind", null_dist_tmp)
    colnames(null_dist) <- comps
    results <- .pcgse(data = data, prcomp.output = list(rotation = W, 
                                                        x = Z), pc.indexes = seq_along(comps), feature.sets = feature.sets, 
                      set.statistic = set.statistic, set.test = "parametric")
    s.foreground <- results$statistic
    xx <- array(unlist(null_dist_tmp), dim = c(nrow(null_dist_tmp[[1]]), 
                                               ncol(null_dist_tmp[[1]]), length(null_dist_tmp)))
    ll <- lapply(seq_len(nperm), function(i) xx[, , i] > 
                   abs(s.foreground))
    results$p.values <- Reduce("+", ll)/nperm
  }
  else {
    results <- .pcgse(data = data, prcomp.output = list(rotation = W, 
                                                        x = Z), pc.indexes = seq_along(comps), feature.sets = feature.sets, 
                      set.statistic = set.statistic, set.test = statistical.test)
  }
  pathways <- rownames(feature.sets)
  colnames(results$p.values) <- colnames(results$statistics) <- colnames(results$feature.statistics) <- comps
  rownames(results$p.values) <- rownames(results$statistics) <- pathways
  rownames(results$feature.statistics) <- colnames(data)
  if (!p.adj.method %in% p.adjust.methods) 
    stop("p.adj.method needs to be an element of p.adjust.methods")
  adj.p.values <- apply(results$p.values, 2, function(lfw) p.adjust(lfw, 
                                                                    method = p.adj.method))
  sigPathways <- lapply(comps, function(j) rownames(adj.p.values)[adj.p.values[, 
                                                                                 j] <= alpha])
  output <- list(feature.sets = feature.sets, pval = results$p.values, 
                 pval.adj = adj.p.values, feature.statistics = results$feature.statistics, 
                 set.statistics = results$statistics, sigPathways = sigPathways)
  return(output)
}
head(rownames(MSigDB_v6.0_C2_mouse), n=3)

