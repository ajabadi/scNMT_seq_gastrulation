create_heatmap <- Altools::create_heatmap
#' @param ... passed to Altools::create_heatmap.
#' low = '#a4eaff',
#' high = "#1C28A3",
#' seg.color = 'grey40',
#' NA_color = 'grey70',
#' seg.shape = c('none', 'cross', 'line'),
#' title = NULL,
#' xtitle = NULL,
#' ytitle = NULL,
#' show.legend = TRUE,
#' legend.title = 'Expression'
#' @examples 
#' create_assay_heatmap(gastru.mae, 'met_promoter', subset = c(20, 100))
create_assay_heatmap <- function(mae, assay, subset = NULL, ...) {
  mat <- t(assay(mae, assay))
  if (!is.null(subset)) {
    rowz <- min(subset[1], dim(mat)[1])
    colz <- min(subset[2], dim(mat)[2])
    mat <- mat[seq_len(rowz), seq_len(colz)]
  }
  create_heatmap(
    mat,
    ...
  )
  
}
## ------------------------------------------------------------------------ ##
#' ggplot an array (reduced dimensions)
#'
#' @param array Any array coercible to data.frame with columns as reduced dims
#' @param col coloring factor
#' @param dims dimensions to plot
#'
#' @return
#' @export
#'
#' @examples
#'  ggplot_redDim(matrix(rnorm(100), ncol = 2), col = sample(1:3, size = 50, replace = TRUE), x='FOO', y = 'BAR', title = 'EGG')

ggplot_redDim <- function(array, dims = c(1,2), col=NULL, shape = NULL, axis = 'UMAP_', size = 1.5, shape.legend = 'Lineage', col.legend = 'Stage', grad_cols =  NULL, ...) {
  library(ggplot2)
  df <- data.frame(array)
  if (!is.null(axis)) {
    colnames(df) <- paste0(axis,seq_along(df))
  }
  
  p <- ggplot(df)
  
  if (!is.null(col)) {
    if (!is.null(shape)) {
      p <- p + geom_point(aes_string(x = paste0(axis, dims[1]), y = paste0(axis, dims[2]), col = factor(col), shape = factor(shape)), size = size)
    } else {
      p <- p + geom_point(aes_string(x = paste0(axis, dims[1]), y = paste0(axis, dims[2]), col = factor(col)), size = size)
    }
    p <- p +  guides(col = guide_legend(col.legend))
    
  } else {
    if (!is.null(shape)) {
      p <- p + geom_point(aes_string(x = paste0(axis, dims[1]), y = paste0(axis, dims[2]), shape = factor(shape)), size = size)
    } else {
      p <- p + geom_point(aes_string(x = paste0(axis, dims[1]), y = paste0(axis, dims[2])), size = size)
    }
  }
  
  if (!is.null(shape)) {
    p <- p + scale_shape_manual(values = seq_len(length(unique(shape)))) + 
      guides(shape = guide_legend(shape.legend))
  }
  if ( !is.null(grad_cols)) {
    cols <- colorRampPalette(grad_cols)(4)
    names(cols) <- sprintf('E%s.5', 4:7)
    p <- p + scale_color_manual(values = cols)
  }
  
  p <- p + labs(...)
  
  return(p)
}
## ------------------------------------------------------------------------ ##
ggplot_redDim_expression <- function(array, dims = c(1,2), col=NULL,  axis = 'UMAP_', size = 1.5, col.legend = 'Expression', low.col = '#37a6f0', high.col = '#121296') {
  library(ggplot2)
  df <- data.frame(array)
  colnames(df) <- paste0(axis,seq_along(df))
  p <- ggplot(df)
  p <- p + geom_point(aes_string(x = paste0(axis, dims[1]), y = paste0(axis, dims[2]), col = col), size = size)
  p <- p + scale_color_gradient(low = low.col, high = high.col)
  p <- p +labs(color = col.legend)
  
  return(p)
}
## ------------------------------------------------------------------------ ##gg_sidebyside <- function(p1, p2) {
gg_sidebyside <- function(p1, p2) {
  library(egg)
  egg::ggarrange(plots = list(p1, ggplot() + theme_void(), p2), nrow = 1, widths = c(5,1,5))
}

gg_sidebyside_dims <- function(arr, dims, coldata, ...) {
  g1 <- ggplot_redDim(arr, col = coldata$stage, dims = dims, grad_cols = c('#37a6f0', '#121296'), ...)
  g2 <- ggplot_redDim(arr, col = coldata$lineage10x_2, dims = dims, grad_cols = NULL, col.legend = 'Lineage', ...)
  gg_sidebyside(g1, g2)
}
## ------------------------------------------------------------------------ ##
boxplot_means <-
  function(X, ## block.spls.res$X[[omic]]
           coldata, ## data.frame(colData(study_mae))
           selected_features,  ## selectVar(block.spls.res)[[omic]]$value
           by_pheno = 'stage', ## colnames(coldata)[i]
           x_lab,
           y_lab,
           sign = c('all', 'positive', 'negative'))
  {
    require(data.table)
    sign <- match.arg(sign)
    selected_features <- data.table(selected_features, keep.rownames = 'feature')
    
    if (sign != 'all') {
      if (sign == 'positive') {
        selected_features <- selected_features[value.var > 0]
      } else {
        selected_features <- selected_features[value.var < 0]
      }
    }
    
    if (dim(selected_features)[1] < 3) {
      stop("not many feature weights with the given sign")
    }
    
    
    selected_features <- selected_features$feature
    
    coldata_dt <- data.table(coldata, keep.rownames = "cell", check.names = FALSE, stringsAsFactors = FALSE)
    
    dt <- data.table(X, keep.rownames = "cell", check.names = FALSE, stringsAsFactors = FALSE)
    dt <- melt.data.table(dt, variable.name = "feature", value.name = "rate")
    dt <- data.table::merge.data.table(dt, coldata_dt)
    dt <- dt[,.(mean_rate = mean(rate, na.rm=TRUE)), by = c(by_pheno, "feature")]
    dt[, selected := FALSE]
    
    
    dt[feature %in% selected_features, selected := TRUE]
    ggplot(data = dt[selected==TRUE], aes_string(x = by_pheno, y = 'mean_rate', fill = by_pheno))  + geom_boxplot(data = dt, aes_string(x = by_pheno, y = 'mean_rate'),  fill = 'grey70', alpha = 0.6, outlier.shape = NA, coef= 0) +  geom_violin(alpha = 0.6)  +  geom_boxplot(alpha = 0.6)  + 
      # geom_jitter(width = 0.4, height = 0) +
      labs(x = x_lab, y = y_lab)  +  theme_bw() 
  }
## ------------------------------------------------------------------------ ##
boxplot_means_diablo <- function(diablo_obj, omic, comp = 1, coldata, by_pheno, x_lab, y_lab, sign) {
  X <- diablo_obj$X[[omic]]
  selected_features <- selectVar(diablo_obj, comp = comp, block = omic)[[1]]$value
  boxplot_means(X = X, coldata = coldata, selected_features = selected_features, by_pheno = by_pheno, x_lab = x_lab, y_lab = y_lab, sign = sign)
}