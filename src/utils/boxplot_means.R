boxplot_means <- function(mae, obj, omic, obj_block=NULL, x_lab, y_lab, by_pheno = 'stage_lineage', background = FALSE, sign = c('all', 'positive', 'negative'), colgroup=NULL) {
  if(is.null(obj_block)) {
    obj_block <- omic
  }
  X <- obj$X[[obj_block]]
  library(data.table)
  coldata_dt <- data.table(coldata, keep.rownames = "cell", check.names = FALSE, stringsAsFactors = FALSE)
  feats <- selectVar(obj)[[obj_block]]$value
  signs <- sign(feats$value.var)
  message(sum(signs == 1), ' out of ', length(signs), ' have positive weights\n')
  sign <- match.arg(sign)
  if (sign == 'positive') {
    feats <- feats[feats$value.var > 0,,drop=FALSE]
  } else if (sign == 'negative') {
    feats <- feats[feats$value.var < 0,,drop=FALSE]
    feats$value.var <- -feats$value.var
  }
  selected_features <- rownames(feats)
  dt <- data.table(X[,selected_features], keep.rownames = "cell", check.names = FALSE, stringsAsFactors = FALSE)
  
  suff <- paste0('_', block)
  selected_features <- stringr::str_remove(selected_features, pattern = suff)
  
  if (length(selected_features) > 100) {
    message('Considering top 100 only')
    selected_features <- selected_features[1:100]
  }
  dt <- melt.data.table(dt, variable.name = "feature", value.name = "rate")
  dt <- data.table::merge.data.table(dt, coldata_dt)
  df <- data.frame(dt[,.(mean_rate = mean(rate, na.rm=TRUE)), by = c(by_pheno, "feature")])
  
  dt_bg <- data.table(X, keep.rownames = "cell", check.names = FALSE, stringsAsFactors = FALSE)
  dt_bg <- melt.data.table(dt_bg, variable.name = "feature", value.name = "rate")
  dt_bg <- data.table::merge.data.table(dt_bg, coldata_dt)
  df_bg <- data.frame(dt_bg[,.(mean_rate = mean(rate, na.rm=TRUE)), by = c(by_pheno, "feature")])
  
  df$stage <- paste0("E", df$day)
  p1 <- ggplot(df, aes_string(x = by_pheno, y = 'mean_rate'))  +  geom_violin(aes_string(fill = by_pheno))  +  geom_boxplot(aes_string(fill = by_pheno), width = 0.1, alpha = 0.4 ) + scale_fill_manual(values = colgroup) 
  # geom_jitter(width = 0.4, height = 0) +
  
  
  if(background) {
    p1 <- p1 +  geom_boxplot(data = df_bg, fill = 'grey60', outlier.shape = NA, coef = 0, alpha = 0.75)
  }
  p1 <- p1 + labs(x = x_lab, y = y_lab)  +  theme_bw()
  p1 <- p1 + theme(axis.text.x = element_text(angle = 45, vjust = 0.92, hjust = 0.95)) 
  feats <- toupper(symbols[selected_features,])
  return(list(p = p1, features = feats))
}
## ------------------------------------------------------------------------ ##
get_plotindiv_cols <- function(plotindiv) {
  pl <- plotindiv
  cols <- pl$df[,c('group', 'col')]
  cols <- cols[!duplicated(cols),]
  colgroup <- cols$col
  names(colgroup) <- cols$group
  colgroup
}