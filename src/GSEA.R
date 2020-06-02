source('/Users/alabadi/Projects/dev/R/my/utils/GenomicUtils628/R/GO.R')

diablo_feas_genes <- function(obj, block = 'rna', comps=1) {
    feas_genes <- colnames(obj$X[[block]])
    feas_genes <- rep(0, length(feas_genes)) %>% set_names(feas_genes)

    for (comp in comps){
      feas_genes[selectVar(obj, comp = comp)[[block]]$name] <- 1 ## is a named vector
    }

  feas_genes
}
