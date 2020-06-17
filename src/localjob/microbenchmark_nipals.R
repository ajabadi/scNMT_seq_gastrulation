library(MultiAssayExperiment)
library(here)
library(ggplot2)
library(knitr)
library(microbenchmark)
library(nipals)

mb_results <- list()
ns <- c(815)
ps <- c(100, 200, 400, 800)
mb_results <- list()

for (n in ns) {
  for (p in ps) {
    mb_results[[sprintf("N%sP%s",n,p)]] <- microbenchmark_nipals(value_mat = met_promoter, wt_mat = wt_met_promoter, subset = c(n, p), times = 3L)
    
  }
}
saveRDS(mb_results, file = here('notebook/savedata/mb_results.rds'))