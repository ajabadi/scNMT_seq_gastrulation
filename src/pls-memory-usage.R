.libPaths('~/R_libs')
library(mixOmics)
library(profmem)

pls_mem_monitor <- function(p = 50, q = 50, n = 20) {
  
  
  set.seed(100)
  create_nb_mat <- function(ncells, mu, size) {
    matrix(rnbinom(n = ncells*length(mu), mu = mu, size = size), ncol=length(mu), byrow = TRUE)
  }
  set.seed(300)
  y <- create_nb_mat(ncells = n, mu = 2^runif(q*n, -5, 5), size = 100)
  set.seed(400)
  x <- create_nb_mat(ncells = n, mu = 2^runif(p*n, -5, 5), size = 100)
  
  # Need dimnames, annoyingly enough.
  rt <- system.time({
    profmem_res <- profmem(pls_res <- tryCatch({mixOmics::pls(X = x, Y = y)}, error = function(e) e))
  })["elapsed"]
  
  if (is(pls_res, "mixo_pls")) {
    pls_res <- "successfull"
  }
  
  mem <- sum(profmem_res$bytes, na.rm = TRUE)/(1024^3)
  
  return(list(pls = pls_res, mem = mem, runtime = rt))
}

pls_loop <- function(pqs = seq(50,100,50), n=50) {
  eg <- expand.grid(p=pqs, q=pqs)
  dup_rows <- function(eg) {
    is_dup <- logical(length = dim(eg)[1])
    is_dup <- c(FALSE, sapply(2:dim(eg)[1], function(i){
      row_i <- eg[i,]
      any(apply(eg[seq_len(i-1),], 1, function(x) sum(rev(x) == row_i) == 2L))
      
    }))
    is_dup
  }
  eg <- eg[!dup_rows(eg),]
  eg$memory <- eg$pls <- eg$runtime <- NA
  
  for (i in seq_len(dim(eg)[1])) {
    pq <- eg[i,1:2, drop=TRUE]
    res <- pls_mem_monitor(p = pq$p, q = pq$q, n = n)
    eg[i, "memory"] <- res$mem
    eg[i, "pls"] <- res$pls
    eg[i, "runtime"] <- res$runtime
  }
  eg
  
}

pls_loop_cells <- function(ncells = c(20, 40), pqs = c(40,60)) {
  names(ncells) <- paste0("ncells_", ncells)
  res <- lapply(ncells, function(x){
    pls_loop(pqs = pqs, n = x)
  })
  saveRDS(res, file = sprintf("ncells-%s.rds", ncells))
  res
}

bar <- pls_loop_cells(ncells = c(200, 1000, 5000), pqs = c(5000, 10000, 20000, 40000))
# bar <- pls_loop_cells(ncells = c(20, 30), pqs = c(80, 100))
saveRDS(bar, file = here("output/pls-memory-usage.rds"))
