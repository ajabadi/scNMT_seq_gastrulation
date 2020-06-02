WPCA <-
  function(X = matrix(rlnorm(9, 3, 3), nrow = 3),
           W = matrix(rlnorm(9, meanlog = 0, sdlog = 1), nrow = 3),
           ncomp = 3,
           niter = 700,
           nrefine = 30,
           xi = 1) {
    require(MASS)

  X <- scale(X, center = TRUE, scale=FALSE)
  nvar <- dim(X)[1]; nobs <- dim(X)[2]
  P <- matrix(0, nrow = nvar, ncol = ncomp)
  C <- matrix(0, nrow = ncomp, ncol = nobs)
  WS <- rowSums(W)
  covar <- tcrossprod(WS)^xi * tcrossprod(X*W)/tcrossprod(W)
  covar[!is.finite(covar) | is.na(covar)] <- 0
  
  for (i in 1:ncomp){
    u <- rep(nvar,x = 1) / nvar
    for (j in 1:niter){
      u <- covar%*%u
      d <- as.double(sqrt(crossprod(u)))
      u <- u/d
    }
    d <- as.double(t(u)%*%covar%*%u)
    for (j in 1:nrefine){
      u <- ginv(covar - d*diag(nvar)) %*% u
      u <- u / as.double(sqrt(crossprod(u)))
      d <- as.double(t(u)%*%covar%*%u)
    }
    covar <- covar - u%*%d%*%t(u)
    P[,i] <- u
  }
  for (i in 1:nobs){
    w <- diag(W[,1])^2
    C[,i] <- ginv(t(P)%*%w%*%P) %*% (t(P)%*%w%*%X[,1])
  }
  return(list(P=P, C=C))
}

N = 3
P = 6
X = matrix(rlnorm(N*P, 3, 3), nrow = P)
W = matrix(rlnorm(N*P, meanlog = 0, sdlog = 1), nrow = P)
W = matrix(rep(0.1,N*P), nrow = P)
## why do they differ
WPCA(X = X, W = matrix(rep(0.1,N*P), nrow = P), niter = 1000)
WPCA(X = X, W = matrix(rep(1,N*P), nrow = P), niter = 1000)
prcomp(x = t(X))
