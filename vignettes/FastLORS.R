## ------------------------------------------------------------------------
n <- 20
p <- 50
q <- 30
k <- 4
set.seed(123)
X <- matrix(rbinom(n*p,1,0.5),n,p)
L <- matrix(rnorm(n*k),n,k) %*% t(matrix(rnorm(q*k),q,k))
B <- matrix(0, ncol(X), ncol(L))
activeSNPs <- sort(sample(c(1:nrow(B)), 20))
for(i in 1:length(activeSNPs)){
genes_influenced <- sort(sample(c(1:ncol(B)),5))
B[activeSNPs[i], genes_influenced] <- 2
}
E <- matrix(rnorm(n*q),n,q)
Y <- X %*% B + L + E

## ------------------------------------------------------------------------
library(FastLORS)
Fast_LORS_Obj <- Run_LORS(Y, X, screening = "LORS-Screening")

## ------------------------------------------------------------------------
library(FastLORS)
LORS_Obj <- Run_LORS(Y, X, method = "LORS", screening = "LORS-Screening", tune_method = "LORS", cross_valid = FALSE)

