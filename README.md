# FastLORS
An R package for joint modeling to perform eQTL mapping

This R package allows users to perform joint modeling of SNPs and the expression of genes. The package provides an implementation of LORS from Yang et al. (2013) as well as Fast-LORS from Rhyne et al. (2018). Fast-LORS uses the proximal gradient method to solve the LORS problem and improves upon the computational time required for joint modeling.

The package also includes two screening methods for reducing the number of SNPs included in the modeling. LORS-Screening from Yang et al. (2013) and HC-Screening from Rhyne et al. (2018) are the two screening methods included in the package. If the data is small enough so that modeling is feasible without screening, the screening option can simply be set to "None."

# Installation
Install FastLORS from GitHub using

```r{echo = FALSE, message = FALSE}
library(devtools)

install_github("jdrhyne2/FastLORS")
```

# Example Use

The following examples demonstrate how to use the package given a matrix of gene expression (Y) and a matrix of SNPs (X).  Note that the data used in the FastLORS paper is included in the "Data" folder on the FastLORS GitHub page.  To learn to use the package using a smaller dataset, smaller data can be simulated using the following code.

```r{echo = FALSE, message = FALSE}
n <- 100
p <- 1000
q <- 500
k <- 7
set.seed(123)
X <- matrix(rbinom(n*p,1,0.33),n,p)
X[sort(sample(which(X == 1), ceiling(length(which(X == 1))/4)))] <- 2
L <- matrix(rnorm(n*k),n,k) %*% t(matrix(rnorm(q*k),q,k))
B <- matrix(0, ncol(X), ncol(L))
activeSNPs <- sort(sample(c(1:nrow(B)), 20))
for(i in 1:length(activeSNPs)){
  genes_influenced <- sort(sample(c(1:ncol(B)),5))
  B[activeSNPs[i], genes_influenced] <- 2
}
E <- matrix(rnorm(n*q),n,q)
Y <- X %*% B + L + E
```

The following example demonstrates how to run Fast-LORS with HC-Screening.  By default, Fast-LORS and two-fold cross-validation will be used in tuning the parameters.

```r{echo = FALSE, message = FALSE}
FL <- Run_LORS(Y, X, method = "Fast-LORS", screening = "HC-Screening", tune_method = "Fast-LORS")
```
The next example demonstrates how to run the original LORS method with LORS-Screening and the original LORS parameter tuning method.

```r{echo = FALSE, message = FALSE}
L0 <- Run_LORS(Y, X, method = "LORS", screening = "LORS-Screening", tune_method = "LORS", cross_valid = FALSE)
```
# Reference

Rhyne, J., Chi, E., Tzeng, J.Y., and Jeng, X.J. (2018) FastLORS: Joint Modeling for eQTL Mapping in R.

Rhyne, J., Tzeng, J.Y., Zhang, T., and Jeng, X.J. (2018) eQTL Mapping via Effective SNP Ranking and Screening.
