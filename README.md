# FastLORS
An R package for joint modeling to perform eQTL mapping

This R package allows users to perform joint modeling of SNPs and the expression of genes. The package provides an implementation of LORS from Yang et al. (2013) as well as FastLORS from Rhyne et al. (2018). FastLORS uses the proximal gradient method to solve the LORS problem and improves upon the computational time required for joint modeling.

The package also includes an implmentation of LORS-Screening from Yang et al. (2013).  If the number of SNPs is large, the screening option in Run_LORS should be set to "LORS-Screening" so that joint modeling will be feasible.  If the data is small enough so that modeling is feasible without screening, the screening option can simply be set to "None."

# Installation
Install FastLORS from GitHub using

```r{echo = FALSE, message = FALSE}
library(devtools)

install_github(repo = "jdrhyne2/FastLORS")
```

# Example Use

The following examples demonstrate how to use the package given a matrix of gene expression (Y) and a matrix of SNPs (X).  Note that the data used in the FastLORS paper is included in the "Data" folder on the FastLORS GitHub page.  To learn to use the package using a smaller dataset, smaller data can be simulated using the following code.

```r{echo = FALSE, message = FALSE}
n <- 20
p <- 50
q <- 30
k <- 4
set.seed(123)
X <- matrix(rbinom(n*p,1,0.5),n,p)
L <- matrix(rnorm(n*k),n,k) %*% t(matrix(rnorm(n*k),q,k))
B <- matrix(0, ncol(X), ncol(L))
activeSNPs <- sort(sample(c(1:nrow(B)), 20))
for(i in 1:length(activeSNPs)){
  genes_influenced <- sort(sample(c(1:ncol(B)),5))
  B[activeSNPs[i], genes_influenced] <- 2
}
E <- matrix(rnorm(n*q),n,q)
Y <- X %*% B + L + E
```

The following example demonstrates how to run FastLORS with LORS-Screening.  FastLORS and two-fold cross-validation will be used in tuning the parameters.

```r{echo = FALSE, message = FALSE}
FL <- Run_LORS(Y, X, method = "FastLORS", screening = "LORS-Screening", tune_method = "FastLORS")
```
The next example demonstrates how to run the original LORS method with LORS-Screening and the original LORS parameter tuning method without cross-validation.

```r{echo = FALSE, message = FALSE}
L0 <- Run_LORS(Y, X, method = "LORS", screening = "LORS-Screening", tune_method = "LORS", cross_valid = FALSE)
```

# Example Use With a Cluster

If the user's data is too large to use the FastLORS package on a local machine, the package contains functions that allow the analysis to be performed in parallel on a cluster.  This section describes an example use for this.

## Screening

If a user has access to a cluster, or more than one machine, the LORS-Screening procedure can be run in parallel using the following code.

```r{echo = FALSE, message = FALSE}
LORS_Screen_Obj <- LORS_Screen_Parallel(Y, X, chunk = 1)
```

The chunk argument determines which columns of the X matrix should be used in LORS-Screening.  The LORS_Screen_Parallel function breaks down the X matrix into batches of 1000 columns.  Thus, chunk = 1 corresponds to columns 1,2,...,1000 and chunk = 2 corresponds to columns 1001, ..., 2000 etc. Since the X matrix in the data folder has 22179 columns, this function must be run for chunk = 1,2,...,23 for this data.

The user should save LORS_Screen_Obj and then combine the estimated coefficent matrices returned by LORS_Screen_Parallel.  Following this, the SNPs selected from LORS-Screening can be found using the following code, where Bhat is the combined coefficient matrix from each run of LORS_Screen_Parallel.

```r{echo = FALSE, message = FALSE}
index_list = c()
for (i in 1:ncol(Bhat)){
  cands = abs(Bhat[,i])
  sorted_cands <- sort(cands, index.return=TRUE, decreasing=TRUE)$ix[1:nrow(X)]
  index_list <- c(index_list, sorted_cands)
}
selectedSNPs <- sort(unique(index_list))
```

## Parameter Tuning

The parameter tuning procedures of both LORS and FastLORS can also be run either on a cluster or on multiple machines.  The main function needed to perform this task is ParamTuneParallel, which outputs the selected value of lambda and 20 possible values for rho.  FastLORS or LORS can then be used to fit models on the training set using the selected lambda and one candidate value of rho. 

To use FastLORS to perform parameter tuning, the following code can be used.

```r{echo = FALSE, message = FALSE}
rho_index <- 1
params <- ParamTuneParallel(Y, X, fold = 1)
lambda <- params[[1]][[1]]
rho <- params[[1]][[2]][rho_index]
Training <- params[[2]]
Validation <- params[[3]]
myL2 <- Fast_LORS_Tuning(Y, X, rho, lambda, Training, Validation)
```

To use LORS for parameter tuning, the following code can be used

```r{echo = FALSE, message = FALSE}
rho_index <- 1
params <- ParamTuneParallel(Y, X, fold = 1)
lambda <- params[[1]][[1]]
rho <- params[[1]][[2]][5]
Training <- params[[2]]
Validation <- params[[3]]
B <- matrix(0, nrow = ncol(X), ncol = ncol(Y))
myL2 <- LORS2(Y = Y, X = X, L = NULL, Omega1 = Training, Omega2 = Validation, B = B, rho = rho, lambda = lambda)
```

Note that the rho_index in the code above can be set from 1,2,...,20.  If one has access to a cluster, 20 different files or one for each rho_index, could be submitted simulateneously, significantly reducing the time required for the parameter tuning procedure.

The user should save the myL2 object so that the optimal value of rho can be chosen based on the residual error returned by LORS or FastLORS.  The value of rho which minimizes the error returned by Fast_LORS_Tuning or LORS2 should be selected to use in the joint modeling stage.

## Joint Modeling After Parallel Parameter Tuning

Suppose the parameter values chosen by the parameter tuning procedure above are lambda_best and rho_best.  Joint modeling can then be performed using FastLORS with the following code.

```r{echo = FALSE, message = FALSE}
FL_Fit <- Fast_LORS(Y, X, rho_best, lambda_best)
```

If LORS is used for joint modeling, the following code could be used.

```r{echo = FALSE, message = FALSE}
LORS_Fit <- LORS0(Y, X, rho_best, lambda_best)
```

# Data

SNP and gene expression data of chromosome 1 for the Asian participants of the third phase of the International HapMap Project (HapMap3) are found in the "Data" folder.

The SNP dataset (X) contains 22179 SNPs for 160 participants of HapMap3 of Asian descent (CHB and JPT).  Note that LD pruning was performed using PLINK to remove SNPs highly correlated with others and SNPs with missing values were also removed.  The original data can be found at ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/hapmap3_r3/plink_format/.

The gene expression dataset (Y) contains the expression of 2010 gene probes for the 160 participants.  The original expression data can be found at https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-264/.  

# Reference

Rhyne, J., Jeng, X.J., Chi, E., and Tzeng, J.Y. (2019) FastLORS: Joint Modeling for eQTL Mapping in R.

Rhyne, J., Tzeng, J.Y., Zhang, T., and Jeng, X.J. (2018) eQTL Mapping via Effective SNP Ranking and Screening.

Yang, C., Wang, L. Zhang, S., and Zhao, H. (2013) Accounting for non-genetic factors by low-rank representation and sparse regression for eQTL mapping. Bioinformatics 29(8) 1026-1034.
