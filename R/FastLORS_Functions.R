#' svd_st
#'
#' \code{svd_st} is a function for performing soft-thresholded singular value decomposition of a matrix X
#'
#' @param X a matrix
#' @param lambda the value to to apply soft thresholding with
#' @return L the soft-thresholded singular value decomposition of X.
#' @export

svd_st <- function(X, lambda){
  mysvd <- svd(X)
  U <- mysvd[["u"]]
  V <- mysvd[["v"]]
  VT <- t(V)
  D <- diag(mysvd[["d"]])
  d <- mysvd[["d"]]
  index_list <- which(d >= rep(lambda, length(d)))
  if(length(index_list) > 1){
    W <- diag(c(d[index_list] - lambda))
    L <- U[,index_list] %*% W %*% VT[index_list,]
  }
  if(length(index_list) == 1){
    W <- matrix(d[index_list] - lambda, nrow = 1, ncol = 1)
    L <- U[,index_list] %*% W %*% VT[index_list,]
  }
  if(length(index_list) == 0){
    L <- matrix(0, nrow(X), ncol(X))
  }

  return(L)
}

#' prox_1
#'
#' \code{prox_1} is the soft-thresholding function.  Let the entries of b be b_j.  The function subtracts tau from b_j for b_j > tau, sets b_j to 0 where abs(b_j) < tau, and adds tau where b_j < -tau.
#'
#' @param b a matrix or vector
#' @param tau the value to to apply soft thresholding with
#' @return prox_b the soft-thresholding function applied to b with threshold tau
#' @export

prox_1 <- function(b, tau){
  prox_b <- sign(b)*(sapply(abs(b) - tau, FUN=function(x) {max(x,0)}))
  return(prox_b)
}


#' Fast_LORS
#'
#' \code{Fast_LORS} is a function for solving the LORS optimization problem in Can Yang et al. (2013) through the proximal gradient method
#'
#' @param Y gene expression matrix
#' @param X matrix of SNPs
#' @param rho parameter for enforcing sparsity of coefficient matrix
#' @param lambda parameter for enforcing low-rank structure of hidden factor matrix
#' @param maxiter maximum number of iterations
#' @param eps constant used when checking the convergence.  Ensures no division by 0.
#' @param tol tolerance level for convergence
#' @param verbose chooses whether details should be printed to console.  Default is FALSE.
#' @param omega_SOR the value of omega to use if applying successive over-relaxation with FastLORS.
#' @importFrom glmnet glmnet
#' @importFrom methods as
#' @return \item{B}{The estimated coefficients} \item{mu}{The estimated intercept} \item{L}{The estimated matrix of hidden factors} \item{f_val_vec}{The objective function values} \item{res_vec}{The relative change in objective function values} \item{iter}{The number of iterations}
#' @export
#' @examples
#'
#' ##Example
#'
#' ## Generate some data
#' n <- 50
#' p <- 200
#' q <- 100
#' k <- 10
#' set.seed(123)
#' X <- matrix(rbinom(n*p,1,0.5),n,p)
#' L <- matrix(rnorm(n*k),n,k) %*% t(matrix(rnorm(q*k),q,k))
#' B <- matrix(0, ncol(X), ncol(L))
#' activeSNPs <- sort(sample(c(1:nrow(B)), 20))
#' for(i in 1:length(activeSNPs)){
#' genes_influenced <- sort(sample(c(1:ncol(B)),5))
#' B[activeSNPs[i], genes_influenced] <- 2
#' }
#' E <- matrix(rnorm(n*q),n,q)
#' Y <- X %*% B + L + E
#'
#' rho <- runif(1,3,5)
#' lambda <- runif(1,3,5)
#'
#' ## Usage
#' Fast_LORS(Y, X, rho, lambda)

Fast_LORS <- function(Y, X, rho, lambda, maxiter = 5000, eps = 2.2204e-16, tol = 1e-4, verbose = FALSE, omega_SOR = 1.999) {

  ### Initial Setup

  n <- nrow(Y)
  q <- ncol(Y)
  p <- ncol(X)
  ones <- matrix(1, nrow = n, ncol = 1)

  B <- matrix(0, nrow = p, ncol = q)
  mu <- matrix(0, nrow = 1, ncol = q)
  L <- matrix(0, nrow = n, ncol = q)

  ### Proximal Mapping Functions

  ## Nuclear Norm:  Soft Thresholded SVD
  svd_st <- function(X, lambda){
    mysvd <- svd(X)
    U <- mysvd[["u"]]
    V <- mysvd[["v"]]
    VT <- t(V)
    D <- diag(mysvd[["d"]])
    d <- mysvd[["d"]]
    index_list <- which(d >= rep(lambda, length(d)))
    if(length(index_list) > 1){
      W <- diag(c(d[index_list] - lambda))
      L <- U[,index_list] %*% W %*% VT[index_list,]
    }
    if(length(index_list) == 1){
      W <- matrix(d[index_list] - lambda, nrow = 1, ncol = 1)
      L <- U[,index_list] %*% W %*% VT[index_list,]
    }
    if(length(index_list) == 0){
      L <- matrix(0, nrow(X), ncol(X))
    }

    return(L)
  }

  ## 1 norm:  Soft Threshold Function
  prox_1 <- function(b, tau){
    prox_b <- sign(b)*(sapply(abs(b) - tau, FUN=function(x) {max(x,0)}))
    return(prox_b)
  }

  fval_old <- 0
  f_val_vec <- c()
  res_vec <- c()

  t_L <- 1
  t_B <- 1/norm(t(X) %*% X, type = "2")
  t_mu <- 1/nrow(Y)

  ones <- matrix(1, nrow = n, ncol = 1)

  if(verbose == TRUE){
    for(iter in 1:maxiter){

      #### Update B, mu, and L
      L <- svd_st(L - t_L * (X %*% B + ones %*% mu + L - Y), t_L * lambda)
      B <- prox_1(B - t_B * t(X) %*% (X %*% B + ones %*% mu + L - Y), t_B * rho)
      mu <- mu - t_mu * t(ones) %*% (X %*% B + ones %*% mu + L - Y)

      ## Update via SOR

      if(iter > 1){
        L_update <- (1-omega_SOR) * L_update + omega_SOR * L
        B_update <- (1-omega_SOR) * B_update + omega_SOR * B
        mu_update <- (1-omega_SOR) * mu_update + omega_SOR * mu

        L <- L_update
        B <- B_update
        mu <- mu_update
      }

      if(iter == 1){
        L_update <- L
        B_update <- B
        mu_update <- mu

        L <- L_update
        B <- B_update
        mu <- mu_update
      }

      #### Check Convergence

      dum <- c(Y - X %*% B - ones%*%mu - L)
      fval <- 0.5 * norm(dum, type = "2")^2 + rho * sum(abs(B)) + lambda*sum(svd(L)[["d"]])

      res = abs(fval-fval_old)/abs(fval_old+eps)

      print(paste('Iter ', iter, 'fval', fval, 'res', res))

      fval_old <- fval
      f_val_vec <- c(f_val_vec, fval)
      res_vec <- c(res_vec, res)

      if (res < tol){
        break
      }

      if(iter > 1){
        if(f_val_vec[iter] > f_val_vec[iter-1]){
          break
        }
      }

    }

    if(iter < maxiter & res >= tol){
      print("Beginning LORS Updates")
      for (iter_LORS in (iter+1):maxiter){

        #### Compute L
        mysvd <- svd(Y - X%*%B-matrix(1, nrow = n, ncol = 1) %*% mu)
        U <- mysvd[["u"]]
        V <- mysvd[["v"]]
        VT <- t(V)
        d <- mysvd[["d"]]
        index_list <- which(d >= lambda)
        W <- diag(d[index_list] - lambda)
        L <- U[,index_list] %*% W %*% VT[index_list,]

        for (j in 1:q){
          fit <- glmnet(X, Y[,j]-L[,j], family = "gaussian", lambda = rho/n, standardize = FALSE)
          a0 <- fit[["a0"]]
          old_beta <- fit[["beta"]]
          my_beta <- as(old_beta, "matrix")
          B[,j] <- my_beta
          mu[,j] <- a0
        }

        dum <- c(Y - X%*%B - matrix(1, nrow = n, ncol = 1)%*%mu - L)
        fval <- 0.5 * t(dum) %*% dum + rho * sum(abs(B)) + lambda*sum(abs(diag(W)))

        res <- abs(fval-fval_old)/abs(fval_old+eps)

        print(paste('Iter ', iter_LORS, 'fval', fval, 'res', res))

        fval_old <- fval
        f_val_vec <- c(f_val_vec, fval)
        res_vec <- c(res_vec, res)

        if (res < tol){
          break
        }

      }
    }
  }

  if(verbose == FALSE){
    for(iter in 1:maxiter){

      #### Update B, mu, and L
      L <- svd_st(L - t_L * (X %*% B + ones %*% mu + L - Y), t_L * lambda)
      B <- prox_1(B - t_B * t(X) %*% (X %*% B + ones %*% mu + L - Y), t_B * rho)
      mu <- mu - t_mu * t(ones) %*% (X %*% B + ones %*% mu + L - Y)

      ## Update via SOR

      if(iter > 1){
        L_update <- (1-omega_SOR) * L_update + omega_SOR * L
        B_update <- (1-omega_SOR) * B_update + omega_SOR * B
        mu_update <- (1-omega_SOR) * mu_update + omega_SOR * mu

        L <- L_update
        B <- B_update
        mu <- mu_update
      }

      if(iter == 1){
        L_update <- L
        B_update <- B
        mu_update <- mu

        L <- L_update
        B <- B_update
        mu <- mu_update
      }

      #### Check Convergence

      dum <- c(Y - X %*% B - ones%*%mu - L)
      fval <- 0.5 * norm(dum, type = "2")^2 + rho * sum(abs(B)) + lambda*sum(svd(L)[["d"]])

      res = abs(fval-fval_old)/abs(fval_old+eps)

      fval_old <- fval
      f_val_vec <- c(f_val_vec, fval)
      res_vec <- c(res_vec, res)

      if (res < tol){
        break
      }

      if(iter > 1){
        if(f_val_vec[iter] > f_val_vec[iter-1]){
          break
        }
      }

    }

    if(iter < maxiter & res >= tol){
      for (iter_LORS in (iter+1):maxiter){

        #### Compute L
        mysvd <- svd(Y - X%*%B-matrix(1, nrow = n, ncol = 1) %*% mu)
        U <- mysvd[["u"]]
        V <- mysvd[["v"]]
        VT <- t(V)
        d <- mysvd[["d"]]
        index_list <- which(d >= lambda)
        W <- diag(d[index_list] - lambda)
        L <- U[,index_list] %*% W %*% VT[index_list,]

        for (j in 1:q){
          fit <- glmnet(X, Y[,j]-L[,j], family = "gaussian", lambda = rho/n, standardize = FALSE)
          a0 <- fit[["a0"]]
          old_beta <- fit[["beta"]]
          my_beta <- as(old_beta, "matrix")
          B[,j] <- my_beta
          mu[,j] <- a0
        }

        dum <- c(Y - X%*%B - matrix(1, nrow = n, ncol = 1)%*%mu - L)
        fval <- 0.5 * t(dum) %*% dum + rho * sum(abs(B)) + lambda*sum(abs(diag(W)))

        res <- abs(fval-fval_old)/abs(fval_old+eps)

        fval_old <- fval
        f_val_vec <- c(f_val_vec, fval)
        res_vec <- c(res_vec, res)

        if (res < tol){
          break
        }


      }
    }
  }

  ### Return the B, L, mu, objective function values, residual values, and total iterates
  return(list("B" = B, "L" = L, "mu" = mu, "f_vals" = f_val_vec, "res_vec" = res_vec, "iter" = length(f_val_vec)))
}

#' LORS0
#'
#' \code{LORS0} is a function for solving the LORS optimization problem through the method described in Can Yang et al. (2013).
#' This function is adapted from the authors MATLAB implementation
#'
#' @param Y gene expression matrix
#' @param X matrix of SNPs
#' @param rho parameter for enforcing sparsity of coefficient matrix
#' @param lambda parameter for enforcing low-rank structure of hidden factor matrix
#' @param maxiter maximum number of iterations
#' @param eps constant used when checking the convergence.  Ensures no division by 0.
#' @param tol tolerance level for convergence
#' @param verbose chooses whether details should be printed to console.  Default is FALSE.
#' @importFrom methods as
#' @return \item{B}{The estimated coefficients} \item{mu}{The estimated intercept} \item{L}{The estimated matrix of hidden factors} \item{f_val_vec}{The objective function values} \item{res_vec}{The relative change in objective function values}
#' @export
#' @examples
#'
#' ##Example
#'
#' #' ## Generate some data
#' n <- 50
#' p <- 200
#' q <- 100
#' k <- 10
#' set.seed(123)
#' X <- matrix(rbinom(n*p,1,0.5),n,p)
#' L <- matrix(rnorm(n*k),n,k) %*% t(matrix(rnorm(q*k),q,k))
#' B <- matrix(0, ncol(X), ncol(L))
#' activeSNPs <- sort(sample(c(1:nrow(B)), 20))
#' for(i in 1:length(activeSNPs)){
#' genes_influenced <- sort(sample(c(1:ncol(B)),5))
#' B[activeSNPs[i], genes_influenced] <- 2
#' }
#' E <- matrix(rnorm(n*q),n,q)
#' Y <- X %*% B + L + E
#'
#' rho <- runif(1,3,5)
#' lambda <- runif(1,3,5)
#' LORS0(Y, X, rho, lambda)

LORS0 <- function(Y, X, rho, lambda, maxiter = 1000, eps = 2.2204e-16, tol = 1e-4, verbose = FALSE){

  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)

  B <- matrix(0, nrow = p, ncol=q)
  mu <- matrix(0, nrow=1, ncol=q)

  fval_old <- 0
  f_val_vec <- c()
  res_vec <- c()
  if(verbose == TRUE){
    for (iter in 1:maxiter){

      #### Compute L
      mysvd <- svd(Y - X%*%B-matrix(1, nrow = n, ncol = 1) %*% mu)
      U <- mysvd[["u"]]
      V <- mysvd[["v"]]
      VT <- t(V)
      d <- mysvd[["d"]]
      index_list <- which(d >= lambda)
      W <- diag(d[index_list] - lambda)
      L <- U[,index_list] %*% W %*% VT[index_list,]

      for (j in 1:q){
        fit <- glmnet(X, Y[,j]-L[,j], family = "gaussian", lambda = rho/n, standardize = FALSE) #### need to check this
        a0 <- fit[["a0"]]
        old_beta <- fit[["beta"]]
        my_beta <- as(old_beta, "matrix")
        B[,j] <- my_beta
        mu[,j] <- a0
      }

      dum <- c(Y - X%*%B - matrix(1, nrow = n, ncol = 1)%*%mu - L)
      fval <- 0.5 * t(dum) %*% dum + rho * sum(abs(B)) + lambda*sum(abs(diag(W)))

      res <- abs(fval-fval_old)/abs(fval_old+eps)

      print(paste('Iter ', iter, 'fval', fval, 'res', res))

      if (res < tol){
        break
      }

      fval_old <- fval
      f_val_vec <- c(f_val_vec, fval)
      res_vec <- c(res_vec, res)
    }
  }

  if(verbose == FALSE){
    for (iter in 1:maxiter){

      #### Compute L
      mysvd <- svd(Y - X%*%B-matrix(1, nrow = n, ncol = 1) %*% mu)
      U <- mysvd[["u"]]
      V <- mysvd[["v"]]
      VT <- t(V)
      d <- mysvd[["d"]]
      index_list <- which(d >= lambda)
      W <- diag(d[index_list] - lambda)
      L <- U[,index_list] %*% W %*% VT[index_list,]

      for (j in 1:q){
        fit <- glmnet(X, Y[,j]-L[,j], family = "gaussian", lambda = rho/n, standardize = FALSE) #### need to check this
        a0 <- fit[["a0"]]
        old_beta <- fit[["beta"]]
        my_beta <- as(old_beta, "matrix")
        B[,j] <- my_beta
        mu[,j] <- a0
      }

      dum <- c(Y - X%*%B - matrix(1, nrow = n, ncol = 1)%*%mu - L)
      fval <- 0.5 * t(dum) %*% dum + rho * sum(abs(B)) + lambda*sum(abs(diag(W)))

      res <- abs(fval-fval_old)/abs(fval_old+eps)

      if (res < tol){
        break
      }

      fval_old <- fval
      f_val_vec <- c(f_val_vec, fval)
      res_vec <- c(res_vec, res)
    }
  }
  return (list("B" = B, "L" = L, "mu" = mu, "f_val_vec" = f_val_vec, "res_vec" = res_vec))
}

#' LORS2
#'
#' \code{LORS2} is a function used in parameter tuning in LORS. See the parameter tuning section described in Can Yang et al. (2013).
#' This function is adapted from the authors MATLAB implementation
#'
#' @param Y gene expression matrix
#' @param X matrix of SNPs
#' @param L matrix of hidden factors
#' @param rho parameter for enforcing sparsity of coefficient matrix
#' @param lambda parameter for enforcing low-rank structure of hidden factor matrix
#' @param Omega1 Boolean matrix for training data
#' @param Omega2 Boolean matrix for validation data
#' @param B a matrix of coefficients for the SNPs
#' @param maxIter the maximum number of iterations
#' @param tol tolerance level for convergence
#' @importFrom glmnet glmnet
#' @importFrom methods as
#' @return \item{B}{The estimated coefficients} \item{mu}{The estimated intercept} \item{L}{The estimated matrix of hidden factors} \item{Err}{The residual sum of squares on the validation set}
#' @export
#' @examples
#'
#' ##Example
#'
#' #' ## Generate some data
#' n <- 50
#' p <- 200
#' q <- 100
#' k <- 10
#' set.seed(123)
#' X <- matrix(rbinom(n*p,1,0.5),n,p)
#' L <- matrix(rnorm(n*k),n,k) %*% t(matrix(rnorm(q*k),q,k))
#' B <- matrix(0, ncol(X), ncol(L))
#' activeSNPs <- sort(sample(c(1:nrow(B)), 20))
#' for(i in 1:length(activeSNPs)){
#' genes_influenced <- sort(sample(c(1:ncol(B)),5))
#' B[activeSNPs[i], genes_influenced] <- 2
#' }
#' E <- matrix(rnorm(n*q),n,q)
#' Y <- X %*% B + L + E
#'
#' Omega0 <- !(is.na(Y))
#' mask <- matrix(runif(nrow(Y)*ncol(Y)) > 0.5, nrow = nrow(Y), ncol = ncol(Y))
#' Omega1 <- Omega0 & mask
#' Omega2 <- Omega0 & !mask
#' rho <- runif(1,3,5)
#' lambda <- runif(1,3,5)
#' tol <- 1e-4
#'
#' ## Usage
#' LORS2(Y, X, L, Omega1, Omega2, B, rho, lambda, tol)

LORS2 <- function(Y, X, L, Omega1, Omega2, B, rho, lambda, tol, maxIter = 1000){

  eps = 2.2204e-16

  n = nrow(X)
  p = ncol(X)
  q = ncol(Y)

  mu = matrix(0, nrow=1, ncol=q)
  fval_old = 0

  maxInnerIts = 50

  #### initialization

  if (is.null(L) == 1){
    L = Y
  }

  energy = Inf

  for (iter in 1:maxIter){

    for (innerIts in 1:maxInnerIts){

      Z = Y - X %*% B - matrix(1, nrow = n, ncol = 1) %*% mu;
      C = Z*Omega1 + L*(!Omega1)

      U = svd(C)[["u"]]
      D = diag(svd(C)[["d"]])
      V = svd(C)[["v"]]
      VT = t(V)
      # soft impute
      d = diag(D)

      index_list = which(d > lambda)

      W = diag(d[index_list] - lambda, nrow = length(index_list), ncol = length(index_list))
      L = U[,index_list] %*% W %*% VT[index_list,]

      Lnorm = sum(d[index_list] - lambda)
      energy_old = energy

      cL = c(L)
      cZ = c(Z)
      newvec = matrix(cL[c(Omega1 == 1)]-cZ[c(Omega1 == 1)], ncol = 1)
      energy = lambda * Lnorm + norm(newvec,'F')/2

      mydiff = abs(energy - energy_old) / energy_old

      if (!is.na(mydiff)){
        if (abs(energy - energy_old) / energy_old < tol){
          break
        }
      }
    }

    ### Compute B using glmnet

    for (j in 1:q){
      fit <- glmnet(X[Omega1[,j],], matrix(Y[Omega1[,j],j]) - matrix(L[Omega1[,j],j]), family = "gaussian", lambda = rho/sum(Omega1[,j]), standardize = FALSE)
      a0 <- fit[["a0"]]
      old_beta <- fit[["beta"]]
      my_beta <- as(old_beta, "matrix")
      B[,j] <- my_beta
      mu[,j] <- a0
    }

    # Convergence
    residual = Y - X%*%B - matrix(1, nrow = n, ncol = 1) %*% mu - L
    dum = residual*Omega1
    dum = c(dum)
    fval = 0.5 * t(dum) %*% dum + rho * sum(abs(c(B))) + lambda*sum(abs(diag(W)))
    res = abs(fval-fval_old)/abs(fval_old+eps)

    if (res < tol){
      break
    }

    fval_old <- fval

  }

  err = residual*Omega2
  Err = t(c(err)) %*% c(err) / sum(Omega2)

  return (list("B" = B, "mu" = mu, "L" = L, "Err" = Err))
}

#' GetMaxRho
#'
#' \code{GetMaxRho} is a function used to determine the maximum value as a candidate for rho.  See the parmeter tuning section of Yang et al. (2013)
#' Note:  This function is adapted from the LORS MATLAB implementation
#'
#' @param Y gene expression matrix
#' @param X matrix of SNPs
#' @param L matrix of hidden factors
#' @param Omega0 Boolean matrix of observed entries
#' @return MaxRho The maximum value of rho to be used in parameter tuning
#' @export
#' @examples
#'
#' ##Example
#'
#' ## Generate some data
#' n <- 50
#' p <- 200
#' q <- 100
#' k <- 10
#' set.seed(123)
#' X <- matrix(rbinom(n*p,1,0.5),n,p)
#' L <- matrix(rnorm(n*k),n,k) %*% t(matrix(rnorm(q*k),q,k))
#' B <- matrix(0, ncol(X), ncol(L))
#' activeSNPs <- sort(sample(c(1:nrow(B)), 20))
#' for(i in 1:length(activeSNPs)){
#' genes_influenced <- sort(sample(c(1:ncol(B)),5))
#' B[activeSNPs[i], genes_influenced] <- 2
#' }
#' E <- matrix(rnorm(n*q),n,q)
#' Y <- X %*% B + L + E
#' Omega0 <- !(is.na(Y))
#'
#' ## Usage
#' GetMaxRho(X, Y, L, Omega0)
GetMaxRho <- function(X, Y, L, Omega0){
  q = ncol(Y)
  maxrho = matrix(0, nrow = q, ncol = 1)

  for(i in 1:q){
    maxrho[i,1] = max(abs(t(X) %*% ((matrix(Y[,i], nrow = nrow(Y), ncol = 1) - matrix(L[,i], nrow = nrow(L), ncol = 1)) * matrix(Omega0[,i], nrow = nrow(Omega0), ncol = 1))))
  }

  MaxRho = max(maxrho)
  return(MaxRho)
}

#' SVT
#'
#' \code{SVT} is a function to perform soft-thresholded singular value decomposition.  It is used to get an initial estimate for L.
#' Note:  This function is adapted from the LORS MATLAB implementation
#'
#' @param Y gene expression matrix
#' @param lambda a tuning parameter
#' @return L the singular value decomposition of Y, soft-thresholded with lambda.
#' @export
#' @examples
#'
#' ##Example
#' set.seed(123)
#'
#' Y <-matrix(rnorm(50*100, 7,1), nrow = 50, ncol = 100)
#' lambda <- runif(1,3,5)
#' SVT(Y, lambda)
SVT <- function(Y, lambda)
{
  mysvd = svd(Y)
  U <- mysvd[["u"]]
  V <- mysvd[["v"]]
  VT <-t(V)
  d <- mysvd[["d"]]
  index_list = which(d >= lambda)
  W = diag(d[index_list] - lambda)
  L <- U[,index_list] %*% W %*% VT[index_list,]
  return(L)
}

#' softImpute
#'
#' \code{softImpute} is a function from Mazudmer et al. (2010). It solves the problem min || X - Z ||_Omega + alpha || Z ||_Nulear and is used in parameter tuning for LORS.
#' Note:  This function is adapted from the LORS MATLAB implementation
#'
#' @param X a (possibly) incomplete matrix
#' @param Z the target matrix
#' @param Omega0 Boolean matrix of observed entries
#' @param Omega1 Boolean matrix of training entries
#' @param Omega2 Boolean matrix of validation entries
#' @param alpha0 initial tuning parameter
#' @param maxRank maximum rank of the solution
#' @return \item{Z}{Estimate of the target matrix} \item{Err}{Squared Error of the difference between X and Z on the validation set} \item{rank_alpha}{The rank of the estimates} \item{znorm}{The sum of the soft-thresholded singular values of the estimates} \item{Alpha}{The tuning parameters used}
#' @export
#'
#'

softImpute <- function(X,Z,Omega0,Omega1,Omega2,alpha0,maxRank){

  X_0 = X*Omega0
  if (is.null(Z) == TRUE) {
    Z = X_0
  }

  if (is.null(alpha0) == TRUE){
    my_svd = svd(X_0)
    UU = my_svd[["u"]]
    DD = diag(my_svd[["d"]])
    VV = my_svd[["v"]]
    alpha0 = DD[2,2]
  }

  if (is.null(maxRank) == TRUE){
    maxRank = -1
  }

  # parameters
  eta = 0.9
  epsilon = 1e-4
  maxInnerIts = 50

  ## soft-impute

  # 1. initialize
  alpha = alpha0

  Err = c()
  rank_alpha = c()
  znorm = c()
  Alpha = c()

  while (TRUE){
    energy = Inf
    for (innerIts in 1:maxInnerIts){
      # (a)i
      C = X*Omega1 + Z*(!Omega1)
      my_svd2 = svd(C)
      U = my_svd2[["u"]]
      D = diag(my_svd2[["d"]])
      V = my_svd2[["v"]]
      VT = t(V)
      # soft impute
      d = my_svd2[["d"]]
      idx = which(d > alpha)
      if (length(idx) > 0){
        Z = matrix(U[,idx], ncol = length(idx)) %*% diag( d[idx] - alpha , nrow = length(idx)) %*% matrix(VT[idx,], nrow = length(idx))
      }

      else{
        Z = matrix(0, nrow = nrow(U), ncol = ncol(VT))
      }

      Znorm = sum(d[idx]-alpha)
      energy_old = energy
      cZ = c(Z)
      cX = c(X)
      newvec = matrix(cZ[c(Omega1 == 1)]-cX[c(Omega1 == 1)], ncol = 1)
      energy = alpha*Znorm + norm(newvec, 'F')/2

      mydiff = abs(energy - energy_old) / energy_old
      if (!is.na(mydiff)){
        if (abs(energy - energy_old) / energy_old < epsilon){
          break
        }
      }
    }

    e = X * Omega2 - Z * Omega2
    err2 = t(c(e)) %*%  c(e)
    Err = cbind(Err, err2)
    znorm = cbind(znorm, Znorm)
    k = length(idx)
    rank_alpha = cbind(rank_alpha, k)
    Alpha = cbind(Alpha, alpha)

    if (k <= maxRank && alpha > 1e-3){
      alpha = alpha * eta
    }

    else{
      break
    }

  }
  return(list("Z" = Z, "Err" = Err, "rank_alpha" = rank_alpha, "znorm" = znorm, "Alpha" = Alpha))
}

#' linspace
#'
#' \code{linspace} is a function to space a sequence linearly from x1 to x2
#'
#' @param x1 a starting point
#' @param x2 an ending point
#' @param n length of sequence
#' @export
#' @return a linear spaced sequence from x1 to x2 of length n
#' @examples
#' linspace(100,10,5)
#'
linspace <- function (x1, x2, n = 100) {
  stopifnot(is.numeric(x1), is.numeric(x2), length(x1) == 1, length(x2) == 1)
  n <- floor(n)
  if (n <= 1)
    x2
  else seq(x1, x2, length.out = n)
}

#' logspace
#'
#' \code{logspace} is a function to space a sequence evenly on the log scale from x1 to x2
#'
#' @param x1 a starting point
#' @param x2 an ending point
#' @param n length of sequence
#' @return a sequence from x1 to x2 of length n spaced evenly on the log scale
#' @export
#' @examples
#' ##Example
#' logspace(100,10,5)
#'
logspace <- function (x1, x2, n = 50) {
  if (x2 == pi)
    x2 <- log10(x2)
  10^linspace(x1, x2, n)
}

#' InitialEst
#' \code{InitialEst} is a function to build an initial estimate for B
#'
#' @param Y gene expression matrix
#' @param X matrix of SNPs
#' @param lambda tuning parameter
#' @importFrom stats runif
#' @return B An initial estimate of the coefficient matrix
#' @export
#' @examples
#'
#' ## Example
#' ## Generate some data
#' n <- 50
#' p <- 200
#' q <- 100
#' k <- 10
#' set.seed(123)
#' X <- matrix(rbinom(n*p,1,0.5),n,p)
#' L <- matrix(rnorm(n*k),n,k) %*% t(matrix(rnorm(q*k),q,k))
#' B <- matrix(0, ncol(X), ncol(L))
#' activeSNPs <- sort(sample(c(1:nrow(B)), 20))
#' for(i in 1:length(activeSNPs)){
#'   genes_influenced <- sort(sample(c(1:ncol(B)),5))
#'   B[activeSNPs[i], genes_influenced] <- 2
#' }
#' E <- matrix(rnorm(n*q),n,q)
#' Y <- X %*% B + L + E
#'
#' Init_est <- InitialEst(Y,X)

InitialEst <- function(Y, X, lambda = NULL){

  n = nrow(Y)
  q = ncol(Y)
  p = ncol(X)

  Omega0 = !(is.na(Y))
  Y[is.na(Y)] <- 0

  if(is.null(lambda) == TRUE){
    # first use soft-impute to select a reasonable lambda

    mask = matrix(runif(n*q) > 0.5, nrow = n, ncol = q)
    Omega1 = Omega0 & mask
    Omega2 = Omega0 & !mask

    maxRank = min(n,q)/2
    mySI= softImpute(Y,NULL,Omega0, Omega1, Omega2, NULL,maxRank)
    Z = mySI[["Z"]]
    Err = mySI[["Err"]]
    rank_alpha = mySI[["rank_alpha"]]
    Znorm = mySI[["znorm"]]
    Alpha = mySI[["Alpha"]]
    bestind_lam = min(Err)==Err;

    lambda = min(Alpha[bestind_lam])
  }

  myB <- c()
  print("Building Initial Estimate")
  for(SNP_col in 1:ncol(X)){
    #print(paste("Building initial estimate: On column", SNP_col,"of",ncol(X)))
    X1 <- matrix(X[,SNP_col], ncol = 1)
    LS <- LORSscreen(Y, X1, lambda, 0.01)
    B_row <- LS$B
    myB <- rbind(myB, B_row)
  }

  return(list("B"  = myB))
}

#' standardizeBhat
#' \code{standardizeBhat} is a function used to standardize a coefficient matrix
#'
#' @param Y gene expression matrix
#' @param X matrix of SNPs
#' @param Bhat a coefficient matrix
#' @importFrom stats var
#' @return A standardized estimate of the coefficient matrix.
#' @export
#'
#'

standardizeBhat <- function(Y, X, Bhat){
  Bhat_standard <- Bhat
  n <- nrow(Y)
  print("Standardizing Initial Estimate")
  for(i in 1:ncol(X)){
    for(j in 1:ncol(Y)){
      var_B_ij <- (t(X[,i]) %*% X[,i])^(-1) * var(Y[,j] - X[,i] * Bhat[i,j])
      Bhat_standard[i,j] <- Bhat[i,j]/sqrt(var_B_ij)
    }
  }
  return(Bhat_standard)
}

#' rankHC
#' \code{rankHC} is a function used to rank a Bhat matrix by higher criticism statistics
#'
#' @param Bhat_standardized a standardized coefficient matrix
#' @return \item{index}{The indices of the sorted Higher Criticism values} \item{HC_vec}{The higher criticism statistics for the standardized Bhat matrix}
#' @export
#'
#'

rankHC <- function(Bhat_standardized){
  HC_vec <- c()
  q <- ncol(Bhat_standardized)
  for (j in 1:nrow(Bhat_standardized)){
    if(length(which(Bhat_standardized[j,] == 0)) != ncol(Bhat_standardized)){  #### for very sparse matrices, some rows may be all 0.  Don't want to calculate HC for these
      t_vec <- sort(abs(Bhat_standardized[j,]), decreasing=FALSE)
      t_vec <- t_vec[1:(length(t_vec) - 2)] ### values at the end of t_vec can be extreme

      HC = max(sqrt(q)*(S(t_vec, Bhat_standardized[j,])/q - 2 * survival(t_vec)) / sqrt(2 * survival(t_vec) * (1 - 2 * survival(t_vec))))
      HC_vec <- c(HC_vec, HC)
    }

    if(length(which(Bhat_standardized[j,] == 0)) == ncol(Bhat_standardized)){  #### for very sparse matrices, some rows may be all 0.  HC = 0 for these
      HC_vec <- c(HC_vec, 0)
    }

  }

  sorted_HC <- sort(HC_vec, index.return=TRUE, decreasing=TRUE)

  return(list("index" = sorted_HC$ix, "HC_vec" = HC_vec))
}



#' S
#' \code{S} is a function used internally in rankHC.  It calculates empirical cdf's.
#'
#' @param t cutoff value of empirical cdf
#' @param my_matrix a coefficient matrix
#' @importFrom stats ecdf
#' @return The empirical distribution function of the coefficients evaluated at t
#' @export
#'
#'
S <- function(t,my_matrix){
  length(my_matrix)*(1 - ecdf(abs(my_matrix))(t))
}

#' survival
#' \code{survival} is a function used internally in rankHC.  It calculates the survival function of the standard normal distribution
#'
#' @param t a cutoff value
#' @importFrom stats pnorm
#' @return s 1 - Pr(Z <= t) where Z is a standard normal random variable
#' @export
#'
#'
survival <- function(t){
  s <- 1 - pnorm(t, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
  return (s)
}


#' LORSscreen
#' \code{LORSscreen} is a function to solve the LORS-Screening optimization problem in Yang et al. (2013)
#'
#' @param Y gene expression matrix
#' @param X a SNP
#' @param lambda tuning parameter
#' @param tol a tolerance level
#' @return B the estimated coefficients for the SNP
#' @return L the estimated hidden factors
#' @return mu the estimate for the intercept
#' @export
#' @examples
#'
#' ##Example
#'
#' ## Generate some data
#' n <- 50
#' p <- 200
#' q <- 100
#' k <- 10
#' set.seed(123)
#' X <- matrix(rbinom(n*p,1,0.5),n,p)
#' L <- matrix(rnorm(n*k),n,k) %*% t(matrix(rnorm(q*k),q,k))
#' B <- matrix(0, ncol(X), ncol(L))
#' activeSNPs <- sort(sample(c(1:nrow(B)), 20))
#' for(i in 1:length(activeSNPs)){
#' genes_influenced <- sort(sample(c(1:ncol(B)),5))
#' B[activeSNPs[i], genes_influenced] <- 2
#' }
#' E <- matrix(rnorm(n*q),n,q)
#' Y <- X %*% B + L + E
#'
#' ## Usage to build initial estimate
#'
#' Bhat_initial <- c()
#' for(SNP_col in 1:ncol(X)){
#'    X1 <- matrix(X[,SNP_col], ncol = 1)
#'    LS <- LORSscreen(Y, X1, lambda = 0.1, 0.01)
#'    B_row <- LS$B
#'    Bhat_initial <- rbind(Bhat_initial, B_row)
#' }

LORSscreen <- function(Y, X, lambda, tol){

  eps = 2.2204e-16

  n = nrow(Y)
  q = ncol(Y)

  B = matrix(0, nrow = 1, ncol=q)
  mu = matrix(0, nrow=1, ncol=q)

  maxIter = 100

  fval_old = 0
  for (iter in 1:maxIter){

    #### Compute L
    mysvd <-svd(Y - X%*%B-matrix(1, nrow = n, ncol = 1) %*% mu)
    U <- mysvd[["u"]]
    V <- mysvd[["v"]]
    VT <- t(V)
    d <- mysvd[["d"]]
    index_list = which(d >= lambda)
    nl = length(index_list)
    W = matrix(0, nrow = nl, ncol = nl)
    for (s in 1:nl){
      W[index_list[s],index_list[s]] = d[index_list[s]] - lambda
    }
    L <- U[,index_list] %*% W %*% VT[index_list,]

    #### Solve the Least Squares Problem

    Z <- Y - L
    A <- cbind(X,1)
    myqr <- qr(A)
    Q <- qr.Q(myqr)
    R <- qr.R(myqr)
    B_QR <- solve(t(R) %*% R) %*% t(A) %*% Z

    B <- matrix(B_QR[1,], nrow = 1, ncol = q)
    mu <- matrix(B_QR[2,], nrow = 1, ncol = q)

    dum = Y - X%*%B - matrix(1, nrow = n, ncol = 1)%*%mu - L
    dum = c(dum)
    fval = 0.5 * t(dum) %*% dum + lambda*sum(abs(diag(W)))

    res = abs(fval-fval_old)/abs(fval_old+eps)

    if (res < tol){
      break
    }

    fval_old <- fval
  }
  return (list("B" = B, "L" = L, "mu" = mu))
}

#' HC_Screening
#' \code{HC_Screening} is a function to apply the HC-Screening screening method of Rhyne et al. (2018) (In progress)
#' Note:  HC-Screening ranks SNPs by their higher criticism statistics and selects the top n, where n is the number of samples
#' @param Y gene expression matrix
#' @param X matrix of SNPs
#' @importFrom utils head
#' @return selectedSNPs The SNPs selected by HC-Screening
#' @export
#' @examples
#'
#' ## Example
#'
#' ## Generate some data
#' n <- 50
#' p <- 200
#' q <- 100
#' k <- 10
#' set.seed(123)
#' X <- matrix(rbinom(n*p,1,0.5),n,p)
#' L <- matrix(rnorm(n*k),n,k) %*% t(matrix(rnorm(q*k),q,k))
#' B <- matrix(0, ncol(X), ncol(L))
#' activeSNPs <- sort(sample(c(1:nrow(B)), 20))
#' for(i in 1:length(activeSNPs)){
#' genes_influenced <- sort(sample(c(1:ncol(B)),5))
#' B[activeSNPs[i], genes_influenced] <- 2
#' }
#' E <- matrix(rnorm(n*q),n,q)
#' Y <- X %*% B + L + E
#'
#' ## Usage
#' HC_Screening(Y, X)
#'
#'
HC_Screening <- function(Y, X){
  Bhat <- (InitialEst(Y, X))[[1]]
  Bhat_standard <- standardizeBhat(Y, X, Bhat)
  rHC <- rankHC(Bhat_standard)
  selectedSNPs <- sort(head(rHC$index,nrow(Y)))
  return(selectedSNPs)
}

#' Run_LORS
#'
#' \code{Run_LORS} is a function used to run either FastLORS or LORS
#'
#' @param Y gene expression matrix
#' @param X matrix of SNPs
#' @param method chooses with modeling method to run
#' @param screening Either "LORS-Screening", "HC-Screening", or "None".  The default method, LORS-Screening, is recommended if the number of SNPs is large. HC-Screening of Rhyne et al. (2018) is under development but is included here as an option.
#' @param tune_method chooses whether FastLORS should be used for parameter tuning or the original LORS procedure should be used.  Default is FastLORS
#' @param seed random seed to be used for setting training and validation set.  Default is 123.
#' @param cross_valid chooses whether cross-validation should be used in parameter tuning.  Default is TRUE.
#' @param maxiter maximum number of iterations
#' @param eps constant used when checking the convergence.  Ensures no division by 0.
#' @param tol tolerance level for convergence
#' @param omega_SOR the value of omega to use if applying successive over-relaxation with FastLORS.
#' @importFrom glmnet glmnet
#' @importFrom stats runif
#' @return \item{LORS_Obj or Fast_LORS_Obj}{A list produced from LORS or FastLORS containing (1) B: estimate of the coefficient matrix (2) L:  estimate of the matrix of hidden factors (3) mu: estiamte of the vector of intercepts (4) f_val_vec: objective function values and (5) res_vec: relative change in objective function values} \item{selectedSNPs}{The SNPs selected by the screening method} \item{screening_time}{The time (in seconds) spent on screening step} \item{param_time}{The time (in seconds) spent on the parameter tuning step} \item{model_time}{The time (in seconds) spent on the joint modeling step} \item{total_time}{The time (in seconds) spent on the screening, parameter tuning, and joint modeling steps}  \item{rho}{The value of rho chosen through parameter tuning} \item{lambda}{The value of lambda chosen through parameter tuning}
#' @export
#' @examples
#'
#' ##Example
#'
#' ## Generate some data
#' n <- 20
#' p <- 50
#' q <- 30
#' k <- 4
#' set.seed(123)
#' X <- matrix(rbinom(n*p,1,0.5),n,p)
#' L <- matrix(rnorm(n*k),n,k) %*% t(matrix(rnorm(q*k),q,k))
#' B <- matrix(0, ncol(X), ncol(L))
#' activeSNPs <- sort(sample(c(1:nrow(B)), 20))
#' for(i in 1:length(activeSNPs)){
#' genes_influenced <- sort(sample(c(1:ncol(B)),5))
#' B[activeSNPs[i], genes_influenced] <- 2
#' }
#' E <- matrix(rnorm(n*q),n,q)
#' Y <- X %*% B + L + E
#'
#' ## Usage
#' Run_LORS(Y, X, method = "FastLORS")

Run_LORS <- function(Y, X, method = "FastLORS", screening = "LORS-Screening", tune_method = "FastLORS", seed = 123,  maxiter = 10000, eps = 2.2204e-16, tol = 1e-4, cross_valid = TRUE, omega_SOR = 1.999){
  start <- proc.time()

  if(screening == "HC-Screening"){
    print("Beginning screening via HC-Screening")
    start_screening <- proc.time()
    selectedSNPs <- HC_Screening(Y, X)
    X <- X[,selectedSNPs]
    end_screening <- proc.time()
    screening_time <- end_screening[3] - start_screening[3]
  }

  if(screening == "LORS-Screening"){
    print("Beginning screening via LORS-Screening")
    start_screening <- proc.time()
    selectedSNPs <- Run_LORS_Screening(Y, X)
    X <- X[,selectedSNPs]
    end_screening <- proc.time()
    screening_time <- end_screening[3] - start_screening[3]
  }

  if(screening == "None"){
    selectedSNPs <- c(1:ncol(X))
    screening_time <- 0
  }

  n <- nrow(Y)
  q <- ncol(Y)
  p <- ncol(X)

  Omega0 <- !(is.na(Y))
  Y[is.na(Y)] <- 0

  # first use soft-impute to select a reasonable lambda
  set.seed(seed)
  mask = matrix(runif(n*q) > 0.5, nrow = n, ncol = q)
  Omega1 = Omega0 & mask
  Omega2 = Omega0 & !mask

  maxRank = min(n,q)/2
  print("Begin Parameter Tuning:  Performing two-fold cross validation")
  ParamTune <- function(Training, Validation, tune_method){
    mySI <- softImpute(Y,NULL,Omega0, Training, Validation, NULL,maxRank)
    Z <- mySI[["Z"]]
    Err <- mySI[["Err"]]
    rank_alpha <- mySI[["rank_alpha"]]
    Znorm <- mySI[["znorm"]]
    Alpha <- mySI[["Alpha"]]
    bestind_lam <- min(Err)==Err;

    lambda <- mean(Alpha[bestind_lam])

    ## Get a good initialization of L
    L <- SVT(Y,lambda)

    ## Set a sequence of rho

    nrho <- 20
    MaxRho <- GetMaxRho(X, Y, L, Omega0)

    rhoseq <- logspace(log10(MaxRho),log10(MaxRho*.05),nrho)

    rhoErr <- matrix(0, nrow = nrho, ncol = 1)
    B <- matrix(0, nrow = p, ncol = q)
    mu <- matrix(0, nrow = 1, ncol = q)
    for( irho in 1:nrho){
      print(paste("On rho",irho,"of 20"))
      rho <- rhoseq[irho]
      if(tune_method == "FastLORS"){
        myL2 <- Fast_LORS_Tuning(Y = Y, X = X, rho = rho, lambda = lambda, Training = Training, Validation = Validation, tol = tol, maxiter = maxiter, B = B, L = L, mu = mu, omega_SOR = omega_SOR)
      }
      if(tune_method == "LORS"){
        myL2 <- LORS2(Y = Y, X = X, L = L, Omega1 = Training, Omega2 = Validation, B = B, rho = rhoseq[irho], lambda = lambda, tol = tol, maxIter = maxiter)
      }
      B <- myL2[["B"]]
      mu <- myL2[["mu"]]
      L <- myL2[["L"]]
      Err <- myL2[["Err"]]
      rhoErr[irho,1] <- Err
    }

    ## use the best rho solve the optimization problem
    bestind_rho <- rhoErr==min(rhoErr)
    rho <- rhoseq[bestind_rho]
    return(list("rho" = rho, "lambda" = lambda))
  }

  start_param <- proc.time()
  print("Beginning parameter tuning:  Fold 1")
  params1 <- ParamTune(Training = Omega1, Validation = Omega2, tune_method)
  if(cross_valid == TRUE){
    print("Parameter tuning: Fold 2")
    params2 <- ParamTune(Training = Omega2, Validation = Omega1, tune_method)
  }
  if(cross_valid == FALSE){
    params2 <- params1
  }
  lambda <- mean(c(params1[["lambda"]], params2[["lambda"]]))
  rho <- mean(c(params1[["rho"]], params2[["rho"]]))
  end_param <- proc.time()
  param_time <- end_param[3] - start_param[3]

  if (method == "FastLORS"){
    print("Running FastLORS")
    start_model <- proc.time()
    Fast_LORS_Obj <- Fast_LORS(Y, X, rho, lambda, maxiter, eps, tol, omega_SOR = omega_SOR)
    end_model <- proc.time()
    model_time <- end_model[3] - start_model[3]
    end <- proc.time()
    total_time <- end[3] - start[3]
    return(list("Fast_LORS_Obj" = Fast_LORS_Obj, "selectedSNPs" = selectedSNPs,
                "screening_time" = screening_time, "param_time" = param_time,
                "model_time" = model_time, "total_time" = total_time, "rho" = rho, "lambda" = lambda))
  }

  if (method == "LORS"){
    print("Running LORS")
    start_model <- proc.time()
    LORS_Obj <- LORS0(Y, X, rho, lambda, maxiter, eps, tol)
    end_model <- proc.time()
    model_time <- end_model[3] - start_model[3]
    end <- proc.time()
    total_time <- end[3] - start[3]
    return(list("LORS_Obj" = LORS_Obj, "selectedSNPs" = selectedSNPs,
                "screening_time" = screening_time, "param_time" = param_time,
                "model_time" = model_time, "total_time" = total_time, "rho" = rho, "lambda" = lambda))
  }

}

#' Fast_LORS_Tuning
#'
#' \code{Fast_LORS_Tuning} is a function used perform parameter tuning using FastLORS instead of LORS
#'
#' @param Y gene expression matrix
#' @param X matrix of SNPs
#' @param rho parameter used to enforce sparsity of B
#' @param lambda parameter used to enforce low-rank structure of L
#' @param Training Boolean matrix for training data
#' @param Validation Boolean matrix for validation data
#' @param maxiter maximum number of iterates
#' @param eps a small constant to prevent dividing by zero when checking relative change in function values.
#' @param tol tolerance threshold for convergence
#' @param B an estimate for matrix of coefficients.  Default is NULL.
#' @param mu an estimate for the intercept.  Default is NULL.
#' @param L an estimate for the hidden factors  Default is NULL.
#' @param omega_SOR the value of omega to use if applying successive over-relaxation with FastLORS.
#' @importFrom glmnet glmnet
#' @importFrom methods as
#' @return \item{B}{matrix of coefficients} \item{L}{matrix of hidden factors} \item{mu}{vector of intercepts} \item{Err}{Residual sum of squares of validation data} \item{f_vals}{Objective function values} \item{res_vec}{Relative change in objective function values} \item{iter}{Number of iterations}
#' @export
#'

Fast_LORS_Tuning <- function(Y, X, rho, lambda, Training, Validation, maxiter = 5000, eps = 2.2204e-16, tol = 1e-4, B = NULL, mu = NULL, L = NULL, omega_SOR = 1.999) {

  ### Initial Setup

  n <- nrow(Y)
  q <- ncol(Y)
  p <- ncol(X)
  ones <- matrix(1, nrow = n, ncol = 1)

  if(is.null(B) == TRUE){
    B <- matrix(0, nrow = p, ncol = q)
  }
  if(is.null(mu) == TRUE){
    mu <- matrix(0, nrow = 1, ncol = q)
  }
  if(is.null(L) == TRUE){
    L <- matrix(0, nrow = n, ncol = q)
  }

  ### Proximal Mapping Functions

  ## Nuclear Norm:  Soft Thresholded SVD
  svd_st <- function(X, lambda){
    mysvd <- svd(X)
    U <- mysvd[["u"]]
    V <- mysvd[["v"]]
    VT <- t(V)
    D <- diag(mysvd[["d"]])
    d <- mysvd[["d"]]
    index_list <- which(d >= rep(lambda, length(d)))
    if(length(index_list) > 1){
      W <- diag(c(d[index_list] - lambda))
      L <- U[,index_list] %*% W %*% VT[index_list,]
    }
    if(length(index_list) == 1){
      W <- matrix(d[index_list] - lambda, nrow = 1, ncol = 1)
      L <- U[,index_list] %*% W %*% VT[index_list,]
    }
    if(length(index_list) == 0){
      L <- matrix(0, nrow(X), ncol(X))
    }

    return(L)
  }

  ## 1 norm:  Soft Threshold Function
  prox_1 <- function(b, tau){
    prox_b <- sign(b)*(sapply(abs(b) - tau, FUN=function(x) {max(x,0)}))
    return(prox_b)
  }

  fval_old <- 0
  f_val_vec <- c()
  res_vec <- c()

  t_L <- 1
  t_B <- 1/norm(t(X) %*% X, type = "2")
  t_mu <- 1/nrow(Y)

  ones <- matrix(1, nrow = n, ncol = 1)

  for(iter in 1:maxiter){

    #### Update B, mu, and L
    L <- svd_st(L - t_L * (Training * (X %*% B + ones %*% mu + L - Y)), t_L * lambda)
    B <- prox_1(B - t_B * t(X) %*% (Training * (X %*% B + ones %*% mu + L - Y)), t_B * rho)
    mu <- mu - t_mu * t(ones) %*% (Training * (X %*% B + ones %*% mu + L - Y))

    ## Update via SOR

    if(iter > 1){
      L_update <- (1-omega_SOR) * L_update + omega_SOR * L
      B_update <- (1-omega_SOR) * B_update + omega_SOR * B
      mu_update <- (1-omega_SOR) * mu_update + omega_SOR * mu

      L <- L_update
      B <- B_update
      mu <- mu_update
    }

    if(iter == 1){
      L_update <- L
      B_update <- B
      mu_update <- mu

      L <- L_update
      B <- B_update
      mu <- mu_update
    }

    #### Check Convergence

    dum <- c(Training * (Y - X %*% B - ones%*%mu - L))
    fval <- 0.5 * norm(dum, type = "2")^2 + rho * sum(abs(B)) + lambda*sum(svd(L)[["d"]])

    res = abs(fval-fval_old)/abs(fval_old+eps)

    fval_old <- fval
    f_val_vec <- c(f_val_vec, fval)
    res_vec <- c(res_vec, res)

    #print(paste('Iter ', iter, 'fval', fval, 'res', res))

    if (res < tol){
      break
    }

    if(iter > 1){
      if(f_val_vec[iter] > f_val_vec[iter-1]){
        break
      }
    }

    residual = Y - X%*%B - matrix(1, nrow = n, ncol = 1) %*% mu - L
    err = residual*Validation
    Err = t(c(err)) %*% c(err) / sum(Validation)

  }

  if(iter < maxiter & res >= tol){
    #print("Beginning LORS Updates")
    energy = Inf

    for (iter_LORS in (iter + 1):maxiter){

      for (innerIts in 1:50){

        Z = Y - X %*% B - matrix(1, nrow = n, ncol = 1) %*% mu;
        C = Z*Training + L*(!Training)

        U = svd(C)[["u"]]
        D = diag(svd(C)[["d"]])
        V = svd(C)[["v"]]
        VT = t(V)
        # soft impute
        d = diag(D)

        index_list = which(d > lambda)

        W = diag(d[index_list] - lambda, nrow = length(index_list), ncol = length(index_list))
        L = U[,index_list] %*% W %*% VT[index_list,]

        Lnorm = sum(d[index_list] - lambda)
        energy_old = energy

        cL = c(L)
        cZ = c(Z)
        newvec = matrix(cL[c(Training == 1)]-cZ[c(Training == 1)], ncol = 1)
        energy = lambda * Lnorm + norm(newvec,'F')/2

        mydiff = abs(energy - energy_old) / energy_old

        if (!is.na(mydiff)){
          if (abs(energy - energy_old) / energy_old < tol){
            break
          }
        }
      }

      ### Compute B using glmnet

      for (j in 1:q){
        fit <- glmnet(X[Training[,j],], matrix(Y[Training[,j],j]) - matrix(L[Training[,j],j]), family = "gaussian", lambda = rho/sum(Training[,j]), standardize = FALSE)
        a0 <- fit[["a0"]]
        old_beta <- fit[["beta"]]
        my_beta <- as(old_beta, "matrix")
        B[,j] <- my_beta
        mu[,j] <- a0
      }

      # Convergence
      residual = Y - X%*%B - matrix(1, nrow = n, ncol = 1) %*% mu - L
      dum = residual*Training
      dum = c(dum)
      fval = 0.5 * t(dum) %*% dum + rho * sum(abs(c(B))) + lambda*sum(abs(diag(W)))
      res = abs(fval-fval_old)/abs(fval_old+eps)

      fval_old <- fval
      f_val_vec <- c(f_val_vec, fval)
      res_vec <- c(res_vec, res)

      #print(paste('Iter ', iter_LORS, 'fval', fval, 'res', res))

      if (res < tol){
        break
      }

    }

    err = residual*Validation
    Err = t(c(err)) %*% c(err) / sum(Validation)
  }

  ### Return the B, L, mu, objective function values, residual values, and total iterates
  return(list("B" = B, "L" = L, "mu" = mu, "Err" = Err, "f_vals" = f_val_vec, "res_vec" = res_vec, "iter" = iter))
}




#' Run_LORS_Screening
#' \code{Run_LORS_Screening} is a function to to apply the LORS-Screening Algorithm in Yang et al. (2013)
#'
#' @param Y gene expression matrix
#' @param X matrix of SNPs
#' @param lambda tuning parameter
#' @importFrom stats runif
#' @return selectedSNPs the SNPs selected by LORS-Screening
#' @export
#' @examples
#'
#' ##Example
#'
#' ## Generate some data
#' n <- 50
#' p <- 200
#' q <- 100
#' k <- 10
#' set.seed(123)
#' X <- matrix(rbinom(n*p,1,0.5),n,p)
#' L <- matrix(rnorm(n*k),n,k) %*% t(matrix(rnorm(q*k),q,k))
#' B <- matrix(0, ncol(X), ncol(L))
#' activeSNPs <- sort(sample(c(1:nrow(B)), 20))
#' for(i in 1:length(activeSNPs)){
#' genes_influenced <- sort(sample(c(1:ncol(B)),5))
#' B[activeSNPs[i], genes_influenced] <- 2
#' }
#' E <- matrix(rnorm(n*q),n,q)
#' Y <- X %*% B + L + E
#'
#' selectedSNPs <- Run_LORS_Screening(Y, X)

Run_LORS_Screening <- function(Y, X, lambda = NULL){

  n = nrow(Y)
  q = ncol(Y)
  p = ncol(X)

  Omega0 = !(is.na(Y))
  Y[is.na(Y)] <- 0

  if(is.null(lambda) == TRUE){
    # first use soft-impute to select a reasonable lambda
    set.seed(123456789)
    mask = matrix(runif(n*q) > 0.5, nrow = n, ncol = q)
    Omega1 = Omega0 & mask
    Omega2 = Omega0 & !mask

    maxRank = min(n,q)/2
    mySI= softImpute(Y,NULL,Omega0, Omega1, Omega2, NULL,maxRank)
    Z = mySI[["Z"]]
    Err = mySI[["Err"]]
    rank_alpha = mySI[["rank_alpha"]]
    Znorm = mySI[["znorm"]]
    Alpha = mySI[["Alpha"]]
    bestind_lam = min(Err)==Err;

    lambda = min(Alpha[bestind_lam])
  }

  myB <- c()
  print("Building Initial Estimate")
  for(SNP_col in 1:ncol(X)){
    #print(paste("Building initial estimate: On column", SNP_col,"of",ncol(X)))
    X1 <- matrix(X[,SNP_col], ncol = 1)
    LS <- LORSscreen(Y, X1, lambda, 0.01)
    B_row <- LS$B
    myB <- rbind(myB, B_row)
  }

  index_list = c()
  for (i in 1:ncol(myB)){
    cands = abs(myB[,i])
    sorted_cands <- sort(cands, index.return=TRUE, decreasing=TRUE)$ix[1:nrow(X)]
    index_list <- c(index_list, sorted_cands)
  }
  selectedSNPs <- sort(unique(index_list))
  return(selectedSNPs)
}

#' ParamTuneParallel
#' \code{ParamTuneParallel} is a function used to run the parameter tuning of either FastLORS or LORS in parallel
#'
#' @param Y gene expression matrix
#' @param X matrix of SNPs
#' @param fold fold in cross-validation
#' @param seed random seed used to create training and validation sets
#' @importFrom glmnet glmnet
#' @importFrom stats runif
#' @return \item{myParams}{lambda and a sequence of rho values to use in parameter tuning} \item{Training}{Training Data} \item{Validation}{Validation Data}
#' @export
#' @examples
#'
#' ##Example
#'
#' ## Generate some data
#' n <- 20
#' p <- 50
#' q <- 30
#' k <- 4
#' set.seed(123)
#' X <- matrix(rbinom(n*p,1,0.5),n,p)
#' L <- matrix(rnorm(n*k),n,k) %*% t(matrix(rnorm(q*k),q,k))
#' B <- matrix(0, ncol(X), ncol(L))
#' activeSNPs <- sort(sample(c(1:nrow(B)), 20))
#' for(i in 1:length(activeSNPs)){
#' genes_influenced <- sort(sample(c(1:ncol(B)),5))
#' B[activeSNPs[i], genes_influenced] <- 2
#' }
#' E <- matrix(rnorm(n*q),n,q)
#' Y <- X %*% B + L + E
#'
#' ## Usage
#' ParamTuneParallel(Y, X, fold = 1)

ParamTuneParallel <- function(Y, X, fold, seed = 123){

  n <- nrow(Y)
  q <- ncol(Y)
  p <- ncol(X)

  Omega0 <- !(is.na(Y))
  Y[is.na(Y)] <- 0

  # first use soft-impute to select a reasonable lambda
  set.seed(seed)
  mask = matrix(runif(n*q) > 0.5, nrow = n, ncol = q)
  Omega1 = Omega0 & mask
  Omega2 = Omega0 & !mask

  maxRank = min(n,q)/2
  print("Begin Parameter Tuning:  Performing two-fold cross validation")
  myParamTune <- function(Training, Validation, tune_method){
    mySI <- softImpute(Y,NULL,Omega0, Training, Validation, NULL,maxRank)
    Z <- mySI[["Z"]]
    Err <- mySI[["Err"]]
    rank_alpha <- mySI[["rank_alpha"]]
    Znorm <- mySI[["znorm"]]
    Alpha <- mySI[["Alpha"]]
    bestind_lam <- min(Err)==Err;

    lambda <- mean(Alpha[bestind_lam])

    ## Get a good initialization of L
    L <- SVT(Y,lambda)

    ## Set a sequence of rho

    nrho <- 20
    MaxRho <- GetMaxRho(X, Y, L, Omega0)

    rhoseq <- logspace(log10(MaxRho),log10(MaxRho*.05),nrho)
    return(list("lambda" = lambda, "rhoseq" = rhoseq))
  }
  if(fold == 1){
    Training <- Omega1
    Validation <- Omega2
    myParams <- myParamTune(Training, Validation)
  }
  if(fold == 2){
    Training <- Omega2
    Validation <- Omega1
    myParams <- myParamTune(Training, Validation)
  }
  return(list("myParams" = myParams, "Training" = Training, "Validation" = Validation))
}

#' LORS_Screen_Parallel
#' \code{LORS_Screen_Parallel} is a function used to run LORS-Screening on a subset of the columns of X.  Can be used to perform LORS-Screening in parallel on a cluster.
#'
#' @param Y gene expression matrix
#' @param X matrix of SNPs
#' @param chunk a group of columns to run the screening on.  Done in batches of 1000.
#' @importFrom glmnet glmnet
#' @importFrom stats runif
#' @return \item{myB}{matrix of coefficients from LORS-Screening} \item{lambda}{tuning parameter used in LORS-Screening}
#' @export
#' @examples
#'
#' ##Example
#'
#' ## Generate some data
#' n <- 20
#' p <- 50
#' q <- 30
#' k <- 4
#' set.seed(123)
#' X <- matrix(rbinom(n*p,1,0.5),n,p)
#' L <- matrix(rnorm(n*k),n,k) %*% t(matrix(rnorm(q*k),q,k))
#' B <- matrix(0, ncol(X), ncol(L))
#' activeSNPs <- sort(sample(c(1:nrow(B)), 20))
#' for(i in 1:length(activeSNPs)){
#' genes_influenced <- sort(sample(c(1:ncol(B)),5))
#' B[activeSNPs[i], genes_influenced] <- 2
#' }
#' E <- matrix(rnorm(n*q),n,q)
#' Y <- X %*% B + L + E
#'
#' ## Usage
#' LORS_Screen_Parallel(Y, X, chunk = 1)

LORS_Screen_Parallel <- function(Y, X, chunk){
  start_col <- (chunk - 1) * 1000 + 1
  end_col <- min(c(ncol(X), (chunk - 1) * 1000 + 1000))

  n = nrow(Y)
  q = ncol(Y)
  p = ncol(X)

  Omega0 = !(is.na(Y))
  Y[is.na(Y)] <- 0

  # first use soft-impute to select a reasonable lambda
  set.seed(123456789)
  mask = matrix(runif(n*q) > 0.5, nrow = n, ncol = q)
  Omega1 = Omega0 & mask
  Omega2 = Omega0 & !mask

  maxRank = min(n,q)/2
  set.seed(123456789)
  mySI= softImpute(Y,NULL,Omega0, Omega1, Omega2, NULL,maxRank)
  Z = mySI[["Z"]]
  Err = mySI[["Err"]]
  rank_alpha = mySI[["rank_alpha"]]
  Znorm = mySI[["znorm"]]
  Alpha = mySI[["Alpha"]]
  bestind_lam = min(Err)==Err;

  lambda = min(Alpha[bestind_lam])

  myB <- c()
  print("Building Initial Estimate")
  for(SNP_col in start_col:end_col){
    print(paste("Building initial estimate: On column", SNP_col,"of",end_col))
    X1 <- matrix(X[,SNP_col], ncol = 1)
    LS <- LORSscreen(Y, X1, lambda, 0.01)
    B_row <- LS$B
    myB <- rbind(myB, B_row)
  }
  mylist <- list("myB" = myB, "lambda" = lambda)
  return(mylist)
}
