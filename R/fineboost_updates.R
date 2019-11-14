
######### Calculate the parameter tau for the kernel smoothing  ###########################

#' @title Exploring results for different values of kernel smoothing parameter
#' @param X   The design matrix X (N times P) with samples/individuals along the rows
#'            and putatively correlated ordered features (SNPs) along the columns.
#' @param Y  The response vector of length N
#' @param maxsize the maximum number of non-zero causal loci based on L0Learn.
#' @return The objective function calculation for different choices of tau.
#' @export
#'
calc_tau = function(X, Y, maxsize=10){
  fit = L0Learn.cvfit(X, Y, maxSuppSize = maxsize, intercept = F)
  lambda = fit$fit$lambda[[1]][which(fit$fit$suppSize[[1]] == maxsize)]
  gamma=0
  betaval = coef(fit, lambda, gamma)
  res = Y - X%*%as.numeric(betaval)
  prod = t(X)%*%res

  tau_func <- function(tau){
    tau*logSumExp(prod/tau) - tau*log(length(prod))
  }
  cat("tau 1e-05:", tau_func(1e-05), ",",
      "tau 0.01:", tau_func(0.01), ",",
      "tau 0.05:", tau_func(0.05), ",",
      "tau 0.5:", tau_func(0.5), "\n")
}


######### A function to update the coefficient vector based on kernel smoothing  #######

#' @title One instance of the kernel-smoothed FS-boost weights update
#' @param X   The design matrix X (N times P) with samples/individuals along the rows
#'            and putatively correlated ordered features (SNPs) along the columns.
#' @param y  The response vector of length N
#' @param LD The LD matrix for all P variables in the X matrix.
#' @param tau The kernel smoothing bandwidth estimate.
#' @return An object with the following items
#'   \item{weights}{A new vector of weights of length P based on kernel smoothing around the top feature.}
#'   \item{objective}{The current value of the objective to be minimized.}
#'   \item{corvals}{The current values of correlations of the P features with the current residual.}
#' @export
#'
update_kernel_weights <- function(X, y, LD, tau, kernel = "L2"){

  cor_vals = apply(X, 2, function(z) return(cor(z, y))) ## corr. of X.columns and z

  abs_cor_vals = abs(cor_vals) ## abs. values of correlations

  idx = which(abs_cor_vals == max(abs_cor_vals))[1] ## which top feature that has highest correlation

  if(kernel == "L2"){
    dist = sqrt(abs(LD[,idx])) ## overall kernel intensities around the top feature
    exp_abs_cor_vals = exp(abs_cor_vals*dist/tau)+1e-10 ## optimal kernel smoothing
  }else if (kernel == "L1"){
    dist = abs(LD[,idx]) ## overall kernel intensities around the top feature
    exp_abs_cor_vals = exp(abs_cor_vals*dist/tau)+1e-10 ## optimal kernel smoothing
  }else if (kernel == "epanechnikov"){
    dist = 0.75*(1 - (1 - abs(LD[,idx]))^2)
    exp_abs_cor_vals = exp(abs_cor_vals*dist/tau)+1e-10 ## optimal kernel smoothing
  }else if (kernel == "prune"){
    dist = rep(0, nrow(LD))
    dist[which(LD[,idx] > 0.5)] = 1
    exp_abs_cor_vals = dist*exp(abs_cor_vals/tau)+1e-10 ## optimal kernel smoothing
  }

  if(any(is.infinite(exp_abs_cor_vals))){
    stop("weights out of bounds, please choose a larger tau")
  }

  weights= exp_abs_cor_vals/sum(exp_abs_cor_vals) ## updated weight at this iteration step

  current_obj = tau*matrixStats::logSumExp(abs_cor_vals/tau) - tau*log(ncol(X)) ## objective function

  ll = list("weights" = weights,
            "kernel_wt" = dist,
            "objective" = current_obj,
            "corvals" = cor_vals)
  return(ll)
}


######### A function to update the coefficient vector without kernel smoothing  #######

#' @title One instance of the no kernel smoothed general FS-boost weights update
#' @param X   The design matrix X (N times P) with samples/individuals along the rows
#'            and putatively correlated ordered features (SNPs) along the columns.
#' @param y  The response vector of length N
#' @return An object with the following items
#'   \item{weights}{A new vector of weights of length P based on kernel smoothing around the top feature.}
#'   \item{objective}{The current value of the objective to be minimized.}
#'   \item{corvals}{The current values of correlations of the P features with the current residual.}
#' @export
#'
update_non_kernel_weights <- function(X, y){

  cor_vals = apply(X, 2, function(z) return(cor(z, y))) ## corr. of X.columns and z

  abs_cor_vals = abs(cor_vals) ## abs. values of correlations

  idx = which(abs_cor_vals == max(abs_cor_vals))[1] ## which top feature that has highest correlation

  weights = rep(0, ncol(X))
  weights[idx] = 1 ## updated weight at this iteration step

  current_obj = sum(weights*abs_cor_vals) ## objective function

  ll = list("weights" = weights,
            "objective" = current_obj,
            "corvals" = cor_vals)
  return(ll)
}


######### Calculate the adaptive step size using line search at each update  ###########################

#' @title Get optimal step sizes for boosting updates
#' @param xx A vector of length N (arbitrary).
#' @param yy A second vector of length N.
#' @param curr_step The stepsize at the start of the update.
#' @return The optimal step size for the next update of the coefficients.
#' @export
#'
calc_step = function(xx,yy, curr_step){
  epsilon_range = c(0.01, 0.05, 0.1, 1, 10, 100)
  tmp = sapply(epsilon_range, function(z) return(sqrt(sum((yy-z*xx)^2))))
  best_eps = epsilon_range[which.min(tmp)]
  if(best_eps != curr_step){
    cat("Changed the step-size to:", best_eps, "\n")
  }
  return(best_eps)
}






