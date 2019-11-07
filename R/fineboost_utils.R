######### CCG : Causal confidence grade is an alternative to Posterior inclusion probability (PIP)  #######

#' @title Compute inclusion confidence grade (CCG) for all variables
#' @param ff an object of class fineboos, generated from `fineboostR::fineboost_smooth()`
#' @param use_csets Whether to use the confidence sets when computing the ICG. Defaults to TRUE.
#' @return a vector inclusion confidence grade for all variables.
#' @export
#'
fineboost_get_ccg <- function (ff, use_csets = TRUE){
  if (class(ff) == "fineboost"){
    if(!use_csets){
      ccg = rep(0, ff$P)
      ccg[which(ff$csets$csets_index > 0)] = apply(ff$csets$cluster_signal[, which(ff$csets$csets_index > 0)],
                                                   2, function(x) return(max(x)))
      return(ccg)
    }
    if(use_csets){
      num_sets = max(ff$csets$csets_index)
      ccg = rep(0, length(ff$csets$csets_index))
      for(jj in 1:num_sets){
        idx=which(ff$csets$csets_index == jj)
        ccg[idx] = apply(ff$weights_path[,idx, drop=F], 2, function(x) return(max(x)))
      }
      return(ccg)
    }
  }
}


######### Confidence sets generated from the output of Fineboost  #######

#' @title Get confidence sets for local signal clusters from FineBoost output
#' @param ff an object of class fineboos, generated from `fineboostR::fineboost_smooth()`
#' @param LD the LD P times P matrix, where P is the number of features, equal to ff$P,
#' @param coverage A number between 0 and 1 (close to 1) specifying the coverage of the estimated signal clusters.
#'                 Default set to 0.95.
#' @param clus_thresh A number between 0 and 1 (close to 0) that is used to filter out local signal clusters with
#'                    depleted number of boosting iterations and high level of uniformity of signal. Default is
#'                    set to 0.1.
#' @param nmf_try Number of initializations of the Non-negative matrix factorization on the trajectory of weights
#'                tested. Default is set to 10.
#' @param min_abs_corr Minimum of absolute value of correlation allowed in a credible set. The default, 0.5,
#'                     corresponds to squared correlation of 0.25, which is a commonly used threshold for
#'                     genotype data in genetics studies.
#' @return An object with the following items
#'   \item{evidence}{A vector of length `Lmax` giving the evidence strength for each of signal clusters.}
#'   \item{unfiltered_sets}{A list of length `Lmax` denoting the confidence sets that contain causal signal, where the
#'                         strength of confidence is provided by the coverage.}
#'   \item{cluster_signal}{The cluster signal for the `unfiltered sets`.}
#'   \item{filtered_sets}{The threshold used to filter away cluster signals of the evidence for that cluster falls below this
#'                        threshold. Default is set at 0.1.}
#'   \item{csets_index}{A vector of length P that labels each variable by the confidence set in the `filtered_sets` it
#'                      belongs to and 0 otherwise}
#' @export
#'
fineboost_get_csets <- function (ff,
                                 coverage = 0.95,
                                 clus_thresh = 0.1,
                                 nmf_try = 10){
  if (class(ff) == "fineboost"){
    vx = vector(mode="list", nmf_try)
    for(rr in 1:nmf_try){
      Lmax = ff$Lmax
      set.seed(rr)
      res = NNLM::nnmf(ff$weights_path, k=Lmax)
      W2 = t(apply(res$W, 1, function(x) return(x/sum(x))))
      clus_med = c()
      presence_prop = c()
      for(k in 1:Lmax){
        idx = which(W2[,k] > 0.5)
        if(length(idx) > 0){
          clus_med = rbind(clus_med, apply(ff$weights_path[idx, , drop=F], 2, median))
        }
        presence_prop = c(presence_prop, length(idx)/nrow(W2))
      }

      confidence_sets = list()
      for(mm in 1:nrow(clus_med)){
        tempp = order(clus_med[mm, ], decreasing=T)
        confidence_sets[[mm]] = tempp[1:min(which(cumsum(clus_med[mm,tempp]) > coverage))]
      }

      evidence_strength = presence_prop
      confidence_sets_1 = confidence_sets[order(presence_prop,
                                                decreasing = T)]
      confidence_sets_2 = confidence_sets[order(presence_prop,
                                                decreasing = T)[1:length(which(evidence_strength > clus_thresh))]]
      csets_index = rep(0, ff$P)
      for(kk in 1:length(confidence_sets_2)){
        csets_index[confidence_sets_2[[kk]]] = kk
      }

      ll = list("evidence" = evidence_strength,
                "unfiltered_sets" = confidence_sets_1,
                "cluster_signal" = clus_med,
                "filtered_sets" = confidence_sets_2,
                "csets_index" = csets_index)
      vx[[rr]] = ll
    }

    temp = matrix(0, nmf_try, ff$P)
    for(rr in 1:nmf_try){
      temp[rr, which(vx[[rr]]$csets_index != 0)] = 1
    }
    best_run = which.max(rowSums(cor(t(temp))))[1]
    ll2 = vx[[best_run]]
    return(ll2)
  }
}

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
update_kernel_weights <- function(X, y, LD, tau){

  cor_vals = apply(X, 2, function(z) return(cor(z, y))) ## corr. of X.columns and z

  abs_cor_vals = abs(cor_vals) ## abs. values of correlations

  idx = which(abs_cor_vals == max(abs_cor_vals))[1] ## which top feature that has highest correlation

  dist = sqrt(abs(LD[,idx])) ## overall kernel intensities around the top feature

  exp_abs_cor_vals = exp(abs_cor_vals*dist/tau)+1e-10 ## optimal kernel smoothing

  if(any(is.infinite(exp_abs_cor_vals))){
    stop("weights out of bounds, please choose a larger tau")
  }

  weights= exp_abs_cor_vals/sum(exp_abs_cor_vals) ## updated weight at this iteration step

  current_obj = sum(weights*dist*abs_cor_vals) ## objective function

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
