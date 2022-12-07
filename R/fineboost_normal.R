#' @title A genetic fine-mapping method using kernel-based Forward Stepwise boosted regression model.
#'
#' @description This function uses a kernel-based FS-boost framework to find causal fine-mapped SNP sets
#' from GWAS or QTL effect size data. It performs a regression  \eqn{Y = Xb + e} where b is a sparse vector
#' of coefficients with local signal clusters. Unlike other fine-mapping methods, our algorithm is model-free
#' on the coefficients b.
#'
#' @param X   The design matrix X (N times P) with samples/individuals along the rows
#'            and putatively correlated ordered features (SNPs) along the columns.
#'
#' @param Y  The response vector of length N
#'
#' @param Lmax The maximum number of local signal clusters fitted.
#'
#' @param M  The maximum number of boosting iterations to run. Default is 1000.
#'
#' @param LD The external LD matrix for the P features of interest. Defaults to NULL, in which case, in-sample LD is used.
#'
#' @param step The stepsize used in boosting iterations. Default set to 0.05.
#'
#' @param kern_tau The smoothing intensity of the kernel averaging at each boosting iteration. Default set to 0.01.
#'
#' @param method The boosting update method- either `LS` or `FS` indicating the LS-Boost and FS-epsilon methods respectively.
#'               Default is set to LS-Boost.
#'
#' @param kernel The nature of the kernel used for smoothing. Can be either `L1`, `L2`, `epanechnikov` or `prune`. `L1`
#'               kernel uses a L-1 norm based kernel, `L2` uses a L-2 norm based kernel, `epanechnikov` uses an
#'               Epanechnikov kernel and `prune` uses a uniform kernel on all SNPs with high LD to the optimal SNP at
#'               each boosting iteration.
#'
#' @param stop_thresh The stopping threshold (small number) for the objective function, when attained,
#'                    the boosting iterations will stop automatically. Default is 0.1.
#'
#' @param na.rm Drop missing samples in y from both y and X inputs. Default set to FALSE.
#'
#' @param intercept Boolean; if there is an intercept in the model to fit. Defaults to TRUE.
#'
#' @param standardize Boolean; if the columns of X need to be standardized. Defaults to TRUE.
#'
#' @param coverage A number between 0 and 1 (close to 1) specifying the coverage of the estimated signal clusters.
#'                 Default set to 0.95.
#'
#' @param min_within_LD The minimum value of LD permitted for SNPs within a local signal cluster. Default is 0.25.
#' @param min_between_LD The minimum value of LD permitted for SNPs across two local signal clusters. Default is it
#'                       cannot exceed 0.25.
#' @param min_cluster_centrality The minimum value of cluster centrailty required for a SNP to make the cut in a
#'                      local signal cluster. Default is set at 0.5.
#'
#' @param nmf_try The number of NMF initiializations to fix the confidence sets. Default is set to 5.
#'
#' @param min_abs_corr Minimum of absolute value of correlation allowed in a credible set. The default, 0.5,
#'                     corresponds to squared correlation of 0.25, which is a commonly used threshold for
#'                     genotype data in genetics studies.
#'
#' @param verbose If \code{verbose = TRUE}, information about the objective and progress at each iteration of the
#'                  kerne-based boosting procedure is returned.
#'
#' @return A \code{"fineboost"} object with the following elements:
#'    \item{N} The number of rows of the input matrix X.
#'    \item{P} The number of columns of the input matrix X.
#'    \item{Lmax} The maximum number of signal clusters to be allowed when postprocessing the FINEBOOST model.
#'    \item{beta} The final estimated coefficient vector b in the regression  \eqn{Y = Xb + e}.
#'    \item{beta_path} The trajectory of all estimates of beta from iteration 1 till stoppage.
#'    \item{weights_path} The trajectory of all kernel smoothing weights used in updated from iteration 1 till stoppage.
#'    \item{profile_loglik} The profile loglikelihood of the successive iterations until stoppage.
#'    \item{obj_path} the trajectory of the actual objective to minimize from iteration 1 till stoppage.
#'    \item{csets} The confidence sets (of specified input `coverage`) for signal clusters
#'
#' @export




fineboost_normal <- function(X, Y, M=1000,
                             Lmax=5, LD = NULL,
                             step = 0.1, kern_tau = 0.01,
                             method = "LS",
                             kernel = "L2",
                             stop_thresh = 1e-04, na.rm=FALSE,
                             intercept=TRUE, standardize=TRUE,
                             coverage = 0.95, 
                             min_within_LD = 0.5,
                             min_between_LD = 0.25,
                             min_clus_centrality = 0.5,
                             nmf_try = 5, verbose=TRUE){

  if (!(is.double(X) & is.matrix(X)) & !inherits(X,"CsparseMatrix") & is.null(attr(X,"matrix.type")))
    stop("Input X must be a double-precision matrix, or a sparse matrix, or a trend filtering matrix.")
  if (any(is.na(X))) {
    stop("Input X must not contain missing values.")
  }
  if (any(is.na(Y))) {
    if (na.rm) {
      samples_kept = which(!is.na(Y))
      Y = Y[samples_kept]
      X = X[samples_kept,]
    } else {
      stop("Input Y must not contain missing values.")
    }
  }
  p = ncol(X)
  n = nrow(X)
  mean_y = mean(Y)

  if(intercept){
    Y = Y-mean_y
  }
  X = set_X_attributes(X,center=intercept, scale=standardize)



  ##################  initialization #######################################

  ff = list("N" = nrow(X),
            "P" = ncol(X),
            "Lmax" = Lmax,
            "beta" = rep(0, dim(X)[2]),
            "beta_path" = c(),
            "weights_path" = c(),
            "profile_loglik" = c(),
            "obj_path" = 9999999)
  class(ff) = "fineboost"

  current_obj=10^4
  mm=1
  res = Y
  step1 = step

  ##################  Fineboost updates   ###################################

  for(m in 1:M){

    newll = update_kernel_weights(attr(X, "scaled"), res, attr(X, "LD"), tau = kern_tau, kernel = kernel)

    ff$obj_path = c(ff$obj_path, newll$objective)

    ff$weights_path = rbind(ff$weights_path, newll$weights)

    beta_grad = newll$weights*sign(newll$corvals) ## gradient of objective function wrt b

    prediction = (attr(X, "scaled") %*% (beta_grad))

    if(method == "LS"){
      if(newll$objective > 0.5){
        step1 = calc_step(prediction, res, step1)
      }else{
        step1 = step
      }
    }else if(method == "FS"){
      if(newll$objective > 0.5){
       step1 = 0.5
      }else{
       step1 = step
      }
    }else if(method != 'FS'){
      stop("The boosting method has to be either LS or FS; denoting LS-Boost and FS-epsilon type updates.")
    }

    res = res - step1*prediction   ## Gradient Descent on residuals

    ff$beta = ff$beta + step1*beta_grad ## Gradient ascent on b

    ff$beta_path = rbind(ff$beta_path, ff$beta)

    ff$profile_loglik = c(ff$profile_loglik, mean((Y - attr(X, "scaled")%*%ff$beta)^2))
    #ff$profile_loglik = c(ff$profile_loglik, max(abs(Y - attr(X, "scaled")%*%ff$beta)))

    if(((ff$obj_path[length(ff$obj_path)-1] - ff$obj_path[length(ff$obj_path)] < stop_thresh &&
        newll$objective < 0.1) || (ff$profile_loglik[length(ff$profile_loglik)-1] -
                                   ff$profile_loglik[length(ff$profile_loglik)] < stop_thresh &&
                                   ff$profile_loglik[length(ff$profile_loglik)]/ff$profile_loglik[1] < 1e-02)))
    {
      cat("Boosting iterations converge after", m, "iterations! \n")
      break;
    }

    if(verbose){
      cat(paste0("objective:", newll$objective, "at iteration:", m, "\n"))
    }
  }

  ff$num_updates = m
  ff$obj_path = ff$obj_path[-1]

  if(m == M){
    warning(paste("Fineboost updates did not converge in", M, "iterations; checkpoint at last iteration"))
    ff$converged = FALSE
  }

  #############  Post-processing of the Fineboost updates  #####################

  ff$csets = fineboost_get_csets(ff, X, coverage=coverage,
                                 min_within_LD = min_within_LD,
                                 min_between_LD = min_between_LD,
                                 min_clus_centrality = min_clus_centrality,
                                 nmf_try = nmf_try)
  ff$ccg =  fineboost_get_ccg(ff, use_csets = FALSE)
  ff$beta = fineboost_get_coef(X, Y, ff)
  return(ff)
}
