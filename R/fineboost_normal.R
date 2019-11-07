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
#' @param clus_thresh A number between 0 and 1 (close to 0) that is used to filter out local signal clusters with
#'                    depleted number of boosting iterations and high level of uniformity of signal. Default is
#'                    set to 0.1.
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
                             step = 0.05, kern_tau = 0.01,
                             stop_thresh = 1e-04, na.rm=FALSE,
                             intercept=TRUE, standardize=TRUE,
                             coverage = 0.95, clus_thresh=0.1,
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
            "beta" = rep(0, dim(Xmat)[2]),
            "beta_path" = c(),
            "weights_path" = c(),
            "profile_loglik" = c(),
            "obj_path" = 9999999)
  class(ff) = "fineboost"

  current_obj=10^4
  mm=1
  res = Y

  ##################  Fineboost updates   ###################################

  for(m in 1:M){

    newll = update_kernel_weights(attr(X, "scaledX"), res, attr(X, "LD"), tau = kern_tau)

    ff$obj_path = c(ff$obj_path, newll$objective)

    if(ff$obj_path[length(ff$obj_path)] > ff$obj_path[length(ff$obj_path) - 1]+ 0.1){
      stop("objective value increases over iterations; possible model mismatch")
    }

    if(ff$obj_path[length(ff$obj_path)-1] - ff$obj_path[length(ff$obj_path)] < stop_thresh){
      cat("Boosting iterations converge after", m, "iterations! \n")
      break;
    }

    ff$weights_path = rbind(ff$weights_path, newll$weights)

    beta_grad = newll$weights*sign(newll$corvals) ## gradient of objective function wrt b

    res = res - step*(attr(X, "scaledX") %*% beta_grad) ## Gradient Descent on residuals

    ff$beta = ff$beta + step*beta_grad ## Gradient ascent on b

    ff$beta_path = rbind(ff$beta_path, ff$beta)

    ff$profile_loglik = c(ff$profile_loglik, sum((Y - attr(X, "scaledX")%*%ff$beta)^2))

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

  ff$csets = fineboost_get_csets(ff, coverage=coverage, clus_thresh = clus_thresh, nmf_try = nmf_try)
  ff$ccg = fineboost_get_ccg(ff, use_csets = FALSE)
  return(ff)
}
