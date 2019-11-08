######### A function to generate the final values of coefficients in the FINEBoost model  #######

#' @title Generate coefficient estimates using the FINEBoost model
#' @param X   The design matrix X (N times P) with samples/individuals along the rows
#'            and putatively correlated ordered features (SNPs) along the columns.
#' @param y  The response vector of length N
#' @param ff  An object of class fineboost, containing the fitted model with ccg and confidence sets information.
#' @return Returns a P dimensional vector b in \eqn{y = Xb + e} based on the FINEBoost model fit
#' @export
#'
fineboost_get_coef <- function(X, y, ff){
  res = y
  beta_coef = rep(0, ff$P)
  for(mm in 1:length(ff$csets$filtered_sets)){
    beta_up = as.numeric(coef(lm (res ~ X[,ff$csets$filtered_sets[[mm]][1]]-1))[1])
    beta_coef[ff$csets$filtered_sets[[mm]]] = beta_up/length(ff$csets$filtered_sets[[mm]])
    res = res - X%*%beta_coef
  }
  return(beta_coef)
}
