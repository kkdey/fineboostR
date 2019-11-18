######### CCG : Causal confidence grade is an alternative to Posterior inclusion probability (PIP)  #######

#' @title Compute inclusion confidence grade (CCG) for all variables
#' @param ff an object of class fineboos, generated from `fineboostR::fineboost_smooth()`
#' @param use_csets Whether to use the confidence sets when computing the ICG. Defaults to TRUE.
#' @return a vector inclusion confidence grade for all variables.
#' @export
#'
fineboost_get_ccg <- function (ff, use_csets = TRUE){
  if (class(ff) == "fineboost"){
      ccg = rep(0, ff$P)
      include_idx = unlist(ff$csets$unfiltered_sets)
      signal = ff$csets$cluster_signal
      refined_signal = matrix(0, nrow(signal), ncol(signal))
      for(kk in 1:nrow(signal)){
        sorted_signal = sort(signal[kk,], decreasing = T)
        pp = 1
        counter = 0
        while (pp <=5 & min(sorted_signal[1:pp]) > 1/(pp+1)){
          counter = pp
          pp = pp + 1
        }
        if(counter > 0){
          refined_signal[kk, order(signal[kk,], decreasing = T)[1:counter]] = 1/counter
        }else{
          refined_signal[kk,] = signal[kk,]
        }
      }
      ccg = as.vector(1 - apply(1 - refined_signal, 2, prod))
      return(ccg)
  }
}


######### Confidence sets generated from the output of Fineboost  #######

#' @title Get confidence sets for local signal clusters from FineBoost output
#' @param ff an object of class fineboost, generated from `fineboostR::fineboost_smooth()`
#' @param X The attributed X matrix (N times P) obtained using `fineboostR::set_X_attributes()`.
#' @param coverage A number between 0 and 1 (close to 1) specifying the coverage of the estimated signal clusters.
#'                 Default set to 0.90.
#' @param min_within_LD The minimum value of LD permitted for SNPs within a local signal cluster. Default is 0.25.
#' @param min_between_LD The minimum value of LD permitted for SNPs across two local signal clusters. Default is it
#'                       cannot exceed 0.25.
#' @param min_cluster_centrality The minimum value of cluster centrailty required for a SNP to make the cut in a
#'                      local signal cluster. Default is set at 0.5.
#' @param nmf_try Number of initializations of the Non-negative matrix factorization on the trajectory of weights
#'                tested. Default is set to 10.
#'
#' @return An object with the following items
#'   \item{unfiltered_sets}{A list of length `Lmax` denoting the local signal clusters before purification.}
#'   \item{cluster_signal}{The cluster signal for the `unfiltered sets`.}
#'   \item{filtered_sets}{The local signal clusters that remain after purification procedure.}
#'   \item{update_time}{A vector of length `Lmax` giving the proportion of boosting updates coming from that local cluster.}
#'   \item{within_purity}{The maximum LD of SNPs in each individual local cluster.}
#'   \item{between_purity}{The maximum LD of SNPs between any two distinct local clusters.}
#'   \item{csets_index}{A vector of length P that labels each variable by the confidence set in the `filtered_sets` it
#'                      belongs to and 0 otherwise}
#' @export
#'
fineboost_get_csets <- function (ff,
                                 X,
                                 coverage = 0.95,
                                 min_within_LD = 0.5,
                                 min_between_LD = 0.25,
                                 min_clus_centrality = 0.5,
                                 nmf_try = 5){
  if (class(ff) == "fineboost"){
    vx = vector(mode="list", nmf_try)
    ff$weights_path = ff$weights_path[1:floor(0.95*nrow(ff$weights_path)),]
    ff$weights_path[which(ff$weights_path < 0.01)] = 0
    for(rr in 1:nmf_try){
      Lmax = ff$Lmax
      set.seed(rr)
      res = NNLM::nnmf(ff$weights_path, k=Lmax, method ="scd", loss = "mse")
      W2 = t(apply(res$W, 1, function(x) return(x/sum(x))))
      clus_med = c()
      presence_prop = c()
      for(k in 1:Lmax){
        idx = which(W2[,k] > 0.5)
        if(length(idx) > ceiling(0.01*nrow(ff$weights_path))){
          cc = apply(ff$weights_path[idx, , drop=F], 2, median)
          clus_med = as.matrix(rbind(clus_med, cc/sum(cc)))
          presence_prop = c(presence_prop, length(idx)/nrow(W2))
        }
      }

      confidence_sets = list()
      for(mm in 1:nrow(clus_med)){
        tempp = order(clus_med[mm, ], decreasing=T)
        confidence_sets[[mm]] = tempp[1:min(which(cumsum(clus_med[mm,tempp]) > coverage))]
      }

      evidence_strength = presence_prop
      confidence_sets = confidence_sets[order(evidence_strength,
                                                decreasing = T)]
      clus_med = clus_med[order(evidence_strength,
                                decreasing = T),, drop=F]
      # within_LD = c()
      # between_LD = c()
      # for(ee in 1:length(evidence_strength)){
      #   within_LD = min(abs(attr(X, "LD"))[confidence_sets[[ee]], confidence_sets[[ee]]])
      #   if(temp2 < 0.5){temp2 = 0}
      #   evidence_strength[ee] = evidence_strength[ee]*temp2
      # }

      within_purity = c()
      for(ee in 1:length(evidence_strength)){
        within_purity = c(within_purity, quantile((abs(attr(X, "LD"))[confidence_sets[[ee]],
                                confidence_sets[[ee]]])%*%clus_med[ee, confidence_sets[[ee]]], 0.90))
      }


      cluster_centrality = c()
      confidence_sets_1 = confidence_sets
      for(ee in 1:length(evidence_strength)){
        cluster_centrality = (abs(attr(X, "LD"))[confidence_sets[[ee]],
                            confidence_sets[[ee]]])%*%clus_med[ee, confidence_sets[[ee]]]
        confidence_sets_1[[ee]] = confidence_sets[[ee]][which(cluster_centrality > min_clus_centrality)]
      }

      between_purity = matrix(0, length(evidence_strength), length(evidence_strength))
      between_purity2 = matrix(0, length(evidence_strength), length(evidence_strength))


      if(length(evidence_strength) > 1){
        for(ee in 2:length(evidence_strength)){
          for(gg in 1:(ee-1)){
            between_purity[ee,gg] = min(abs(attr(X, "LD"))[confidence_sets[[ee]], confidence_sets[[gg]]])
          }
        }

        for(ee in 2:length(evidence_strength)){
          for(gg in 1:(ee-1)){
            between_purity2[ee,gg] = max(abs(attr(X, "LD"))[confidence_sets[[ee]], confidence_sets[[gg]]])
          }
        }
        tmp2 = apply(between_purity, 1, max)
        tmp3 = apply(between_purity2, 1, max)
        purity_score = rep(1, length(evidence_strength))
        purity_score[which(tmp2 > min_between_LD | tmp3 > 0.90)] = 0
        purity_score[which(within_purity < min_within_LD)] = 0
      }



      confidence_sets_1 = confidence_sets_1[which(purity_score == 1)]
      evidence_strength = sort(evidence_strength, decreasing = T)

      csets_index = rep(0, ff$P)
      for(kk in 1:length(confidence_sets_1)){
        csets_index[confidence_sets_1[[kk]]] = kk
      }

      ll = list("unfiltered_sets" = confidence_sets,
                "cluster_signal" = clus_med,
                "filtered_sets" = confidence_sets_1,
                "update_time" = evidence_strength,
                "within_purity" = within_purity,
                "between_purity" = between_purity,
                "csets_index" = csets_index)
      vx[[rr]] = ll
    }

    ############  Find NMF run that has highest aggr. correlation across multiple runs ##################
    temp = matrix(0, nmf_try, ff$P)
    for(rr in 1:nmf_try){
      temp[rr, which(vx[[rr]]$csets_index != 0)] = 1
    }
    best_run = which.max(rowSums(cor(t(temp))))[1]
    ll2 = vx[[best_run]]

    for(mm in 1:length(ll2$filtered_sets)){
      ll2$filtered_sets[[mm]] = ll2$filtered_sets[[mm]][order(ll2$cluster_signal[mm, ll2$filtered_sets[[mm]]],
                                                              decreasing = T)]
    }

    return(ll2)
  }
}

