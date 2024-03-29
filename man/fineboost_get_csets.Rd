% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fineboost_postprocess.R
\name{fineboost_get_csets}
\alias{fineboost_get_csets}
\title{Get confidence sets for local signal clusters from FineBoost output}
\usage{
fineboost_get_csets(ff, X, coverage = 0.95, min_within_LD = 0.5,
  min_between_LD = 0.25, min_clus_centrality = 0.5, nmf_try = 5)
}
\arguments{
\item{ff}{an object of class fineboost, generated from `fineboostR::fineboost_smooth()`}

\item{X}{The attributed X matrix (N times P) obtained using `fineboostR::set_X_attributes()`.}

\item{coverage}{A number between 0 and 1 (close to 1) specifying the coverage of the estimated signal clusters.
Default set to 0.90.}

\item{min_within_LD}{The minimum value of LD permitted for SNPs within a local signal cluster. Default is 0.25.}

\item{min_between_LD}{The minimum value of LD permitted for SNPs across two local signal clusters. Default is it
cannot exceed 0.25.}

\item{nmf_try}{Number of initializations of the Non-negative matrix factorization on the trajectory of weights
tested. Default is set to 10.}

\item{min_cluster_centrality}{The minimum value of cluster centrailty required for a SNP to make the cut in a
local signal cluster. Default is set at 0.5.}
}
\value{
An object with the following items
  \item{unfiltered_sets}{A list of length `Lmax` denoting the local signal clusters before purification.}
  \item{cluster_signal}{The cluster signal for the `unfiltered sets`.}
  \item{filtered_sets}{The local signal clusters that remain after purification procedure.}
  \item{update_time}{A vector of length `Lmax` giving the proportion of boosting updates coming from that local cluster.}
  \item{within_purity}{The maximum LD of SNPs in each individual local cluster.}
  \item{between_purity}{The maximum LD of SNPs between any two distinct local clusters.}
  \item{csets_index}{A vector of length P that labels each variable by the confidence set in the `filtered_sets` it
                     belongs to and 0 otherwise}
}
\description{
Get confidence sets for local signal clusters from FineBoost output
}
