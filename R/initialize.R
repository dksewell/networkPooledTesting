#' @importFrom cluster pam
initialize = function(Dists,K_vec){
  P = length(K_vec)
  inf_index = which(Dists == Inf, arr.ind=T)
  if(length(inf_index) > 0) Dists[inf_index] = nrow(Dists) - 1
  km = cluster::pam(Dists,P,diss=TRUE)
  
  dist_to_clusters = matrix(0,nrow(Dists),P)
  median_dists = numeric( nrow(Dists))
  for(i in 1: nrow(Dists)){
    dist_to_clusters[i,] = Dists[i,km$medoids]
    lowest_dist = min(dist_to_clusters[i,])
    median_dists[i] = median((dist_to_clusters[i,] - lowest_dist)^2)
  }
  
  Ord = order(median_dists, decreasing = TRUE)
  cluster_assignments = integer(nrow(Dists))
  cluster_counts = K_vec
  dist_to_clusters1 = dist_to_clusters
  for(i in Ord){
    best_cluster = which.min(dist_to_clusters1[i,])
    cluster_assignments[i] = best_cluster
    cluster_counts[best_cluster] = cluster_counts[best_cluster] - 1
    if(cluster_counts[best_cluster] == 0) {
      dist_to_clusters1[,best_cluster] = Inf
    }
  }
  
  return(cluster_assignments)
}
