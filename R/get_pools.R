#' Obtain optimal pooling assignments
#' 
#' This functions takes as an input a network and a matrix of
#' data generated from some transmission model (e.g., SIS model),
#' and outputs the optimal pooling assignment as determined 
#' through a simulated annealing algorithm.
#' 
#' @param A igraph object.  Network conveying transmission opportunities between
#' individuals
#' @param sensitivity numeric.
#' @param specificity numeric.
#' @param K integer. Number of individuals in each pool
#' @param y_draws matrix. (# draws) x (# individuals) matrix.  Sparse format
#' (dgCMatrix) is highly recommended for computational efficiency.
#' @param Z_initial (optional) matrix giving pooling assignments to initialize
#' the simulated annealing algorithm. If left blank (recommended), a constrained 
#' PAM clustering algorithm will be used.
#' @param iter_per_temp integer. Number of iterations per SA temperature
#' @param maxIterConstant integer.  Number of iterations with no change before 
#' stopping the SA algorithm.
#' @param showProgress logical.  Should a progress bar be displayed?
#' @param temperatures positive monotonically decreasing numeric vector.  SA 
#' cooling schedule
#' @param extra_individuals character describing what to do with the remainder of
#' n/P. Either "combine" to evenly distribute the remainder across existing pools, 
#' or "new_cluster" to create a new pool of less than K individuals.
#' 
#' @return get_pools returns a list with named elements Z, the nxP matrix giving
#' the optimal pool assignments, and fvalues, the vector of the objective function
#' over the SA iterations.
#' @import igraph
#' @import Matrix
#' @importFrom intergraph asIgraph
#' @importFrom dplyr near
#' @importFrom graphics par
#' @importFrom methods as
#' @importFrom stats median rbinom runif
#' @importFrom utils txtProgressBar setTxtProgressBar
get_pools = function(A,
                     sensitivity = 0.75,
                     specificity = 0.995,
                     K = 10,
                     y_draws,
                     Z_initial,
                     iter_per_temp = 100,
                     maxIterConstant = 1e3,
                     showProgress = TRUE, 
                     temperatures = 2*0.95^c(1:500-1),
                     extra_individuals = c("combine","new_cluster")[1]){
  
  # Get A as an igraph object, and A_mat as an adjacency matrix
  if(class(A) != "igraph"){
    if(class(A) == "network"){ A <- intergraph::asIgraph(A); A_mat = as_adjacency_matrix(A)}
    if( (class(as.matrix(A)) %in% "matrix") ){A_mat <- A; A <- graph_from_adjacency_matrix(A)}
  }else{
    A_mat = as_adjacency_matrix(A)
  }
  
  # Get basic objects
  n = vcount(A)
  n_temps = length(temperatures)
  nsim = nrow(y_draws)
  Dists = igraph::distances(A)
  DD = 1/Dists
  diag(DD) = 0
  
  
  # Make sure the y draws are in a sparse matrix
  if(class(y_draws) != "dgCMatrix") y_draws = as(y_draws,"dgCMatrix")
  
  # Either combine extra individuals into larger pools or set remainder in a new pool
  if(extra_individuals == "combine"){
    P = floor(n / K)
    Z = Diagonal(n = P) %x% matrix(1,K,1)
    row_difference = n - nrow(Z)
    while(row_difference > 0){
      Z = rbind(Z,
                Diagonal(n = P)[1:min(P,row_difference),])
      row_difference = n - nrow(Z)
    }
    rm(row_difference)
  }
  if(extra_individuals == "new_cluster"){
    P = ceiling(n / K)
    Z = Diagonal(n = P) %x% matrix(1,K,1)
    Z = Z[1:n,]
  }
  K_p = colSums(Z)
  
  # Initialize Z
  if(missing(Z_initial)){
    inits = initialize(Dists,colSums(Z))
    Z[] = 0
    Z[cbind(1:n,inits)] = 1L
  }else{
    Z = Z_initial
  }
  Z = as(Z,"dgCMatrix")
  Z_best = Z
  
  
  # Initialize the pool-to-pool matrix for the proposal 
  SS = matrix(0,P,P)
  clust_assign = lapply(1:P,function(k)which(dplyr::near(Z[,k],1)))
  for(ell in 1:(P-1)){
    for(m in (ell+1):P){
      SS[ell,m] = mean(DD[clust_assign[[ell]],clust_assign[[m]]])
    }
  }
  SS_index = which(upper.tri(SS),arr.ind=T)
  
  # Objective function
  objective_fun = function(Z){
    yZ = y_draws %*% Z
    mu_p = colMeans(yZ)
    Pr0 = (nrow(yZ) - diff(yZ@p)) / nrow(yZ)
    
    num_function = function(k){
      K_p[k] * sensitivity^2 +
        (K_p[k] - mu_p[k]) * (sensitivity * specificity + 1 - sensitivity - sensitivity^2) +
        K_p[k] * Pr0[k] * (1 - specificity) * (specificity + sensitivity - 1)
    }
    denom_function = function(k){
      1 +
        K_p[k] * sensitivity +
        K_p[k] * Pr0[k] * (1 - specificity - sensitivity)
    }
    
    log( sum(sapply(1:P,num_function)) ) - 
      log( sum(sapply(1:P,denom_function)) )
  }
  
  # Set algorithm values
  constant_counter = 0
  value_tracker = numeric(n_temps*iter_per_temp + 1) 
  value_tracker[1] = 
    objective_fun(Z_best)
  best_value = value_tracker[1]
  
  
  # Begin Simulated Annealing
  if(showProgress) pb = txtProgressBar(0,n_temps*iter_per_temp,style=3)
  for(temper in 1:n_temps){
    if(constant_counter >= maxIterConstant) break
    
    for(it in 1:iter_per_temp){
      # Proposal
      ell_m = sample(P*(P-1)/2,1,prob = SS[upper.tri(SS)])
      p1 = SS_index[ell_m,1]
      p2 = SS_index[ell_m,2]
      i1 = sample(clust_assign[[p1]],1)
      i2 = sample(clust_assign[[p2]],1)
      Z_tilde = Z
      Z_tilde[c(i1,i2),] = Z[c(i2,i1),]
      
      new_obj_value = objective_fun(Z_tilde)
      value_difference = 
        new_obj_value - value_tracker[iter_per_temp*(temper - 1) + it]
      
      accProb = exp(value_difference/temperatures[temper])
      
      if(runif(1) <= accProb){
        # if(accProb > 1){
        
        clust_assign[[SS_index[ell_m,1]]] = 
          c(setdiff(clust_assign[[SS_index[ell_m,1]]],i1),i2)
        clust_assign[[SS_index[ell_m,2]]] = 
          c(setdiff(clust_assign[[SS_index[ell_m,2]]],i2),i1)
        
        for(j in setdiff(1:P,SS_index[ell_m,1])){
          SS[min(j,SS_index[ell_m,1]),max(j,SS_index[ell_m,1])] = 
            mean(DD[clust_assign[[SS_index[ell_m,1]]],
                    clust_assign[[j]] ])
        }
        for(j in setdiff(1:P,SS_index[ell_m,])){
          SS[min(j,SS_index[ell_m,2]),max(j,SS_index[ell_m,2])] = 
            mean(DD[clust_assign[[SS_index[ell_m,2]]],
                    clust_assign[[j]] ])
        }
        
        Z = Z_tilde
        
        value_tracker[1 + iter_per_temp*(temper - 1) + it] = 
          new_obj_value
        
        if(value_tracker[1 + iter_per_temp*(temper - 1) + it] > best_value ){
          best_value = value_tracker[1 + iter_per_temp*(temper - 1) + it]
          Z_best = Z
        }
        constant_counter = 0
        
      }else{
        
        value_tracker[1 + iter_per_temp*(temper - 1) + it] = 
          value_tracker[iter_per_temp*(temper - 1) + it]
        constant_counter = constant_counter + 1
        
      }
      
      if(constant_counter >= maxIterConstant) break
      
      if(showProgress) setTxtProgressBar(pb,iter_per_temp*(temper - 1) + it)
    }
    if(showProgress) plot(value_tracker[1:(1 + iter_per_temp*(temper - 1) + it)],type='l')
    
  }
  
  
  par(mar=c(4,4,1,1))
  plot(value_tracker[1:(1 + iter_per_temp*(temper - 1 - (constant_counter == maxIterConstant)) + it)],
       type='l',
       ylab = "fvalues",
       xlab = "Iteration")
  
  return(list(Z = Z_best,
              fvalues = value_tracker[1:(1 + iter_per_temp*(temper - 1 - (constant_counter == maxIterConstant)) + it)]))
}
