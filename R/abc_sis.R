#' ABC sampler for the SIS model
#' 
#' This function runs an approximate Bayesian
#' computation (ABC) algorithm to obtain samples
#' from the SIS model, yielding samples of y
#' given prevalence within a small window of
#' that specified.
#' 
#' @param A_mat nxn adjacency matrix
#' @param prevalence scalar.
#' @param days_infectious integer.
#' @param beta_scalar_range vector with lower and upper bounds for the uniform prior 
#' on the transmission parameter beta
#' @param n_draws How many draws to obtain
#' @param sis_nsims How long should the SIS sampler run (in days)?
#' @param showProgress logical.  Should a progress bar be displayed?
#' @param cl If running in parallel using the parallel package, pass 
#' in object from parallel::makeCluster()
#' @return abc_sis returns a n_draws x (1 + n) matrix.  The first 
#' column is the draws of beta, the remaining columns correspond 
#' to \eqn{y_1,\ldots,y_n}.
#' @import Matrix
#' @import parallel
#' @importFrom rARPACK eigs
abc_sis = function(A_mat,
                   prevalence = 0.025,
                   days_infectious = 7,
                   beta_scalar_range = c(1.25,1.5),
                   n_draws = 5e3,
                   sis_nsims = 300,
                   showProgress = TRUE,
                   cl = NULL){
  
  # Declare objects
  n = nrow(A_mat)
  A_eigs = rARPACK::eigs(A_mat,k=1)
  initial_probs = A_eigs$vectors[,1]/sum(A_eigs$vectors[,1]) * prevalence * nrow(A_mat) # See Newman (2010) p.670
  
  # Draw from uniform prior
  epi_thresh = 1 / A_eigs$values # See, e.g., DOI 10.1007/s00607-011-0155-y
  beta_min = (1 / days_infectious) * epi_thresh
  beta_draws = runif(n_draws,
                     beta_min*beta_scalar_range[1],
                     beta_min*beta_scalar_range[2])
  
  
  # Rejection sampler
  if(is.null(cl)){
    if(showProgress) pb = txtProgressBar(0,n_draws,style=3)
    y_draws = Matrix(0L,n_draws,n)
    for(it in 1:n_draws){
      y_draws[it,] = 
        sis_simulator(max_time = sis_nsims,
                      initial = rbinom(nrow(A_mat),1,initial_probs), 
                      A_mat = A_mat,
                      beta = beta_draws[it],
                      days_infectious = days_infectious)[sis_nsims,]
      
      if(showProgress) setTxtProgressBar(pb,it)
    }
  }else{
    wrapper_fun = function(it){
      temp_data = 
        sis_simulator(max_time = sis_nsims,
                      initial = rbinom(nrow(A_mat),1,initial_probs), 
                      A_mat = A_mat,
                      beta = beta_draws[it],
                      days_infectious = days_infectious)[sis_nsims,]
      return(temp_data)
    }
    
    parallel::clusterEvalQ(cl,{
      library(networkPooledTesting)
    })
    parallel::clusterExport(cl,
                            c("sis_nsims","initial_probs","prevalence","A_mat",
                              "beta_draws","days_infectious","seed_sequence"),
                            envir = environment())
    
    y_draws = parallel::parSapply(cl,1:n_draws,wrapper_fun) %>% t()
  }
  
  beta_y_draws = cbind(beta = beta_draws,
                       y_draws = y_draws)
  return(beta_y_draws)
}

