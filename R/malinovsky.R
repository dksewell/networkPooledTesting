#' Obtain E(C)/E(T) without network information
#' 
#' malinovsky computes the objective function in 
#' Malinovsky et al. (2015) <https://doi.org/10.1111/biom.12385>
#' This gives the expected number of correctly classified
#' individuals divided by the expected number of tests.
#' 
#' @param q numeric between 0 and 1, giving 1 minus the prevalence.
#' @param Sp numeric between 0 and 1, giving the specificity
#' @param Se numeric between 0 and 1, giving the sensitivity
#' @param k integer The number of individuals in each pool
#' @return named list with value, the objective function, 
#' the numerator, the expected number of correctly classified 
#' individuals, and the denominator, the expected number
#' of tests.
malinovsky = function(q, Sp, Se, k){
  numerator = 
    q^k * ( (1-Sp) * (Sp + Se - 1) ) +
    q * ( 1 - Se + Sp*Se - Se^2 ) + Se^2
  denominator = 
    Se - q^k * ( Se + Sp - 1 ) + 1/k
  
  return(list(value = numerator / denominator,
              numerator = numerator,
              denominator = denominator)
  )
}
