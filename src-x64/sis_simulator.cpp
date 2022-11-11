#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' SIS Simulator
//' 
//' This function simulates data according to the 
//' susceptible-infected-susceptible (SIS) 
//' compartmental model.
//' 
//' @param max_time Integer. Number of days to simulate.
//' @param initial 1xn binary vector indicating initially infected individuals.
//' @param A_mat Sparse adjacency matrix.
//' @param beta positive double. Transmission parameter.
//' @param days_infectious integer.  Number of days an individual remains infectious.
//' @return 
//' @export T x n matrix providing which individuals were infected at each 
//' time point.
// [[Rcpp::export]]
arma::sp_mat sis_simulator(const int & max_time,
                        const arma::rowvec & initial,
                        const arma::sp_mat & A_mat,
                        const double & beta,
                        const int & days_infectious){
  int n = initial.n_cols;
  
  // Set up objects and initialize days of infection 
  arma::sp_mat sis_draws(max_time, n);
  sis_draws.row(0) = initial;
  IntegerVector days_til_susc(n);
  IntegerVector temp1 = seq_len(days_infectious);
  IntegerVector temp2(1);
  arma::sp_mat num_infected_neighbors(1,n);
  double infect_prob = 0.0;
  for(int i=0;i<n;i++){
    if(initial(i) == 1){
      temp2 = RcppArmadillo::sample(temp1,1,FALSE,NumericVector::create());
      days_til_susc(i) = temp2(0);
    }else{
      days_til_susc(i) = days_infectious;
    }
  }
  
  for(int tt = 1;tt < max_time; tt++){
    // Get the number of infected neighbors for each individual
    num_infected_neighbors =
      sis_draws.row(tt - 1) * A_mat;
    
    for(int i=0; i< n; i++){
      if(sis_draws(tt - 1,i) == 0){ 
        // Susceptible to infected
        infect_prob = 1.0 - pow(1.0 - beta,num_infected_neighbors(i));
        if(arma::randu() < infect_prob){
          sis_draws(tt, i) = 1;
        }
      }else{
        // Infected to susceptible
        if(days_til_susc(i) == 1){
          days_til_susc(i) = days_infectious;
          sis_draws(tt, i) = 0;
        }else{
          days_til_susc(i) = days_til_susc(i) - 1;
          sis_draws(tt, i) = 1;
        }
      }
    }
  }
  
  
  return sis_draws;
}

