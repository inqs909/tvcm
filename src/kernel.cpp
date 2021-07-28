#include <RcppArmadillo.h> 

using namespace Rcpp;
using namespace arma;


double epa_kernel(double time, double time_zero, double h){
  double pre_val = sqrt(pow(time - time_zero,2));
  double val = pre_val / h;
  double pre_epa = 3 * (1 - pow(val, 2)) / 4;
  double epa = pre_epa /h;
  if (val <= 1){
    if (epa > 0){
      double post = epa;
      return post;
    } else {
      double post = 0;
      return post;
    }
  } else {
    double post = 0;
    return post;
  }
}


double normal_kernel(double time, double time_zero, double h){
  double pre_val = -pow(time-time_zero,2)/2;
  double pre_exp = exp(pre_val);
  double divisor = sqrt(2 * M_PI);
  double post = pre_exp/divisor;
  return post;
}


double triweight_kernel(double time, double time_zero, double h){
  double pre_val = sqrt(pow(time - time_zero,2));
  double val = pre_val / h;
  double pre_triw = 35 * pow(1 - pow(val, 2),3) / 32;
  double triw = pre_triw /h;
  if (val <= 1){
    if (triw > 0){
      double post = triw;
      return post;
    } else {
      double post = 0;
      return post;
    }
  } else {
    double post = 0;
    return post;
  }
}


double quartic_kernel(double time, double time_zero, double h){
  double pre_val = sqrt(pow(time - time_zero,2));
  double val = pre_val / h;
  double pre_q = 15 * pow(1 - pow(val, 2),2) / 16;
  double q = pre_q /h;
  if (val <= 1){
    if (q > 0){
      double post = q;
      return post;
    } else {
      double post = 0;
      return post;
    }
  } else {
    double post = 0;
    return post;
  }
}


//' kernel_function
//'
//' Computes the weight using a kernel function
//' 
//' @param time time point
//' @param time_zero grid point
//' @param h bandwidth
//' @param type kernel function
//' 
//' @export
//'
// [[Rcpp::export]]
double kernel_function(double time, double time_zero, double h, int type){
  if (type==1){
    double post = epa_kernel( time, time_zero, h);
    return post;  
  }  else if (type == 2){
    double post = normal_kernel( time, time_zero, h);
    return post;  
  } else if (type == 3){
    double post = triweight_kernel( time, time_zero, h);
    return post;
  } else {
    double post = quartic_kernel( time, time_zero, h);
    return post;
  }
}


//' kernel_vec
//'
//' Vectorized kernel function
//' 
//' @param time time point as a vector
//' @param time_zero grid point
//' @param h bandwidth
//' @param type kernel function
//' 
//' @export
//'
// [[Rcpp::export]]
arma::vec kernel_vec(arma::vec time, double time_zero, double h, int type){
  int nn = time.size();
  arma::vec post(nn);
  for(int i=0; i < nn; ++i){
    post[i] = kernel_function(time[i], time_zero, h, type);
  }
  return post;
}

