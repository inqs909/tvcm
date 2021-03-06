#include <RcppArmadillo.h> 

using namespace Rcpp;
using namespace arma;


double epa_kernel_pois(double time, double time_zero, double h){
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



double normal_kernel_pois(double time, double time_zero, double h){
  double pre_val = -pow(time-time_zero,2)/2;
  double pre_exp = exp(pre_val);
  double divisor = sqrt(2 * M_PI);
  double post = pre_exp/divisor;
  return post;
}

double triweight_kernel_pois(double time, double time_zero, double h){
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


double quartic_kernel_pois(double time, double time_zero, double h){
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

double kernel_function_pois(double time, double time_zero, double h, int type){
  if (type==1){
    double post = epa_kernel_pois( time, time_zero, h);
    return post;  
  }  else if (type == 2){
    double post = normal_kernel_pois( time, time_zero, h);
    return post;  
  } else if (type == 3){
    double post = triweight_kernel_pois( time, time_zero, h);
    return post;
  } else {
    double post = quartic_kernel_pois( time, time_zero, h);
    return post;
  }
}


double fact (int x){
  if (x == 0){
    return 1;
  } else if (x == 1){
    return 1;
  } else {
    int ide = x - 1;
    arma::vec ff = linspace(1, x, x);
    arma::vec ccpp = cumprod(ff);
    double post = ccpp[ide];
    return post;
  }
}

//' pois_cv_loglik cpp
//'
//' LogLikelihood function for cv
//' 
//' @param x design matrix
//' @param beta vcm vector
//' @param y response vector
//' @param time vector
//' 
//' @export
//'
// [[Rcpp::export]]

double pois_cv_loglik(arma::mat x, arma::mat beta, arma::vec y, arma::vec time){ // ith observation loglik value
  int nn = time.size();
  arma::vec res(nn);
  for (int ii=0; ii < nn; ++ii){
    arma::rowvec bb = beta.row(ii);
    arma::rowvec pre_x = x.row(ii);
    arma::rowvec zz(bb.size()/2, fill::zeros);
    arma::rowvec xx = join_rows(pre_x, zz);
    double yy = y[ii];
    double tt = time[ii];
    double linear = arma::as_scalar(xx * trans(bb));
    res[ii] = yy * linear - exp(linear) - log(fact(yy));
  }
  double post = - sum(res);
  return post;
}

//' pois_log_loglik cpp
//'
//' LogLikilihood function for vcm
//' 
//' @param x design matrix
//' @param beta vcm vector
//' @param y response vector
//' @param time vector
//' @param time_zero grid point double
//' @param h bandwidth double
//' @param type kernel function int
//' 
//' @export
//'
//[[Rcpp::export]]

double pois_log_loglik(arma::mat x, arma::vec beta, arma::vec y,
                       arma::vec time, double time_zero,
                        double h, int type){ // ith observation loglik value
  int nn = y.size();
  arma::vec res(nn);
  for (int ii=0; ii < nn; ++ii){
    arma::rowvec pre_x = x.row(ii);
    double yy = arma::as_scalar(y.row(ii));
    double ff = fact(yy);
    double tt = arma::as_scalar(time.row(ii));
    double tt_zero = tt - time_zero;
    arma::rowvec pre_xt = pre_x * tt_zero;
    arma::rowvec xx = join_rows(pre_x, pre_xt);
    double linear = arma::as_scalar(xx * beta);
    res[ii] = (yy * linear - exp(linear) - log(ff)) * kernel_function_pois(tt, time_zero, h, type);
  }
  double post = - mean(res);
  return post;
}

//' gr_pois_log_loglik cpp
//'
//' Gradient function for vcm
//' 
//' @param x design matrix
//' @param beta vcm vector
//' @param y response vector
//' @param time vector
//' @param time_zero grid point double
//' @param h bandwidth double
//' @param type kernel function int
//' 
//' @export
//' 
//[[Rcpp::export]]

arma::vec gr_pois_log_loglik(arma::mat x, arma::vec beta, arma::vec y,
                       arma::vec time, double time_zero,
                        double h, int type){ // ith observation loglik value
  int nn = y.size();
  int oo = beta.size();
  arma::mat out = mat(oo, nn);
  
  for (int ii=0; ii < nn; ++ii){
    arma::rowvec pre_x = x.row(ii);
    double yy = arma::as_scalar(y.row(ii));
    double tt = arma::as_scalar(time.row(ii));
    double tt_zero = tt - time_zero;
    arma::rowvec pre_xt = pre_x * tt_zero;
    arma::rowvec xx = join_rows(pre_x, pre_xt);
    double linear = arma::as_scalar(xx * beta);
    double res = (yy - exp(linear)) * kernel_function_pois(tt, time_zero, h, type);
    out.col(ii) =  trans(xx) * res;
  }
  
  arma::vec pre = vectorise(sum(out,1));
  arma::vec post =  - pre / nn;
  
  return post;
}

//' hs_pois_log_loglik cpp
//'
//' Hessian function for vcm
//' 
//' @param x design matrix
//' @param beta vcm vector
//' @param y response vector
//' @param time vector
//' @param time_zero grid point double
//' @param h bandwidth double
//' @param type kernel function int
//' 
//' @export
//' 
//[[Rcpp::export]]

arma::mat hs_pois_log_loglik(arma::mat x, arma::vec beta, arma::vec y,
                       arma::vec time, double time_zero,
                        double h, int type){ // ith observation loglik value
  int nn = y.size();
  int oo = beta.size();
  arma::mat hess(oo, oo, fill::zeros);
  for (int ii=0; ii < nn; ++ii){
    arma::rowvec pre_x = x.row(ii);
    double tt = arma::as_scalar(time.row(ii));
    double tt_zero = tt - time_zero;
    arma::rowvec pre_xt = pre_x * tt_zero;
    arma::rowvec xx = join_rows(pre_x, pre_xt);
    double linear = arma::as_scalar(xx * beta);
    arma::mat xx_mat = trans(xx) * xx;
    double res = - exp(linear) * kernel_function_pois(tt, time_zero, h, type);
    hess = hess + xx_mat *  res ;
  }
  
  arma::mat post = - hess / nn;
  return post;
}

