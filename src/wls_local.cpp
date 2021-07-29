#include <RcppArmadillo.h> 

using namespace Rcpp;
using namespace arma;


double epa_kernel_wls(double time, double time_zero, double h){
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



double normal_kernel_wls(double time, double time_zero, double h){
  double pre_val = -pow(time-time_zero,2)/2;
  double pre_exp = exp(pre_val);
  double divisor = sqrt(2 * M_PI);
  double post = pre_exp/divisor;
  return post;
}

double triweight_kernel_wls(double time, double time_zero, double h){
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


double quartic_kernel_wls(double time, double time_zero, double h){
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

double kernel_function_wls(double time, double time_zero, double h, int type){
  if (type==1){
    double post = epa_kernel_wls( time, time_zero, h);
    return post;  
  }  else if (type == 2){
    double post = normal_kernel_wls( time, time_zero, h);
    return post;  
  } else if (type == 3){
    double post = triweight_kernel_wls( time, time_zero, h);
    return post;
  } else {
    double post = quartic_kernel_wls( time, time_zero, h);
    return post;
  }
}

arma::vec kernel_vec_wls(arma::vec time, double time_zero, double h, int type){
  int nn = time.size();
  arma::vec post(nn);
  for(int i=0; i < nn; ++i){
    post[i] = kernel_function_wls(time[i], time_zero, h, type);
  }
  return post;
}


//' vcm_wls cpp
//'
//' Weighted Least Squares Function
//' 
//' @param x design matrix
//' @param y response vector
//' @param time vector
//' @param time_zero grid point double
//' @param h bandwidth double
//' @param type kernel function int
//' 
//' @export
//'
//[[Rcpp::export]]

arma::vec vcm_wls_local(arma::mat x, arma::vec y, arma::vec time,
             double time_zero, double h, int type){ // ith observation loglik value
  int nn = y.size();
  arma::vec tt_zero = time - time_zero;
  arma::mat xt = diagmat(tt_zero) * x;
  arma::mat xx = join_rows(x, xt);
  arma::vec pre_ww = kernel_vec_wls(time, time_zero, h, type);
  arma::mat ww = diagmat(pre_ww);
  arma::mat xwx = trans(xx) * ww * xx;
  arma::mat xwy = trans(xx) * ww * y;
  arma::vec post = inv(xwx) * xwy;
  return post;
}


//' tvcm_wls cpp
//'
//' Weighted Least-Squares Function For Longitudinal Data
//' 
//' @param x design matrix
//' @param y response vector
//' @param time vector
//' @param id grouping vector
//' @param time_zero grid point double
//' @param h bandwidth double
//' @param type kernel function int
//' 
//' @export
//'
//[[Rcpp::export]]

arma::mat tvcm_wls_local(arma::mat x, arma::vec y, arma::vec time, arma::vec id,
             double time_zero, double h, int type){ // ith observation loglik value
  int nn = max(id);
  int oo = 2*x.n_cols;
  arma::mat xwx(oo, oo, fill::zeros);
  arma::mat xwy(oo, 1, fill::zeros);
  for (int ii = 0; ii < nn; ++ii){
    arma::uvec indi = find(id==(ii+1));
    arma::mat wk_mat = x.rows(indi);
    arma::vec tt = time.elem(indi);
    arma::vec tt_zero = tt - time_zero;
    arma::mat xt = diagmat(tt_zero) * wk_mat;
    arma::mat xx = join_rows(wk_mat, xt);
    arma::mat yy = y.rows(indi);
    arma::vec pre_ww = kernel_vec_wls(tt, time_zero, h, type);
    arma::mat ww = diagmat(pre_ww);
    arma::mat pxwx = trans(xx) * ww * xx;
    arma::mat pxwy = trans(xx) * ww * yy;
    xwx = xwx + pxwx;
    xwy = xwy + pxwy;
  }
  arma::vec post = inv(xwx) * xwy;
  return post;
}
