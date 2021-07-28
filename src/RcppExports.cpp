// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// bin_cv_loglik
double bin_cv_loglik(arma::mat x, arma::mat beta, arma::vec y, arma::vec time);
RcppExport SEXP _tvcm_bin_cv_loglik(SEXP xSEXP, SEXP betaSEXP, SEXP ySEXP, SEXP timeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type time(timeSEXP);
    rcpp_result_gen = Rcpp::wrap(bin_cv_loglik(x, beta, y, time));
    return rcpp_result_gen;
END_RCPP
}
// bin_logit_loglik
double bin_logit_loglik(arma::mat x, arma::vec beta, arma::vec y, arma::vec time, double time_zero, double h, int type);
RcppExport SEXP _tvcm_bin_logit_loglik(SEXP xSEXP, SEXP betaSEXP, SEXP ySEXP, SEXP timeSEXP, SEXP time_zeroSEXP, SEXP hSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type time(timeSEXP);
    Rcpp::traits::input_parameter< double >::type time_zero(time_zeroSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(bin_logit_loglik(x, beta, y, time, time_zero, h, type));
    return rcpp_result_gen;
END_RCPP
}
// gr_bin_logit_loglik
arma::vec gr_bin_logit_loglik(arma::mat x, arma::vec beta, arma::vec y, arma::vec time, double time_zero, double h, int type);
RcppExport SEXP _tvcm_gr_bin_logit_loglik(SEXP xSEXP, SEXP betaSEXP, SEXP ySEXP, SEXP timeSEXP, SEXP time_zeroSEXP, SEXP hSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type time(timeSEXP);
    Rcpp::traits::input_parameter< double >::type time_zero(time_zeroSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(gr_bin_logit_loglik(x, beta, y, time, time_zero, h, type));
    return rcpp_result_gen;
END_RCPP
}
// hs_bin_logit_loglik
arma::mat hs_bin_logit_loglik(arma::mat x, arma::vec beta, arma::vec y, arma::vec time, double time_zero, double h, int type);
RcppExport SEXP _tvcm_hs_bin_logit_loglik(SEXP xSEXP, SEXP betaSEXP, SEXP ySEXP, SEXP timeSEXP, SEXP time_zeroSEXP, SEXP hSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type time(timeSEXP);
    Rcpp::traits::input_parameter< double >::type time_zero(time_zeroSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(hs_bin_logit_loglik(x, beta, y, time, time_zero, h, type));
    return rcpp_result_gen;
END_RCPP
}
// kernel_function
double kernel_function(double time, double time_zero, double h, int type);
RcppExport SEXP _tvcm_kernel_function(SEXP timeSEXP, SEXP time_zeroSEXP, SEXP hSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type time(timeSEXP);
    Rcpp::traits::input_parameter< double >::type time_zero(time_zeroSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(kernel_function(time, time_zero, h, type));
    return rcpp_result_gen;
END_RCPP
}
// kernel_vec
arma::vec kernel_vec(arma::vec time, double time_zero, double h, int type);
RcppExport SEXP _tvcm_kernel_vec(SEXP timeSEXP, SEXP time_zeroSEXP, SEXP hSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type time(timeSEXP);
    Rcpp::traits::input_parameter< double >::type time_zero(time_zeroSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(kernel_vec(time, time_zero, h, type));
    return rcpp_result_gen;
END_RCPP
}
// pois_cv_loglik
double pois_cv_loglik(arma::mat x, arma::mat beta, arma::vec y, arma::vec time);
RcppExport SEXP _tvcm_pois_cv_loglik(SEXP xSEXP, SEXP betaSEXP, SEXP ySEXP, SEXP timeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type time(timeSEXP);
    rcpp_result_gen = Rcpp::wrap(pois_cv_loglik(x, beta, y, time));
    return rcpp_result_gen;
END_RCPP
}
// pois_log_loglik
double pois_log_loglik(arma::mat x, arma::vec beta, arma::vec y, arma::vec time, double time_zero, double h, int type);
RcppExport SEXP _tvcm_pois_log_loglik(SEXP xSEXP, SEXP betaSEXP, SEXP ySEXP, SEXP timeSEXP, SEXP time_zeroSEXP, SEXP hSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type time(timeSEXP);
    Rcpp::traits::input_parameter< double >::type time_zero(time_zeroSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(pois_log_loglik(x, beta, y, time, time_zero, h, type));
    return rcpp_result_gen;
END_RCPP
}
// gr_pois_log_loglik
arma::vec gr_pois_log_loglik(arma::mat x, arma::vec beta, arma::vec y, arma::vec time, double time_zero, double h, int type);
RcppExport SEXP _tvcm_gr_pois_log_loglik(SEXP xSEXP, SEXP betaSEXP, SEXP ySEXP, SEXP timeSEXP, SEXP time_zeroSEXP, SEXP hSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type time(timeSEXP);
    Rcpp::traits::input_parameter< double >::type time_zero(time_zeroSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(gr_pois_log_loglik(x, beta, y, time, time_zero, h, type));
    return rcpp_result_gen;
END_RCPP
}
// hs_pois_log_loglik
arma::mat hs_pois_log_loglik(arma::mat x, arma::vec beta, arma::vec y, arma::vec time, double time_zero, double h, int type);
RcppExport SEXP _tvcm_hs_pois_log_loglik(SEXP xSEXP, SEXP betaSEXP, SEXP ySEXP, SEXP timeSEXP, SEXP time_zeroSEXP, SEXP hSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type time(timeSEXP);
    Rcpp::traits::input_parameter< double >::type time_zero(time_zeroSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(hs_pois_log_loglik(x, beta, y, time, time_zero, h, type));
    return rcpp_result_gen;
END_RCPP
}
// vcm_wls
arma::vec vcm_wls(arma::mat x, arma::vec y, arma::vec time, double time_zero, double h, int type);
RcppExport SEXP _tvcm_vcm_wls(SEXP xSEXP, SEXP ySEXP, SEXP timeSEXP, SEXP time_zeroSEXP, SEXP hSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type time(timeSEXP);
    Rcpp::traits::input_parameter< double >::type time_zero(time_zeroSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(vcm_wls(x, y, time, time_zero, h, type));
    return rcpp_result_gen;
END_RCPP
}
// tvcm_wls
arma::mat tvcm_wls(arma::mat x, arma::vec y, arma::vec time, arma::vec id, double time_zero, double h, int type);
RcppExport SEXP _tvcm_tvcm_wls(SEXP xSEXP, SEXP ySEXP, SEXP timeSEXP, SEXP idSEXP, SEXP time_zeroSEXP, SEXP hSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type time(timeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type id(idSEXP);
    Rcpp::traits::input_parameter< double >::type time_zero(time_zeroSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(tvcm_wls(x, y, time, id, time_zero, h, type));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_tvcm_bin_cv_loglik", (DL_FUNC) &_tvcm_bin_cv_loglik, 4},
    {"_tvcm_bin_logit_loglik", (DL_FUNC) &_tvcm_bin_logit_loglik, 7},
    {"_tvcm_gr_bin_logit_loglik", (DL_FUNC) &_tvcm_gr_bin_logit_loglik, 7},
    {"_tvcm_hs_bin_logit_loglik", (DL_FUNC) &_tvcm_hs_bin_logit_loglik, 7},
    {"_tvcm_kernel_function", (DL_FUNC) &_tvcm_kernel_function, 4},
    {"_tvcm_kernel_vec", (DL_FUNC) &_tvcm_kernel_vec, 4},
    {"_tvcm_pois_cv_loglik", (DL_FUNC) &_tvcm_pois_cv_loglik, 4},
    {"_tvcm_pois_log_loglik", (DL_FUNC) &_tvcm_pois_log_loglik, 7},
    {"_tvcm_gr_pois_log_loglik", (DL_FUNC) &_tvcm_gr_pois_log_loglik, 7},
    {"_tvcm_hs_pois_log_loglik", (DL_FUNC) &_tvcm_hs_pois_log_loglik, 7},
    {"_tvcm_vcm_wls", (DL_FUNC) &_tvcm_vcm_wls, 6},
    {"_tvcm_tvcm_wls", (DL_FUNC) &_tvcm_tvcm_wls, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_tvcm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
