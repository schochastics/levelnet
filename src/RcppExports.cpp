// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// getA_cpp
IntegerMatrix getA_cpp(NumericVector x, NumericVector y);
RcppExport SEXP _levelnet_getA_cpp(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(getA_cpp(x, y));
    return rcpp_result_gen;
END_RCPP
}
// mse
IntegerMatrix mse(List adjList, IntegerVector deg);
RcppExport SEXP _levelnet_mse(SEXP adjListSEXP, SEXP degSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type adjList(adjListSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type deg(degSEXP);
    rcpp_result_gen = Rcpp::wrap(mse(adjList, deg));
    return rcpp_result_gen;
END_RCPP
}
// isBipartite
bool isBipartite(IntegerMatrix G, int V, int src);
RcppExport SEXP _levelnet_isBipartite(SEXP GSEXP, SEXP VSEXP, SEXP srcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type G(GSEXP);
    Rcpp::traits::input_parameter< int >::type V(VSEXP);
    Rcpp::traits::input_parameter< int >::type src(srcSEXP);
    rcpp_result_gen = Rcpp::wrap(isBipartite(G, V, src));
    return rcpp_result_gen;
END_RCPP
}
// scobit_loglike_cpp
double scobit_loglike_cpp(NumericVector x1, NumericVector x2, NumericVector y, NumericVector params);
RcppExport SEXP _levelnet_scobit_loglike_cpp(SEXP x1SEXP, SEXP x2SEXP, SEXP ySEXP, SEXP paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(scobit_loglike_cpp(x1, x2, y, params));
    return rcpp_result_gen;
END_RCPP
}
// scobit_loglike_gr_cpp
NumericVector scobit_loglike_gr_cpp(NumericVector x1, NumericVector x2, NumericVector y, NumericVector params);
RcppExport SEXP _levelnet_scobit_loglike_gr_cpp(SEXP x1SEXP, SEXP x2SEXP, SEXP ySEXP, SEXP paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(scobit_loglike_gr_cpp(x1, x2, y, params));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_levelnet_getA_cpp", (DL_FUNC) &_levelnet_getA_cpp, 2},
    {"_levelnet_mse", (DL_FUNC) &_levelnet_mse, 2},
    {"_levelnet_isBipartite", (DL_FUNC) &_levelnet_isBipartite, 3},
    {"_levelnet_scobit_loglike_cpp", (DL_FUNC) &_levelnet_scobit_loglike_cpp, 4},
    {"_levelnet_scobit_loglike_gr_cpp", (DL_FUNC) &_levelnet_scobit_loglike_gr_cpp, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_levelnet(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
