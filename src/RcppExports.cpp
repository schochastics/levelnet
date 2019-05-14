// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// arcDist
double arcDist(NumericVector x, NumericVector y, double r);
RcppExport SEXP _levelnet_arcDist(SEXP xSEXP, SEXP ySEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(arcDist(x, y, r));
    return rcpp_result_gen;
END_RCPP
}
// arcDistMat
NumericMatrix arcDistMat(NumericMatrix X, double r);
RcppExport SEXP _levelnet_arcDistMat(SEXP XSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(arcDistMat(X, r));
    return rcpp_result_gen;
END_RCPP
}
// eigenMatMult
SEXP eigenMatMult(Eigen::SparseMatrix<double> A);
RcppExport SEXP _levelnet_eigenMatMult(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(eigenMatMult(A));
    return rcpp_result_gen;
END_RCPP
}
// choose
int choose(int n, int k);
RcppExport SEXP _levelnet_choose(SEXP nSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(choose(n, k));
    return rcpp_result_gen;
END_RCPP
}
// hypergeom
NumericMatrix hypergeom(IntegerMatrix P, IntegerVector deg, int m);
RcppExport SEXP _levelnet_hypergeom(SEXP PSEXP, SEXP degSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type deg(degSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(hypergeom(P, deg, m));
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
    {"_levelnet_arcDist", (DL_FUNC) &_levelnet_arcDist, 3},
    {"_levelnet_arcDistMat", (DL_FUNC) &_levelnet_arcDistMat, 2},
    {"_levelnet_eigenMatMult", (DL_FUNC) &_levelnet_eigenMatMult, 1},
    {"_levelnet_choose", (DL_FUNC) &_levelnet_choose, 2},
    {"_levelnet_hypergeom", (DL_FUNC) &_levelnet_hypergeom, 3},
    {"_levelnet_isBipartite", (DL_FUNC) &_levelnet_isBipartite, 3},
    {"_levelnet_scobit_loglike_cpp", (DL_FUNC) &_levelnet_scobit_loglike_cpp, 4},
    {"_levelnet_scobit_loglike_gr_cpp", (DL_FUNC) &_levelnet_scobit_loglike_gr_cpp, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_levelnet(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
