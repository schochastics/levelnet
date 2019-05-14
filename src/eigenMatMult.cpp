// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>

// [[Rcpp::export]]
SEXP eigenMatMult(Eigen::SparseMatrix<double> A){
  Eigen::SparseMatrix<double> C = A * A.adjoint();

  return Rcpp::wrap(C);
}
