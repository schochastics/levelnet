#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int choose(int n, int k)
{
  if (k==0 || k==n)
    return 1;
  return  choose(n-1, k-1) + choose(n-1, k);
}

// [[Rcpp::export]]
NumericMatrix hypergeom(IntegerMatrix P,
                        IntegerVector deg,
                        int m) {
  int n = deg.length();
  NumericMatrix H(n,n);
  IntegerMatrix D(n,n);
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      if(deg[i]<deg[j]){
        D(i,j) = deg[i];
      } else{
        D(i,j) = deg[j];
      }
    }
  }

  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      if(i!=j){
        for(int k=P(i,j); k<=D(i,j);k++){
          // int k = P(i,j);
          // Rcout << i << "," << j << "," << k << "\n";
          // Rcout << deg[i] << "," << deg[j] << "\n";
          H(i,j)+= double(choose(deg[i],k) * choose(m-deg[i],deg[j]-k)) / double(choose(m,deg[j]));
        }
      }
    }
  }
  return H;
}

