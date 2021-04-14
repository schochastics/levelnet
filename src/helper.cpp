#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix getA_cpp(NumericVector x,NumericVector y){
  int n = x.length();
  IntegerMatrix A(n,n);
  for(int i=0;i<n;++i){
    for(int j=0;j<n;++j){
      if( ((x[i] >= x[j]) & (x[i]<=y[j])) | ((x[i] <= x[j]) & (y[i]>=x[j])) ){
        A(i,j) = 1;
      }
    }
  }
  return A;
}

IntegerMatrix getxy_cpp(List N,IntegerVector perm){
  int n = perm.length();
  IntegerMatrix xy(n,2);
  for(int i=0;i<n;++i){
    IntegerVector v1 = N[i];
    int m = v1.length();
    for(int j=0;j<m;++j){
      v1[j] = perm[v1[j]-1];
    }
    xy(i,0) = min(v1);
    xy(i,1) = perm[i];
  }
  return xy;
}

// [[Rcpp::export]]
IntegerMatrix mse(List adjList, IntegerVector deg) {
  int n=deg.size();
  IntegerVector marked(n);
  IntegerVector t(n);
  IntegerMatrix dom(n,n);
  for(int v = 0; v < n; ++v) {
    Rcpp::checkUserInterrupt();
    // IntegerVector Nv = as<IntegerVector>(adjList[v]);

    std::vector<int> Nv = as<std::vector<int> >(adjList[v]);
    //check if isolate
    if (Nv.empty()) {
      for(int j = 0; j < n; ++j){
        dom(v,j)=1;
      }
      dom(v,v)=0;
    }

    for(std::vector<int>::size_type j = 0; j!=Nv.size(); ++j){
      // for(int j = 0; j < NvSize; ++j) {
      int u = Nv[j];
      std::vector<int> Nu=adjList[u];
      Nu.push_back(u);
      for(std::vector<int>::size_type i = 0; i!=Nu.size(); ++i){
        // for(int i = 0; i<Nu.size(); ++i){
        int w=Nu[i];
        if(w!=v){
          if(marked[w]!=v){
            marked[w]=v;
            t[w]=0;
          }
          t[w]+=1;
          if(t[w]==deg[v]){
            dom(v,w)=1;
          }
        }
      }
    }
  }
  return dom;
}
