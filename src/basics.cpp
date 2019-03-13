#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector  sdRcpp(NumericMatrix x) {
  int n = x.nrow();
  int p = x.ncol();
  NumericVector out (p);
  for(int i = 0; i < p; ++i) {
    NumericVector col_i = x(_,i);
    double sd_i = sd(col_i);
    out[i] = sd_i*sqrt((n-1.0)/n);
  }
  return out;
}

// [[Rcpp::export]]
NumericMatrix scaleRcpp(NumericMatrix x) {
  int nrow = x.nrow(), ncol = x.ncol();
  NumericMatrix out = clone(x);
  for (int j = 0; j < ncol; j++) {
    NumericVector col_j = x(_,j);
    double sd_out=sd(col_j)*sqrt((nrow-1.0)/nrow);
    double mu=mean(col_j);
    for(int i = 0; i < nrow; ++i) {
      out(i,j) = (col_j[i]-mu)/sd_out;
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericMatrix mmultC(NumericMatrix m1, NumericMatrix m2){
  NumericMatrix out(m1.nrow(),m2.ncol());
  NumericVector rm1, cm2;
  for (size_t i = 0; i < m1.nrow(); ++i) {
    rm1 = m1(i,_);
    for (size_t j = 0; j < m2.ncol(); ++j) {
      cm2 = m2(_,j);
      out(i,j) = std::inner_product(rm1.begin(), rm1.end(), cm2.begin(), 0.);
    }
  }
  return out;
}
