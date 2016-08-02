#include <Rcpp.h>

using namespace Rcpp;

// Learn more about how to use Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//
// and browse examples of code using Rcpp at:
//
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector aMatrix(NumericMatrix A, IntegerVector s, IntegerVector d) {
  int Nrow = A.nrow(); int nrow = Nrow -1;

  for (int i = 0; i < nrow; i++) {
    A(i,i) = 1 + (A(s(i), d(i)))/2;

    for (int j = i+1; j < Nrow; j++) {
      if (j >= nrow)
        break;
      A(i,j) = (A(i, s(j)) + A(i,d(j)))/2;
      A(j,i) = A(i,j);
    }
  }

  NumericVector inbF(nrow);

  for (int a = 0; a < nrow; a++){
    inbF(a) = A(a,a)-1;
  }

  return(inbF);
}

// [[Rcpp::export]]
IntegerVector calcGen(IntegerVector s, IntegerVector d, int nrow, IntegerVector genVec) {
  for (int i = 0; i < nrow; i++) {
    if (s(i) != 0 & d(i) != 0) {
      int x = genVec(s(i)-1); int y = genVec(d(i)-1);
      if (x < 0){x=0;}; if (y<0){y=0;};
      genVec(i) = std::max(x, y) + 1;

    }
  }
  return(genVec);
}
