#include "misc.h"
#include "likelihood.h"
#include <cmath>
#include <limits>
#include <RcppArmadillo.h>

double lik_ee(double y, const double& eta,const uint& order)
{
  double cdf = 0.5*(1+erfc(-eta));
  double pdf = normal_pdf(eta);
  double out = 0;
  if (order == 0){
      out = y*std::log(cdf) + (1-y)*std::log(1-cdf);
    }
  else if(order == 1){
      out = -y*pdf/cdf + (1-y)*pdf/(1-cdf);
    }
  return out;
}

arma::vec lik_ee(arma::vec y, const arma::vec& eta,const uint& order)
{
  for(uint ii = 0; ii < y.n_elem; ii++){
    y(ii) = lik_ee(y(ii), eta(ii),order);
  }
  return y;
}

double obj_fun_ee(arma::vec y, const arma::vec& eta,const uint& order,
const arma::vec& b, const arma::vec& lam1, const arma::vec& lam2)
{
  double obj = -arma::mean(lik_ee(y, eta,0));
  obj += arma::sum(lam1 % arma::abs(b)) + 0.5 * arma::sum(lam2 %
  arma::square(b));
  return obj;
}

// [[Rcpp::export]]
Rcpp::List obj_diff_cpp(const arma::vec& y, const arma::mat& X, const arma::vec& b, const uint& order, const arma::vec& lam1, const arma::vec& lam2)
{
  const uint p = X.n_cols;
  arma::vec eta = X * b;
  double obj = obj_fun_ee(y, eta, order, b, lam1, lam2);
  arma::vec sub_grad(p, arma::fill::zeros);
  //Rcpp::Rcout << -arma::mean(X.each_col() % lik_ee(y, yupp, eta, 1), 0) <<std::endl;
  sub_grad = -arma::mean(X.each_col() % lik_ee(y, eta, order), 0).t();
  sub_grad += lam2 % b + lam1 % arma::sign(b);

  return Rcpp::List::create(Rcpp::Named("obj") = obj, Rcpp::Named("sub_grad") =
  sub_grad);
}
