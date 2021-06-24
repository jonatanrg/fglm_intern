#ifndef LIK_H
#define LIK_H
#include <RcppArmadillo.h>

double lik_ee(double y, const double& yupp, const double& eta, const uint& order);

arma::vec lik_ee(arma::vec y, const arma::vec& eta, const
uint& order);

double obj_fun_ee(arma::vec y, const arma::vec& eta, const arma::vec& b,
           const arma::vec& lam1, const arma::vec& lam2);

Rcpp::List neg_ll_exp_cpp(arma::vec y, arma::mat X, arma::vec b, arma::vec yupp,
                          uint order, double const pen);

Rcpp::List obj_diff_cpp(const arma::vec& y, const arma::mat& X, const arma::vec& b, const
arma::vec& yupp, const arma::vec& lam1, const arma::vec& lam2, const uint& order);

#endif
