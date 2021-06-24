#ifndef MISC_H
#define MISC_H

#include <RcppArmadillo.h>

double soft_t(double, const double&);

arma::vec soft_t(arma::vec, const arma::vec&);

double log1mexp(double);

arma::vec log1mexp(arma::vec);

float normal_pdf(float);

arma::vec normal_pdf(arma::vec);

#endif
