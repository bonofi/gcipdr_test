#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double mccovx1x2_cpp(double rz, Rcpp::Function Gx1, Rcpp::Function Gx2,
                     double rx, arma::vec meanx, arma::vec sdx,
                     bool pNorm_1, bool pNorm_2, int K)
{
  // Generate bivariate normal samples in C++ (Armadillo)
  arma::mat Sigmaz = arma::eye(2, 2);
  Sigmaz(0,1) = Sigmaz(1,0) = rz;
  arma::mat z = arma::mvnrnd(arma::zeros(2), Sigmaz, K);
  
  // Apply pnorm in C++ (Armadillo normcdf)
  arma::vec p1 = arma::normcdf(z.row(0).t());
  arma::vec p2 = arma::normcdf(z.row(1).t());
  
  // If pNorm flags are false, use z directly instead of pnorm(z)
  if (!pNorm_1) p1 = z.row(0).t();
  if (!pNorm_2) p2 = z.row(1).t();
  
  // Call R marginal inverse CDFs (vectorized)
  Rcpp::NumericVector p1_r = Rcpp::wrap(p1);
  Rcpp::NumericVector p2_r = Rcpp::wrap(p2);
  Rcpp::NumericVector x1_r = Rcpp::as<Rcpp::NumericVector>(Gx1(p1_r));
  Rcpp::NumericVector x2_r = Rcpp::as<Rcpp::NumericVector>(Gx2(p2_r));
  arma::vec x1 = Rcpp::as<arma::vec>(x1_r);
  arma::vec x2 = Rcpp::as<arma::vec>(x2_r);
  
  // Compute covariance in C++
  double covx1x2 = arma::sum(x1 % x2) / K;
  double res = (covx1x2 - arma::prod(meanx)) / arma::prod(sdx);
  
  return res - rx;
}

// [[Rcpp::export]]
double mccovx1x2prime_cpp(double rz, Rcpp::Function Gx1, Rcpp::Function Gx2,
                          arma::vec sdx, bool pNorm_1, bool pNorm_2, int K)
{
  // Generate bivariate normal samples in C++
  arma::mat Sigmaz = arma::eye(2, 2);
  Sigmaz(0,1) = Sigmaz(1,0) = rz;
  arma::mat z = arma::mvnrnd(arma::zeros(2), Sigmaz, K);
  
  // Apply pnorm in C++
  arma::vec p1 = arma::normcdf(z.row(0).t());
  arma::vec p2 = arma::normcdf(z.row(1).t());
  
  if (!pNorm_1) p1 = z.row(0).t();
  if (!pNorm_2) p2 = z.row(1).t();
  
  // Call R marginal inverse CDFs
  Rcpp::NumericVector p1_r = Rcpp::wrap(p1);
  Rcpp::NumericVector p2_r = Rcpp::wrap(p2);
  Rcpp::NumericVector x1_r = Rcpp::as<Rcpp::NumericVector>(Gx1(p1_r));
  Rcpp::NumericVector x2_r = Rcpp::as<Rcpp::NumericVector>(Gx2(p2_r));
  arma::vec x1 = Rcpp::as<arma::vec>(x1_r);
  arma::vec x2 = Rcpp::as<arma::vec>(x2_r);
  
  // hLeftX_vec equivalent, computed directly in C++
  double g = 2.0 * M_PI * sqrt(1.0 - rz*rz);
  double gprime = 2.0 * M_PI * 0.5 * (-2.0*rz) * pow(1.0 - rz*rz, -0.5);
  double h = 2.0 * (1.0 - rz*rz);
  double hsquare = h * h;
  double hprime = -4.0 * rz;
  
  arma::vec z1 = z.row(0).t();
  arma::vec z2 = z.row(1).t();
  
  arma::vec l = -(arma::square(z1) - 2.0*rz*z1%z2 + arma::square(z2));
  arma::vec lprime = 2.0 * z1 % z2;
  arma::vec numint = lprime*h - hprime*l;
  arma::vec numout = (numint/hsquare)*g - gprime;
  arma::vec hx = numout / g;
  
  // Compute derivative
  double covprime = arma::sum(x1 % x2 % hx) / K;
  double const_term = 1.0 / arma::prod(sdx);
  
  return const_term * covprime;
}