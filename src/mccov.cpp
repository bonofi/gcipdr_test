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

// for test purposes only: must equal output of rmvnorm

// [[Rcpp::export]]
arma::mat test_rmvnorm_cpp_cholesky(int K, double rz)
{
  Rcpp::RNGScope rngScope;
  
  // Draw 2*K values — same total draws, but different reshaping
  Rcpp::NumericVector z_all = Rcpp::rnorm(2 * K, 0.0, 1.0);
  
  // Row-major fill: z[i,0] = z_all[2*i], z[i,1] = z_all[2*i+1]
  arma::vec z_col1(K);
  arma::vec z_col2(K);
  
  for (int i = 0; i < K; i++) {
    z_col1[i] = z_all[2*i];       // odd positions (0-indexed: 0,2,4,...)
    z_col2[i] = z_all[2*i + 1];   // even positions (0-indexed: 1,3,5,...)
  }
  
  // Apply Cholesky (upper triangular): [[1, rz],[0, sqrt(1-rz^2)]]
  double L_12 = rz;
  double L_22 = sqrt(1.0 - rz*rz);
  
  arma::vec z1 = z_col1;
  arma::vec z2 = z_col1 * L_12 + z_col2 * L_22;
  
  arma::mat out(K, 2);
  out.col(0) = z1;
  out.col(1) = z2;
  return out;
}

// [[Rcpp::export]]
arma::mat test_rmvnorm_cpp_eigen(int K, double rz)
{
  // Sigma = [[1, rz], [rz, 1]]
  arma::mat sigma(2, 2);
  sigma(0,0) = 1.0; sigma(0,1) = rz;
  sigma(1,0) = rz;  sigma(1,1) = 1.0;
  
  // Eigen decomposition (matches rmvnorm default method = "eigen")
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, sigma);
  
  // R matrix: t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0))))
  // = t(eigvec %*% diag(sqrt(eigval)) %*% t(eigvec))
  arma::mat R = arma::trans(
    eigvec * arma::diagmat(arma::sqrt(arma::clamp(eigval, 0.0, arma::datum::inf))) * arma::trans(eigvec)
  );
  
  // Draw n*p values row-major (byrow = TRUE)
  Rcpp::NumericVector z_all = Rcpp::rnorm(K * 2, 0.0, 1.0);
  
  // Fill matrix row-by-row (byrow = TRUE)
  arma::mat z(K, 2);
  for (int i = 0; i < K; i++) {
    z(i, 0) = z_all[2*i];
    z(i, 1) = z_all[2*i + 1];
  }
  
  // Apply R matrix: retval = z %*% R
  arma::mat retval = z * R;
  
  return retval;
}