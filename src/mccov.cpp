#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Internal helper: replicates rmvnorm(K, sigma=Sigmaz, method="eigen") exactly
// Uses R's RNG stream — respects set.seed()
void rmvnorm_eigen_cpp(int K, double rz, arma::vec& z1, arma::vec& z2)
{
  // Build Sigma = [[1, rz], [rz, 1]]
  arma::mat sigma(2, 2);
  sigma(0,0) = 1.0; sigma(0,1) = rz;
  sigma(1,0) = rz;  sigma(1,1) = 1.0;
  
  // Eigen decomposition — matches rmvnorm method="eigen"
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, sigma);
  
  // R = t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0))))
  arma::mat R = arma::trans(
    eigvec *
      arma::diagmat(arma::sqrt(arma::clamp(eigval, 0.0, arma::datum::inf))) *
      arma::trans(eigvec)
  );
  
  // Draw K*2 values row-major (byrow=TRUE matches rmvnorm)
  Rcpp::NumericVector z_all = Rcpp::rnorm(K * 2, 0.0, 1.0);
  
  // Fill matrix row-by-row
  arma::mat z(K, 2);
  for (int i = 0; i < K; i++) {
    z(i, 0) = z_all[2*i];
    z(i, 1) = z_all[2*i + 1];
  }
  
  // Apply R matrix: retval = z %*% R
  arma::mat retval = z * R;
  
  // Extract columns
  z1 = retval.col(0);
  z2 = retval.col(1);
}

// [[Rcpp::export]]
double mccovx1x2_cpp(double rz, Rcpp::Function Gx1, Rcpp::Function Gx2,
                     double rx, arma::vec meanx, arma::vec sdx,
                     bool pNorm_1, bool pNorm_2, int K)
{
  // Generate correlated normals matching rmvnorm exactly
  arma::vec z1, z2;
  rmvnorm_eigen_cpp(K, rz, z1, z2);
  
  // Apply pnorm or identity
  arma::vec p1 = pNorm_1 ? arma::normcdf(z1) : z1;
  arma::vec p2 = pNorm_2 ? arma::normcdf(z2) : z2;
  
  // Call R marginal inverse CDFs
  arma::vec x1 = Rcpp::as<arma::vec>(Gx1(Rcpp::wrap(p1)));
  arma::vec x2 = Rcpp::as<arma::vec>(Gx2(Rcpp::wrap(p2)));
  
  // Compute covariance
  double covx1x2 = arma::sum(x1 % x2) / K;
  double res = (covx1x2 - arma::prod(meanx)) / arma::prod(sdx);
  
  return res - rx;
}

// [[Rcpp::export]]
double mccovx1x2prime_cpp(double rz, Rcpp::Function Gx1, Rcpp::Function Gx2,
                          arma::vec sdx, bool pNorm_1, bool pNorm_2, int K)
{
  // Generate correlated normals matching rmvnorm exactly
  arma::vec z1, z2;
  rmvnorm_eigen_cpp(K, rz, z1, z2);
  
  // Apply pnorm or identity
  arma::vec p1 = pNorm_1 ? arma::normcdf(z1) : z1;
  arma::vec p2 = pNorm_2 ? arma::normcdf(z2) : z2;
  
  // Call R marginal inverse CDFs
  arma::vec x1 = Rcpp::as<arma::vec>(Gx1(Rcpp::wrap(p1)));
  arma::vec x2 = Rcpp::as<arma::vec>(Gx2(Rcpp::wrap(p2)));
  
  // hLeftX computation — vectorized, matches hLeftX_vec exactly
  double g       = 2.0 * M_PI * sqrt(1.0 - rz*rz);
  double gprime  = 2.0 * M_PI * 0.5 * (-2.0*rz) * pow(1.0 - rz*rz, -0.5);
  double h       = 2.0 * (1.0 - rz*rz);
  double hsquare = h * h;
  double hprime  = -4.0 * rz;
  
  arma::vec l      = -(arma::square(z1) - 2.0*rz*(z1%z2) + arma::square(z2));
  arma::vec lprime = 2.0 * z1 % z2;
  arma::vec numint = lprime*h - hprime*l;
  arma::vec numout = (numint/hsquare)*g - gprime;
  arma::vec hx     = numout / g;
  
  // Compute derivative
  double covprime   = arma::sum(x1 % x2 % hx) / K;
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