#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List newtrap_one_cpp(Rcpp::Function fdist_and_fprime,
                           Rcpp::Nullable<Rcpp::Function> safecheck_arg,
                           double start, double tol = 0.01, int maxit = 50)
{
  double x = start;
  bool adjusted = false;
  int i = 0;
  int j = 0;
  
  bool has_safecheck = !safecheck_arg.isNull();
  
  // Initial evaluation — get both objective and derivative
  SEXP both_sexp = fdist_and_fprime(x);
  Rcpp::List both = Rcpp::as<Rcpp::List>(both_sexp);
  double out = Rcpp::as<double>(both["objective"]);
  
  while (out > tol || out < -tol)
  {
    // Use derivative from current evaluation
    double den = Rcpp::as<double>(both["derivative"]);
    
    // Safeguard: escape bottleneck if step is degenerate
    double ratio = out / den;
    if (ratio > 1e5 || ratio < -1e5)
    {
      x = x - 2.0 * out;
      i = i + 1;
      if (i > maxit)
      {
        Rcpp::warning("Search could not escape bottleneck: procedure has failed");
        break;
      }
      both_sexp = fdist_and_fprime(x);
      both = Rcpp::as<Rcpp::List>(both_sexp);
      out = Rcpp::as<double>(both["objective"]);
      continue;
    }
    
    // Newton-Raphson update
    x = x - (out / den);
    
    // Safecheck step (only if provided)
    if (has_safecheck)
    {
      Rcpp::Function safecheck = Rcpp::as<Rcpp::Function>(safecheck_arg);
      SEXP x_sexp = safecheck(x);
      x = Rcpp::as<double>(x_sexp);
      
      SEXP modified_sexp = Rf_getAttrib(x_sexp, Rf_install("modified"));
      bool modified = Rcpp::as<bool>(modified_sexp);
      
      if (modified)
      {
        j = j + 1;
        
        SEXP break_loop_sexp = Rf_getAttrib(x_sexp, Rf_install("break.loop"));
        bool break_loop = Rcpp::as<bool>(break_loop_sexp);
        
        if (break_loop)
        {
          adjusted = true;
          break;
        }
      }
    }
    
    // Evaluate at new x — gives us both objective and derivative for next iteration
    both_sexp = fdist_and_fprime(x);
    both = Rcpp::as<Rcpp::List>(both_sexp);
    out = Rcpp::as<double>(both["objective"]);
    
    i = i + 1;
    
    if (i > maxit)
    {
      Rcpp::warning("no zero found after max iteration");
      break;
    }
  }
  
  if (adjusted)
    Rcpp::warning("newtrap module: solution needed be adjusted to escape loop");
  
  Rcpp::NumericVector value(1);
  value[0] = x;
  value.attr("adjusted") = adjusted;
  
  return Rcpp::List::create(
    Rcpp::Named("value") = value,
    Rcpp::Named("iter")  = i
  );
}