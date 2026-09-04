// // order_minimal.cpp contains C++ code used in order_cppWrapper.R.
// // Copyright (C) 2018 Federico Bonofiglio

// // This file is part of gcipdr.

//     // gcipdr is free software: you can redistribute it and/or modify
//     // it under the terms of the GNU General Public License as published by
//     // the Free Software Foundation, either version 3 of the License, or
//     // (at your option) any later version.

//     // gcipdr is distributed in the hope that it will be useful,
//     // but WITHOUT ANY WARRANTY; without even the implied warranty of
//     // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     // GNU General Public License for more details.

//     // You should have received a copy of the GNU General Public License
//     // along with gcipdr.  If not, see <https://www.gnu.org/licenses/>.

 
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>


using namespace Rcpp;
using namespace arma;




ivec rankC(const vec& X){  // no silly while loops !!!
  
  NumericVector x = Rcpp::wrap(X);
  NumericVector sorted = clone(x).sort();
  IntegerVector minranking = match(x, sorted);

  ivec ranked(minranking); 
  
  return ranked; 

}




// should be equivalent to rjm min in R

double rjmC( const vec& x,  const vec& y){

  return sum( rankC(x)%rankC(y) );

}




// should be equivalent to spcor ties min in R


double spcorC( const vec& x, const vec& y){ 

  vec rx = conv_to<vec>::from( rankC(x)) ;
  vec ry = conv_to<vec>::from( rankC(y)) ;
  return as_scalar(cor(rx,ry)) ;

}



// should be equivalent to joint mom in R

double jmC( const vec& x, const vec& y){

  return sum( x%y );
 
}

// pearson corr
double pearcorrC( const vec& x, const vec& y ){ 
   return as_scalar( cor(x ,y ) );   

}


// absolute value for double 

double absv(double x){

  return sqrt( pow( x, 2 )  );

}


// function pointer syntax from Rcpp Gallery Dirk Edelbeutel:  
// http://gallery.rcpp.org/articles/passing-cpp-function-pointers/ 

// define a pointer ...to be later referenced to: here funcpointer is a deferenced function.

typedef double (*funcpointer)( const vec& x, const vec& y);


// next I use XPtr to assign the above function pointer to the adress (&) of my functions. In practice I redirect the pointer to my functions, while assigning them labels.

XPtr<funcpointer> adressmyfuncs(std::string flabel){

  if (flabel == "rank.corr")
    return( XPtr<funcpointer>( new funcpointer ( &spcorC  ) )     );
  else if (flabel == "joint.rank")
    return( XPtr<funcpointer>( new funcpointer ( &rjmC  ) )     );
  else if (flabel == "joint.moment")
    return( XPtr<funcpointer>( new funcpointer ( &jmC  ) )     );
  else if (flabel == "moment.corr")
    return( XPtr<funcpointer>( new funcpointer ( &pearcorrC  ) )     );
  else
    return XPtr<funcpointer>(R_NilValue);

}



// zero finder in Rüschendorf function

// [[Rcpp::export]]
List findzerocorr(const arma::vec& xell, 
                  const arma::vec& x, 
                  const double test, 
                  const std::string rf_type){

  // call type of function to apply

  XPtr<funcpointer> xpfun = adressmyfuncs(rf_type);
  funcpointer fun = *xpfun;
  mat zero(x.n_elem, x.n_elem);
  mat corr(x.n_elem, 1);
  mat out(x.n_elem, 2);
  out.col(0) = xell;
  double outcorr ;
  double check;
 
  arma::uword i = 0;  // instead of int i = 0;

  check = 1;  // intial check value

  while( absv(check) > test ){
    vec z = shuffle(x);
    check = fun( xell, z );

   // update output
    out.col(1) =  z ;
    outcorr  =  check ; 

    // save intermediate results
      
    zero.col(i) = z;
    corr.row(i) = check;

    i = i + 1 ;

    if( i==x.n_elem ){
      mat tmat(x.n_elem, 1);
      tmat.fill(test);
      vec diffs =  sqrt( square( corr - tmat ) );
      uvec pick = find( diffs == min(diffs) ) ;
      out.col(1) = zero.col( pick(0)  );
      mat res  =  corr.row( pick(0) ) ;
      outcorr = as_scalar(res);
      break;
    }
  }

  return List::create(
		      Named("vec") =  out,  
		      Named("cor") =  outcorr,
		      Named("iter") =  i
			   
                           ); 


}








