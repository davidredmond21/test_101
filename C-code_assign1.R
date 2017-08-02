# ## Assignment1 C-programming
# part 2
# David Redmond
# ID 14207377
#
#  MLE calcualtion for a Cauchy distribution using 3 numerical optimization algorithms
#
#
library(Rcpp)
library(inline)
library(rbenchmark)
#
# incl contains 4 functions which will be used repeatedly by each algorithm
# ll_1 : is the 1-st deravitive log-likelihood of the Cauchy PDF taking "xx" input vector, and "t" as rate parameter
# ll_2 : is the 2-nd deravitive log-likelihood of the Cauchy PDF taking "xx" input vector, and "t" as  rate parameter
# minval & maxval: these are use to find initial estimates or intervals for the algorithms - the return the min/max-value of vector
# these functions use iterators
#

ll_r= function(xx,t){ 
  return( 2* sum(xx-t)/sum(1 + (xx-t)^2))
}

ll_r2= function(xx,t){ 
  return( 2* sum( (xx-t)^2-2)/sum((1 + (xx-t)^2)^2))
}
# incl contains 4 functions which will be used repeatedly by each algorithm
# ll_1 : is the 1-st deravitive log-likelihood of the Cauchy PDF taking "xx" input vector, and "t" as rate parameter
# ll_2 : is the 2-nd deravitive log-likelihood of the Cauchy PDF taking "xx" input vector, and "t" as  rate parameter
# minval & maxval: these are use to find initial estimates or intervals for the algorithms - the return the min/max-value of vector
# these functions use iterators
#
incl<-'
  double ll_1(   NumericVector & xx, const double & t)
  {
  NumericVector::iterator iter;
  double numer = 0;
  double denom = 1;
  for ( iter = xx.begin(); iter < xx.end(); iter++) {
    numer +=( *iter - t );
    denom +=( 1 + pow( (*iter - t),2));
    }
  return ( 2*numer/denom );
  }
  double ll_2(   NumericVector & xx, const double & t )
  {
  NumericVector::iterator iter;
  double numer = 0;
  double denom = 1;
  for ( iter = xx.begin(); iter < xx.end(); iter++) {
    numer = pow( ( *iter -t),2)-2 ;
    denom = pow( (1 + pow( ( *iter - t),2)),2); 
    }
  return ( 2*numer/denom );
  }
  double minval ( NumericVector & xx )
  {
    double minval = xx[0];
    NumericVector::iterator iter;
    for ( iter = xx.begin(); iter < xx.end(); iter++) {
      if(*iter < minval )
      minval = *iter;
    }
  return ( minval );
  }
double maxval ( NumericVector & xx )
  {
    double maxval = xx[0];
    NumericVector::iterator iter;
    for ( iter = xx.begin(); iter < xx.end(); iter++) {
      if(*iter > maxval )
      maxval = *iter;
    }
  return ( maxval );
  }
void check_valid_inputs( int & itr,const double & err)
  {
  if (itr < 3)
    stop("Error : stopping : Iteration < 3");
  if (err < 0 )
    Rcout << "Warning : Tolerance is less then 0" << std::endl;
  else if (err > 2)
    Rcout << "Warning : Tolerance is greater than 2" << std::endl;
  return;
}
'
# ######################## Biscet function ################
# inputs required are the vector xx, the mazimum number of loop iteration
# error tolerance for convergence
# Output is a single value
# Calls the ll_1, minval, &maxval functions from incl library
#
body_bisectCpp <- '
  NumericVector  xx(vec) ;    //create NumericVector xx from input vec
  int  itr   = as<int>(iter); //create integer for loop iterations iter
  double err = as<double>(tol); // create err : error tolerence from tol
  double a = as<double>(init_a); // create a : initial value
  double b = as<double>(init_b); // create a : initial value
  double c =  0 ;
  check_valid_inputs(itr, err); // Check inputs are valie, R checks the numeric vector
// check if 0 sent in and set to max/min in vector
  if (a==0) a = maxval(xx);   // Call maxval from incl
  if (b==0) b = minval(xx);  // Call minval from incl
  double f_a ;              // Dummy variable
  double f_b ;              // dummy variable
  bool converged = false;
  for ( int k = 1; (k < itr || !converged) ; k++) {
    f_a = ll_1(xx, a);
    f_b = ll_1(xx, b);
    c = ( a + b )/2 ;
    if ( f_a * ll_1(xx,c ) > 0)
        a = c ;
      else
        b = c ;
    converged = (fabs(a-b) < err) ;
  }
 // Rcout << "converged = " << converged << std::endl;
 return( wrap( c )  );
 '
#
#compile,and link the above C++ function
bisectC <- cxxfunction( signature( vec = "numeric",
                                   iter= "integer",
                                   tol = "numeric",
                                   init_a = "numeric",
                                   init_b = "numeric"),
                       body = body_bisectCpp,
                       includes = incl,
                       plugin = "Rcpp")
# Functional block level test########
# remove comments to apply block tests
## Parameters for testbench
##itr=2000
##tol = 1e-9
##bisectC(vec=yy, iter=itr,init_a=inita, init_b = initb, tol)
#
#################################Secant method##########################
# Secant method
#  x[n] = x[n-1] - f(x[n-1]) * (x[n-1] - x[n-2])/(f(x[n-1] - f(x[n-2]));
# inputs required are the vector xx, the maximum number of loop iteration
# error tolerance for convergence
# Output is a single value for theta
# Calls the ll_1, minval, &maxval functions from incl library
# #######################################################################
body_secantCpp <- '
// inputs are the vector xx, the tolerence for convergence tol, the range min, max
NumericVector  xx(vec) ;    //create NumericVector xx from input vec
int  itr   = as<int>(iter); //create integer for loop iterations iter
double err = as<double>(tol); // create err : error tolerence from tol
double x_0 = as<double>(init_a); // create a : initial value
check_valid_inputs(itr, err); // Check inputs are valie, R checks the numeric vector
double x_1 = x_0 + 0.01;
double root = 0;
double denom =1.0; 
double converged = 100;  // large iniitial value
for ( int k = 1; ( (k < itr) && (converged > err) && (denom !=0)) ; k++) {
     denom = ll_1(xx,x_1)- ll_1(xx,x_0);
     // Rcout << "denom  " << k << " " << denom  << std::endl;
     root = x_1 -  (x_1 - x_0) * ll_1(xx,x_1) / denom ;
     x_0 = x_1;
     x_1 = root;
     converged = fabs(x_0 - x_1) ;
    }
if( fabs(converged) < 1) {
    // Rcout << "convergence OK " << converged << std::endl;
    return(wrap(root));
}    
else
    Rcout << "Secant algorithm failed to convergence -> no root"  << std::endl;

'
#
#compile the above C++ function
secantC <- cxxfunction( signature( vec = "numeric",
                                   iter= "integer",
                                   tol = "numeric",
                                   init_a = "numeric"
                                  ),
                        body = body_secantCpp,
                        includes = incl,
                        plugin = "Rcpp")
##
#
#itr=20
#tol = 1e-7
#secantC(yy, itr,tol)
##################Newton Raphson #########################################
# NewtonRaphson
#     x(n+1) = x(n) - f'(x(n))/ f"(x(n))
# ll_1 is the 1st deravitive f'(c)
# ll_2 is the 2nd deravitive f"(c)
# inputs required are the vector xx, the maximum number of loop iteration
# error tolerance for convergence
# Output is a single value for theta
# Calls the ll_1, ll_2, minval, &maxval functions from incl library
#

body_NewRCpp <- '
// inputs are the vector xx, the tolerence for convergence tol, the range min, max
NumericVector  xx(vec) ;    //create NumericVector xx from input vec
int  itr   = as<int>(iter); //create integer for loop iterations iter
double err = as<double>(tol); // create err : error tolerence from tol
double x_0 = as<double>(init_a); // create a : initial value
check_valid_inputs(itr, err); // Check inputs are value, R checks the numeric vector
double root = 0.124; // dummy variable
double denom =1.0 ;  // dummy variable
double converged = 20;  // initial large value for converged

for ( int k =1; ( (k < itr) && (converged > err) && (denom !=0 )) ; k++) {
      denom = ll_2(xx,x_0) ;
      root = x_0 - ll_1(xx,x_0)/ll_2(xx,x_0) ;
      converged = fabs(x_0 - root) ;
      x_0 = root ; 
    }
    if( fabs(converged) < 1) {
      return(wrap(root));
    }
    else
      Rcout << "NewtonRaphson algorithm failed to convergence -> no root"  << std::endl;
'
#compile the above C++ function
NewRapC <- cxxfunction( signature( vec = "numeric",
                            iter= "integer",
                            tol = "numeric",
                            init_a = "numeric" ),
                            body = body_NewRCpp,
                            includes = incl,
                            plugin = "Rcpp")
#
# supplied vector given below as part of testbench
#
xx <- c( 12.262307 , 10.281078 , 10.287090 , 12.734039 ,
         11.731881 , 8.861998 , 12.246509 , 11.244818 ,
         9.696278 , 11.557572 , 11.112531 , 10.550190 ,
         9.018438 , 10.704774 , 9.515617 , 10.003247 ,
         10.278352 , 9.709630 , 10.963905 , 17.314814)
## Parameters for testbench
itr=200
tol = 1e-9
inita = 0.0
initb = 11.0
bisectC(xx, itr, tol, inita, initb)
inita = bisectC(xx, itr, tol, inita, initb)
secantC(xx, itr, tol, inita)
NewRapC(xx, itr, tol, inita)
#
benchmark(bisectC(xx, itr, tol, inita, initb),
          secantC(xx, itr, tol,inita),
          NewRapC(xx, itr, tol,inita),
          order = "relative")[,1:4]
# another vector type Cauchy
## using a combination of Bisecnt and Newton

yy <- dcauchy(seq(from =-3, to =3, by= 0.001),6, 1)
approx= bisectC(yy, 4, 1e-5, 0, +1)
secantC(yy, 300, 1e-9, approx)
NewRapC(yy, 300, 1e-9, approx)
