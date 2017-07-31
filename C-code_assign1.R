# ## Assignment1 C-programming
# part 2
# David Redmond
# ID 14207377
#
#  MLE calcualtion fora Cauchy distribution using 3 numerical optimization algorithms
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
incl<-'
  double ll_1(  NumericVector xx, double t )
  {
  NumericVector::iterator iter;
  double ans = 0;
  for ( iter = xx.begin(); iter < xx.end(); iter++) {
    ans += 2 * ( *iter - t ) / ( 1 + pow(( *iter - t ),2));
    }
  return ( ans );
  }
  double ll_2(  NumericVector xx, double t )
  {
  NumericVector::iterator iter;
  double ans = 0;
  for ( iter = xx.begin(); iter < xx.end(); iter++) {
    ans += 2 * (pow(( *iter - t ),2) - 2) / pow(1 + pow(( *iter - t ),2), 2) ;
    }
  return ( ans );
  }
  double minval ( NumericVector xx )
  {
    double minval = xx[0];
    NumericVector::iterator iter;
    for ( iter = xx.begin(); iter < xx.end(); iter++) {
      if(*iter < minval )
      minval = *iter;
    }
  return ( minval );
  }
double maxval ( NumericVector xx )
  {
    double maxval = xx[0];
    NumericVector::iterator iter;
    for ( iter = xx.begin(); iter < xx.end(); iter++) {
      if(*iter > maxval )
      maxval = *iter;
    }
  return ( maxval );
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
  double a  = maxval(xx);    // Call maxval from incl
  double b  = minval(xx);    // Call minval from incl
  double f_a ;              // Dummy variable
  double f_b ;              // dummy variable
  int j =0;
//
  while ((fabs(a-b) > err) && j <itr) {   // checking exit conditions
    f_a = ll_1(xx, a);
    f_b = ll_1(xx, b);
    if ( f_a * f_b < 0)
        b = (a + b )/2;
      else
        a = (a + b )/2 ;
      j++;
  }
  return( wrap((a+b)/2) );
 '
#
#compile,and link the above C++ function
bisectC <- cxxfunction( signature( vec = "numeric",
                                   iter= "integer",
                                   tol = "numeric" ),
                       body = body_bisectCpp,
                       includes = incl,
                       plugin = "Rcpp")
# Functional test########
# remove comments to apply block tests
yy <- dcauchy(runif(1:1000),1.5, 1)
## Parameters for testbench
itr=20000
tol = 1e-9
bisectC(yy, itr,tol)
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
double x[itr] = {0};
double root = 0;
x[1]  = maxval(xx);
x[0]  = minval(xx);
bool converged = false;
while (!converged) {
  for ( int k = 2; k < itr; k++) {
    if (!(ll_1(xx,x[k-1]) == ll_1(xx,x[k-2])))  // check for div by 0
    {
      x[k] = x[k-1] - ll_1(xx,x[k-1]) * (x[k-1] - x[k-2]) /  (ll_1(xx,x[k-1]) - ll_1(xx,x[k-2]));
      root = x[k];
    }
    converged = (fabs(x[k-2] - x[k-1]) < err) ;
    Rcout << "outside  if " << converged << std::endl;
  }
return(wrap(root));
}
'
#
#compile the above C++ function
secantC <- cxxfunction( signature( vec = "numeric",
                                   iter= "integer",
                                   tol = "numeric" ),
                        body = body_secantCpp,
                        includes = incl,
                        plugin = "Rcpp")
##
#
itr=20
tol = 1e-7
secantC(yy, itr,tol)
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
double x[itr] = {0};
x[1]  = maxval(xx);
x[0]  = minval(xx);
double root = 0.124; // dummy variable
double denom =0.1 ;  // dummy variable
for ( int k = 2; k < itr; k++)
      {
      denom =ll_2(xx,x[k-1]) ;
      if (!(denom==0) &&  (fabs(x[k] - x[k-1]) > err ) )
        x[k] = x[k-1] - (ll_1(xx,x[k-1])/ll_2(xx,x[k-1])) ;
      //Rcout << "outside  if " << root << std::endl;
      root = x[k];
      }
return(wrap(root));
'
#compile the above C++ function
NewRapC <- cxxfunction( signature( vec = "numeric",
                            iter= "integer",
                            tol = "numeric" ),
                            body = body_NewRCpp,
                            includes = incl,
                            plugin = "Rcpp")
#
#NewRapC(yy, itr,tol)

# supplied vector given below as part of testbench
#
xx <- c( 12.262307 , 10.281078 , 10.287090 , 12.734039 ,
         11.731881 , 8.861998 , 12.246509 , 11.244818 ,
         9.696278 , 11.557572 , 11.112531 , 10.550190 ,
         9.018438 , 10.704774 , 9.515617 , 10.003247 ,
         10.278352 , 9.709630 , 10.963905 , 17.314814)
## Parameters for testbench
itr=2000
tol = 1e-7
bisectC(xx, itr,tol)
secantC(xx, itr,tol)
NewRapC(xx, itr,tol)
library(rbenchmark)
benchmark( bisectC(xx, itr,tol), secantC(xx, itr,tol) ,NewRapC(xx, itr,tol), order = "relative")[,1:4]

#
# supplied vector given below as part of testbench#2
#
yy <- dcauchy(runif(1:1000),0.25, 1)
## Parameters for testbench
itr=2000
tol = 1e-8
bisectC(yy, itr,tol)
secantC(yy, itr,tol)
NewRapC(yy, itr,tol)
benchmark( bisectC(xx, itr,tol), secantC(xx, itr,tol) ,NewRapC(xx, itr,tol), order = "relative")[,1:4]
##
