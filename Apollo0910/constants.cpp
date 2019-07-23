/*
 * List of constants for the Apollo Retrieval Code
 * Alex R. Howe
 */
#include <math.h>
#include <stdio.h>
#include "constants.h"

using namespace cons;

// all in cgs units, so don't screw it up
// except wavelengths all in microns
const double cons::G       = 6.67428e-8;
const double cons::c       = 2.99792458e10;
const double cons::pi      = 3.1415926535897932;
const double cons::h       = 6.62607004e-27;
const double cons::hc      = 1.98644568e-16;
const double cons::k       = 1.38064852e-16;
const double cons::N_A     = 6.0221409e23;
const double cons::M_EARTH = 5.9736e27;
const double cons::R_EARTH = 6.371e8;

double cons::gauss(double x, double m, double s)
{
  double in2pi = 0.39894228;
  double a = (x-m)/s;
  return in2pi / s * exp(-0.5*a*a);
}

double cons::blackbodyN(double T, double freq)
{
  return (2*h*pow(freq,3)/c/c) / (exp(h*freq/k/T)-1);
}

double cons::blackbodyL(double T, double wave)
{
  return (2*h*c*c/pow(wave,5)) / (exp(h*c/wave/k/T)-1);
}

double cons::expint(int n, double x)
{
  int maxit = 100;
  double euler = 0.57721566349;
  double fpmin = 1.e-30;
  double eps = 1.e-7;
  
  void nrerror(char error_text[]);
  int i,ii,nm1;
  double a,b,c,d,del,fact,h,psi,ans;

  nm1 = n-1;
  if(n<0 || x<0. || (x==0. && (n==0 || n==1))){
    printf("Exponential integral failed.");
    return 0.;
  }
  else{
    if(n==0) ans = exp(-x)/x;
    else{
      if(x==0.) ans = 1.0/nm1;
      else{
	if(x>1.0){
	  b = x+n;
	  c = 1./fpmin;
	  d = 1./b;
	  h = d;
	  for(i=1; i<maxit; i++){
	    a = -i*(nm1+i);
	    b += 2.0;
	    d = 1./(a*d+b);
	    c = b+a/c;
	    del = c*d;
	    h *= del;
	    if(fabs(del-1.)<eps){
	      ans = h*exp(-x);
	      return ans;
	    }
	  }
	  printf("Exponential integral failed.");
	  return 0.;
	}
	else{
	  ans = (nm1!=0 ? 1./nm1 : -log(x)-euler);
	  fact = 1.;
	  for(i=1; i<=maxit; i++){
	    fact *= -x/i;
	    if(i!=nm1) del = -fact/(i-nm1);
	    else{
	      psi = -euler;
	      for(ii=1; ii<=nm1; ii++) psi += 1./ii;
	      del = fact*(-log(x)+psi);
	    }
	    ans += del;
	    if(fabs(del)<fabs(ans)*eps) return ans;
	  }
	  printf("Exponential integral failed.");
	  return 0.;
	}
      }
    }
  }
  return ans;
}
