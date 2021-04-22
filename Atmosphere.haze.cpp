#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
#include "constants.h"
#include "Atmosphere.haze.h"

using namespace cons;

Atmosphere::Atmosphere(int hazeint, vector<double> hazeparams)
{
  
  hazetype = hazeint;
  haze = hazeparams;
  std::string hazeFile;
  
  if(hazetype!=0){
    switch (hazetype){
      case 0: hazeFile = "Opacities/blank.r.dat"; break;
      case 1: hazeFile = "Opacities/zns.indices.dat"; break;
      case 2: hazeFile = "Opacities/polyacet.indices.dat"; break;
      case 3: hazeFile = "Opacities/tholin.indices.r.dat";
    }
    // read in indices of refraction
    std::ifstream hFilePrep(hazeFile.c_str());
    if(!hFilePrep)
      std::cout << "Haze Cross Section File Not Found" << std::endl;

    string line;
    hazeNumVals = 0;
    while(getline(hFilePrep,line)){
      hazeNumVals++;
    }

    elam = new double[hazeNumVals];
    einn = new double[hazeNumVals];
    eink = new double[hazeNumVals];
    hFilePrep.close();

    std::ifstream hFile(hazeFile.c_str());
    for(int i=0; i<hazeNumVals; i++){
      hFile >> elam[i] >> einn[i] >> eink[i];
    }
    ar = new double[20000];
    ai = new double[20000];
    br = new double[20000];
    bi = new double[20000];
  
    rsr = new double[20000];
    rsi = new double[20000];
    rsx = new double[20000];
    px = new double[20000];

    asf=0.;
    qe=0.;
    qa=0.;
    qs=0.;
    qp=0.;
    
    s1=0.;
    s2=0.;
    pext=0.;
    qext=0.;
    hFile.close();
  }
}

// General 2nd-order interpolation method used to interpolate complex indicies
// ***Need to fix this***
double Atmosphere::value(double x, double* sigma, double* e)
{
  double value;
  int n = hazeNumVals;
  if(x<e[0])
    {
      value = sigma[0] + (log(x)-log(e[0])) * (sigma[1]-sigma[0])/(log(e[1])-log(e[0]));
      return value;
    }
  else if(x>e[n-2])
    {
      value = sigma[n-1] + (log(x)-log(e[n-1])) * (sigma[n-1]-sigma[n-2])/(log(e[n-1])-log(e[n-2]));
      return value;
    }
  else
    {
      int ie=0;
      for(int i=0; i<n; i++)
	{
	  if(x<e[i])
	    {
	      ie = i-1;
	      break;
	    }
	}
      double e0 = log(e[ie]);
      double e1 = log(e[ie+1]);
      double e2 = log(e[ie+2]);
      double sig0 = log(sigma[ie]);
      double sig1 = log(sigma[ie+1]);
      double sig2 = log(sigma[ie+2]);
      double ee = log(x);
      
      double a0 = (ee-e1)*(ee-e2)/((e0-e1)*(e0-e2));
      double b0 = (ee-e0)*(ee-e2)/((e1-e0)*(e1-e2));
      double c0 = (ee-e1)*(ee-e0)/((e2-e1)*(e2-e0));
      double prob = a0*sig0 + b0*sig1 + c0*sig2;
      value = exp(prob);
      return value;
    }
}

/*
 * Calculates absorption and scattering efficiencies of spherical grains
 * Maximum size parameter: 1.e5
 * Inputs: size, wavelength, complex indices of refraction
 * Outputs:
 * qe = extinction efficiency = qa + qs
 * qs = scattering efficiency
 * qa = absorption efficiency
 * qp = radiation pressure efficiency
 * asf = asymmetrical factor
 */
void Atmosphere::MieLargeSzpara(double size, double wavelength, double xm, double ymm)
{
  int nx = 20000;
  double drc = 0.017453292519943;
  double rdc = 57.295779513082;
  double thb = 0.0;
  double thf = 180.0;
  double thi = 10.0;
  // xm = Re(refractive index)
  // ym = -Im(refractive index)
  double ym = -ymm;
  double x = 2.*pi*size/wavelength;  // size parameter

  double xn = sqrt(xm*xm + ym*ym)*x;
  //int nx = (int)(1.10*xn + 10.0);           // dwn recursion nx for all sizes
  double ro = 2.0*x*(xm-1.0);
  int i = smi1xs(x,nx,xm,ym);
  sm5msx(x,i);
  asf = (qe-qp)/qs;
  double xstops = x + 4.0*pow(x,1./3.) + 2.0; // cf Bohen-Huffman P480 (1983)
  int nv = (int)(std::max(xstops,xn)+15.);             // limit smi3dm array
  smi3sm(i,180.0);
  double qb = 2.*(s1*s1 + s2*s2)/x/x;         // back-scatter efficiency
  sm5pqs(x,i);
  int nthd = (int)((thf-thb)/thi + 1.0);

  double fi1=0.;
  double fi2=0.;
  for(int n=0; n<nthd; n++)
    {
      int thd = (int)(thb + (double)(n-1)*thi);
      smi3sm(i,thd);
      if(s1i==0.0) fi1 = pi/2.0;
      if(s1i!=0.0)
	{
	  fi1 = atan(s1r/fabs(s1i));
	  if(s1i<0.0) fi1 = pi-fi1;
	  if(fi1<0.0) fi1 = 2.0*pi - fabs(fi1);
	}
      double fi1d = rdc*fi1;
      if(s1i==0.0) fi2 = pi/2.0;
      if(s2i!=0.0)
	{
	  fi2 = atan(s2r/fabs(s2i));
	  if(s2i<0.0) fi2 = pi-fi2;
	  if(fi2<0.0) fi2 = 2.0*pi - fabs(fi2);
	}
      double fi2d = rdc*fi2;
      double s1s = s1*s1;
      double s2s = s2*s2;
      double sns = 0.50*(s1s+s2s);
      double gn = 4.0*sns/x/x;
      double pol = (s1s-s2s)/(s1s+s2s);
    }
  return;
}

// no explanation given
int Atmosphere::smi1xs(double x, int nx, double xm, double ym)
{
  double scale;
  double ctc = 1.e-14;              // edited from ctc=1.e-8
  double xmx = xm*x;                // m'x
  double ymx = ym*x;                // m"x
  double rp2 = xmx*xmx + ymx*ymx;   // (m'x)^2 + (m"x)^2
  double pnx = x/(2.*nx+3.);        // asymptotic value of pn(x) at n=nx+1
  double pnr = xmx/(2.*nx+3.);      // value of Re(pn(z)) at n=nx+1
  double pni = ymx/(2.*nx+3.);      // value of Im(pn(z)) at n=nx+1
  // dwn recursion of pn(z)=s(n)/s(n-1)=(pnr,pni) and an(z)=s'(n)/s(n)=(rsr(n),rsi(n))
  // also for real arguments: px(n) and rsx(n)
  for(int kk=0; kk<nx; kk++)
    {
      int n = nx-kk-1;
      double cn = (double)n+1.;
      double aln = (2.*cn+1.)*xmx/rp2-pnr; // Re(1./pn(z))
      double ben = (2.*cn+1.)*ymx/rp2+pni; // -Im(1./pn(z))
      ar[n] = -cn*xmx/rp2 + aln;   // Re(an(z))
      ai[n] = cn*ymx/rp2 - ben;    // Im(an(z))
      double pzd = aln*aln + ben*ben;
      pnr = aln/pzd;                // Re(pn(z))
      pni = ben/pzd;                // Im(pn(z))
      br[n] = (cn+1.)/x - pnx;      // s'(n)/s(n) for real x
      if(n==0) break;               // avoid computing px(1) irregular at x=n*pai
      pnx = x/(2.*cn+1.-pnx*x);
      bi[n] = pnx;                  // p(n)=s(n)/s(n-1) for real x
    }
  // generation of Mie coefficient to such n as determined by ctc (the convergence test constant)
  double snm1x = sin(x);
  double cnm1x = cos(x);
  double snx = 0.0;
  // for Rayleigh sizes x<0.1, power series employed for snx = sin(x)/x - cos(x)
  if(x<0.1) snx = x*x/3. - pow(x,4)/30. + pow(x,6)/840. - pow(x,8)/45360.;
  if(x>0.1) snx = snm1x/x - cnm1x;
  double cnx = cnm1x/x + snm1x;
  int nf = 0;
  for(int n=0; n<nx-1; n++)
    {
      double cn = (double)n+1.;
      double dcnx = cnm1x - cn*cnx/x;
      double dsnx = br[n]*snx;
      double annr = ar[n]*snx - xm*dsnx;    // Re(numerator of a(n))
      double anni = ai[n]*snx - ym*dsnx;    // Im(numerator of a(n))
      double ta1 = ar[n]*snx - ai[n]*cnx;
      double ta2 = ai[n]*snx + ar[n]*cnx;
      double andr = ta1 - xm*dsnx + ym*dcnx; // Re(denominator of a(n))
      double andi = ta2 - xm*dcnx - ym*dsnx; // Im(denominator of a(n))
      double andc = andr*andr + andi*andi;
      double bnnr = (xm*ar[n] - ym*ai[n]) * snx - dsnx; // Re(numerator of b(n))
      double bnni = (xm*ai[n] + ym*ar[n]) * snx;        // Im(numerator of b(n))
      double tb1 = ar[n]*snx - ai[n]*cnx;
      double tb2 = ar[n]*cnx + ai[n]*snx;
      double bndr = xm*tb1 - ym*tb2 - dsnx;               // Re(denominator of b(n))
      double bndi = xm*tb2 + ym*tb1 - dcnx;               // Re(denominator of b(n))
      double bndc = bndr*bndr + bndi*bndi;
      ar[n] = (annr*andr + anni*andi)/andc; // Re(a(n))
      ai[n] = (anni*andr - annr*andi)/andc; // Im(a(n))
      br[n] = (bnnr*bndr + bnni*bndi)/bndc; // Re(b(n))
      bi[n] = (bnni*bndr - bnnr*bndi)/bndc; // Im(b(n))
      double ti = ar[n]*ar[n] + ai[n]*ai[n] + br[n]*br[n] + bi[n]*bi[n];
      ti /= ar[0]*ar[0] + ai[0]*ai[0] + br[0]*br[0] + bi[0]*bi[0];
      if(nx!=3)
	{
	  if(ti-ctc<=0.)
	    {
	      nf = n;
	      break;
	    }
	}
      if(n-nx>=0)
	{
	  nf = n;
	  break;
	}
      snx *= bi[n+1];
      double cnm2x = cnm1x;
      cnm1x = cnx;
      cnx = (2.*cn+1.)*cnm1x/x - cnm2x;
      nf = n;
    }
  return nf;
}

// Sub sm5msx for Mie Q factors from smi1xs
void Atmosphere::sm5msx(double x, int i)
{
  double s1r = 0.0;
  double ssc = 0.0;
  double spr = 0.0;
  // i is number of iterations (I think)
  for(int n=0; n<i-1; n++)
    {
      double cn = (double)n+1.;
      double fs1r = 2.*cn+1.;
      double fssc = fs1r;
      double fspr1 = 2.*cn*(cn+2.)/(cn+1.);
      double fspr2 = 2.*(2.*cn+1.)/(cn*(cn+1.));
      s1r += fs1r * (ar[n]+br[n]);
      ssc += fssc * (ar[n]*ar[n] + ai[n]*ai[n] + br[n]*br[n] + bi[n]*bi[n]);
      spr += fspr1 * (ar[n]*ar[n+1] + ai[n]*ai[n+1] + br[n]*br[n+1] + bi[n]*bi[n+1]);
      spr += fspr2 * (ar[n]*br[n] + ai[n]*bi[n]);
    }

  double fc = 2.0/x/x;
  qe = fc*s1r;
  qs = fc*ssc;
  qa = qe-qs;
  qp = qe-fc*spr;
}

// smi3sm for Mie angular functions
void Atmosphere::smi3sm(int i, double thd)
{
  double scale;

  double th = thd * 0.0174532925;
  double z = cos(th);
  double pn1 = 0.0;
  double pn = 1.0;
  double s1r = 0.0;
  double s1i = 0.0;
  double s2r = 0.0;
  double s2i = 0.0;

  // i seems to be the number of iterations
  for(int n=0; n<i; n++)
    {
      double cn = (double)n+1.;
      double tn = cn*z*pn - (cn-1.)*pn1;
      double f = (2.*cn+1.)/(cn*cn+cn);
      s1r += f * (ar[n]*pn + br[n]*tn);
      s1i += f * (ai[n]*pn + bi[n]*tn);
      s2r += f * (ar[n]*tn + br[n]*pn);
      s2i += f * (ai[n]*tn + bi[n]*pn);
      double pn2 = pn1;
      pn1 = pn;
      pn = (z*(2.*cn+1.)*pn1 - (cn+1.)*pn2)/cn;
    }

  double s1s = s1r*s1r + s1i*s1i;
  double s2s = s2r*s2r + s2i*s2i;
  s1 = sqrt(s1s);
  s2 = sqrt(s2s);
}

// Sub sm5pqs for Mie P, Q comparison by smi1sx coefficients
void Atmosphere::sm5pqs(double x, int nm)
{
  double scale;

  double s1r = 0.0;
  double s1i = 0.0;

  for(int i=0; i<nm-1; i++)
    {
      double fs = 2.*i+1.;
      s1r += fs*(ar[i]+br[i]);
      s1i += fs*(ai[i]+bi[i]);
    }

  double fe = 2./x/x;
  qext = fe*s1r;
  pext = fe*s1i;
}
