#include <omp.h>
#include "specincld.h"
#include "constants.h"
#include "Atmosphere.h"
#include "Planet_layer.h"
//#include <valgrind/memcheck.h>

using namespace cons;

// Constructor; includes all the variables that don't change between models
Planet::Planet(int modenum, vector<double> waves, int cloudnum, int hazenum, vector<int> mollist)
{
  npress = 18;
  ntemp = 35;
  nwave = 21205;
  nspec = mollist.size();
  mode = modenum;
  hazetype = hazenum;
  cloudmod = cloudnum;
  wavens = waves;
  
  degrade = round(10000.*log(wavens[0]/wavens[1]));
  ntable = (int)(nwave/degrade);
  if(nspec==0){
    mastertable = vector<vector<vector<vector<double> > > >(npress,vector<vector<vector<double> > >(ntemp,vector<vector<double> >(ntable,vector<double>(1))));
  }
  else{
    mastertable = vector<vector<vector<vector<double> > > >(npress,vector<vector<vector<double> > >(ntemp,vector<vector<double> >(ntable,vector<double>(nspec))));
  }
  
  lmin = 0.6;
  lmax = 5.0;
  res = (lmax-lmin)/nwave;

  rp = 1.0;
  teff = 0.;
  
  vector<double> temphaze(4,0);
  atmosphere = new Atmosphere(hazetype,temphaze);
  readopac(mollist);
}
// end constructor

// Reads in the variables that define a particular model and sets the parameters
void Planet::setParams(vector<double> plparams, vector<double> abund, vector<double> tpprofile)
{
  rp = plparams[0]*R_EARTH;
  grav = pow(10,plparams[1]);
  double pressure = pow(10,plparams[2]);
  tstar = plparams[3];
  rs = plparams[4]*6.96e10;
  mu = plparams[5]/N_A;
  rxsec = plparams[6];
  minP = plparams[7];
  maxP = plparams[8];
  sma = plparams[9]*1.496e13;
  
  nlevel = tpprofile.size();
  nlayer = nlevel-1;      // Layers between the T-P points.
  tpprof = vector<double>(nlevel,0);
  hprof = vector<double>(nlevel,0);
  
  // Flipped around to work with Mark Marley's radiative transfer algorithm.
  for(int i=0; i<nlevel; i++){
    tpprof[i] = tpprofile[nlevel-1-i];
  }
  
  haze = vector<double>(4,0);
  if(hazetype!=0){
    if(cloudmod==2){
      for(int i=10; i<14; i++){
	haze[i-10] = pow(10,plparams[i]);
      }
    }
    if(cloudmod==3){
      haze[0] = pow(10,plparams[10]);
      haze[1] = pow(10,plparams[11]);
      haze[2] = plparams[12];
      haze[3] = plparams[13];
    }
  }
  
  mp = grav/G*rp*rp;
  getProfile(tpprof);
  
  // hmax changed to the 1 mubar level; can change this
  hmin = getH(pressure);
  hmax = hprof[0];
  setWave(nwave,rxsec,wavens,abund);
  return;
}
// end setParams

// Calls the appropriate method given the spectrum type
vector<double> Planet::getSpectrum(int streams){
  vector<double> tdepth;
  if(mode<=1 && streams==2) tdepth = getFlux(wavens);
  else if(mode<=1 && streams==1) tdepth = getFluxOneStream(wavens);
  else if(mode==2) tdepth = transFlux(rs,wavens);
  else printf("Error: invalid settings in getSpectrum.");
  return tdepth;
}

// Reads in the opacities tables
void Planet::readopac(vector<int> mollist){
  printf("Reading in cross sections.\n");
  string specfile;

  string gaslist[13] = {"h2","h2o","ch4","co","co2","nh3","h2s","Burrows_alk","Lupu_alk","crh","feh","tio","vo"};
  
  // Gray atmosphere for testing purposes.
  if(nspec==0){
    int x;
    for(int j=0; j<npress; j++){
      for(int k=0; k<ntemp; k++){
	for(int l=0; l<nwave; l++){
	  x = (int)(l/degrade);
	  mastertable[j][k][x][0] = 6.629e-24;
	}
      }
    }
  }

  else{
    for(int i=0; i<nspec; i++){
      int index = mollist[i];
      specfile = "Opacities/" + gaslist[index] + "wide.dat";
      printf("%s\n",gaslist[index].c_str());
    
      ifstream opacin(specfile.c_str());
      if(!opacin) cout << "Opacity File Not Found" << std::endl;
      int x;
      double val;

      for(int j=0; j<npress; j++){
	for(int k=0; k<ntemp; k++){
	  for(int l=0; l<nwave; l++){
	    x = (int)(l/degrade);
	    opacin >> val;
	    if(isnan(val) || isinf(val) || val < 1.e-50) val = 1.e-50;
	    if(x<ntable) mastertable[j][k][x][i] += val/(double)degrade;
	  }
	}
      }
    }
  }
  
  printf("Finished reading cross sections.\n");
  return;
}

// Calls the methods to create the optical depth tables for this atmosphere
void Planet::setWave(int npoints, double rxsec, vector<double> wavens, vector<double> abund){

  for(int i=1; i<nspec; i++){
    abund[i] = pow(10,abund[i]);
  }
  getOpacProf(rxsec, wavens, abund);
  if(mode<=1) getTauProf(wavens);
  if(mode==2) transTauProf(wavens);
}

// Retrieves the effective temperature for the output file
double Planet::getTeff()
{
  return teff;
}

// Retrieves the pressure at a given altitude
double Planet::getP(double height)
{
  int length = hprof.size();
  if(height>hprof.back()){
    double pr = prprof.back() * exp(G*mp*mu/k/tprof.back() * (1/(rp+height)-1/(rp+hprof.back())));
    return pr;
  }
  
  int jh=0;
  double checkh=-1.;
  for(int i=0; i<length; i++){
    checkh = hprof[i] - height;
    if(checkh>0.){
      jh = i-1;
      break;
    }
  }
  double dhi = (height-hprof[jh])/(hprof[jh+1]-hprof[jh]);
  if(jh<0){
    jh=0.;
    dhi=0.;
  }
  double pr = prprof[jh] + dhi*(prprof[jh+1]-prprof[jh]);
  return pr;
}
// end getP

// retrieves the temperature at a given altitude
double Planet::getT(double height)
{
  int length = hprof.size();
  if(height>hprof[0]) return tprof[0];
  
  int jh=0;
  double checkh=-1.;
  for(int i=0; i<length; i++){
    checkh = hprof[i] - height;
    if(checkh<0.){
      jh = i-1;
      break;
    }
  }
  double dhi = (height-hprof[jh])/(hprof[jh+1]-hprof[jh]);
  if(jh<0){
    jh=0.;
    dhi=0.;
  }
  double T = tprof[jh] + dhi*(tprof[jh+1]-tprof[jh]);
  return T;
}
// end getT

// retrieves the altitude for a given pressure
double Planet::getH(double pr)
{
  int length = prprof.size();
  double pl = log10(pr);
  double deltap = (pl-minP)/(maxP-minP)*(length-1);
  int jp = (int)deltap;
  double dpi = (pl-log10(prprof[jp]))/(log10(prprof[jp+1])-log10(prprof[jp]));
  if(jp<0){
    jp=0.;
    dpi=0.;
  }
  if(jp>=length-1){
    return hprof.back();
  }
  return hprof[jp] + dpi*(hprof[jp+1]-hprof[jp]);
}
// end getH


// Two-stream radiative transfer function.
// Version used by Ben Burningham.
vector<double> Planet::getFlux(vector<double> wavens)
{
  // 8 Gauss points
  int nGauss = 8;
  vector<double> gaussPoints(nGauss,0);
  vector<double> gaussWeights(nGauss,0);

  // Cosines of angles to calculate the streams
  gaussPoints[0] = 0.0446339553;
  gaussPoints[1] = 0.1443662570;
  gaussPoints[2] = 0.2868247571;
  gaussPoints[3] = 0.4548133152;
  gaussPoints[4] = 0.6280678354;
  gaussPoints[5] = 0.7856915206;
  gaussPoints[6] = 0.9086763921;
  gaussPoints[7] = 0.9822200849;

  // Weights for the streams
  gaussWeights[0] = 0.0032951914;
  gaussWeights[1] = 0.0178429027;
  gaussWeights[2] = 0.0454393195;
  gaussWeights[3] = 0.0791995995;
  gaussWeights[4] = 0.1060473494;
  gaussWeights[5] = 0.1125057995;
  gaussWeights[6] = 0.0911190236;
  gaussWeights[7] = 0.0445508044;

  double ubari = 0.5;    // Not sure what this does.
  int nl2 = 2*nlayer;    // Needed for the DSOLVER subroutine.
  double rsf = 0.;       // "surface" reflectivity, can set to zero

  vector<double> fdown(nlevel,0);
  vector<double> fup(nlevel,0);

  vector<double> alpha(nlayer,0);
  vector<double> lamda(nlayer,0);
  vector<double> gama(nlayer,0);
  vector<double> b0(nlayer,0);
  vector<double> b1(nlayer,0);
  vector<double> bdiff(nlayer,0);

  vector<double> cp(nlayer,0);
  vector<double> cm(nlayer,0);
  vector<double> cpm1(nlayer,0);
  vector<double> cmm1(nlayer,0);

  vector<double> e1(nlayer,0);
  vector<double> e2(nlayer,0);
  vector<double> e3(nlayer,0);
  vector<double> e4(nlayer,0);

  vector<double> gg(nlayer,0);
  vector<double> hh(nlayer,0);
  vector<double> xj(nlayer,0);
  vector<double> xk(nlayer,0);

  vector<double> af(nl2,0);
  vector<double> bf(nl2,0);
  vector<double> cf(nl2,0);
  vector<double> df(nl2,0);
  vector<double> as(nl2,0);
  vector<double> ds(nl2,0);
  vector<double> xki(nl2,0);
  vector<double> xk1(nlayer,0);
  vector<double> xk2(nlayer,0);

  vector<double> alpha1(nlayer,0);
  vector<double> alpha2(nlayer,0);
  vector<double> sigma1(nlayer,0);
  vector<double> sigma2(nlayer,0);

  vector<double> em1(nlevel,0);
  vector<double> em2(nlevel,0);
  vector<double> em3(nlevel,0);
  vector<double> epp(nlevel,0);
  vector<double> fpt(nlevel,0);
  vector<double> fmt(nlevel,0);

  vector<double> totalflux(wavens.size(),0);

  // cosbar is layer asymmetry parameter with layer 0 on top
  // zero for most circumstances; may be nonzero for high cloud opacity
  vector<double> cosbar(nlayer,0);

  for(int i=0; i<wavens.size(); i++){
    double wavelength = 1./wavens[i];

    // computes blackbody emission by layer from top to bottom
    for(int j=0; j<nlayer; j++){
      alpha[j] = sqrt( (1.-w0[i][j])/(1.-w0[i][j]*cosbar[j]) );
      lamda[j] = alpha[j]*(1.-w0[i][j]*cosbar[j])/ubari;
      gama[j] = (1.-alpha[j])/(1+alpha[j]);
      double term = 0.5/(1.-w0[i][j]*cosbar[j]);
      b0[j] = blackbodyL(tprof[j],wavelength);
      b1[j] = blackbodyL(tprof[j+1],wavelength);
      bdiff[j] = (b1[j]-b0[j])/taulayer[i][j];
      
      if(taulayer[i][j]<1.e-6){
	b0[j] = 0.5*(blackbodyL(tprof[j],wavelength)+blackbodyL(tprof[j+1],wavelength));
	bdiff[j] = 0.;
      }
      
      cp[j] = b0[j] + bdiff[j]*taulayer[i][j] + bdiff[j]*term;
      cm[j] = b0[j] + bdiff[j]*taulayer[i][j] - bdiff[j]*term;
      cpm1[j] = b0[j] + bdiff[j]*term;
      cmm1[j] = b0[j] - bdiff[j]*term;
    } // End of blackbody loop.

    // Computes e-coefficients (whatever those are), top to bottom
    for(int j=0; j<nlayer; j++){
      double ep = exp( min(lamda[j]*taulayer[i][j], 35.) );
      e1[j] = ep + gama[j]/ep;
      e2[j] = ep - gama[j]/ep;
      e3[j] = gama[j]*ep + 1./ep;
      e4[j] = gama[j]*ep - 1./ep;
    }

    double tautop = taulayer[i][0];
    double btop = (1. - exp(-tautop/ubari))*blackbodyL(tprof[0],wavelength);
    double bsurf = blackbodyL(tprof[nlevel-1],wavelength);
    double bottom = bsurf + b1[nlayer-1]*ubari;

    // DSolver subroutine to compute xk1 and xk2.
    // Computes a,b,c,d coefficients first, top to bottom
    // Then as and ds, *bottom to top*
    // The xk coefficients, top to bottom
    af[0] = 0.;
    bf[0] = gama[0] + 1.;
    cf[0] = gama[0] - 1.;
    df[0] = btop - cmm1[0];
    
    int nn=0;
    int lm1 = nl2-1;

    // even indices
    for(int ii=2; ii<lm1; ii+=2){
      af[ii] = 2.*(1.-gama[nn]*gama[nn]);
      bf[ii] = (e1[nn]-e3[nn])*(1.+gama[nn+1]);
      cf[ii] = (e1[nn]+e3[nn])*(gama[nn+1]-1.);
      df[ii] = e3[nn]*(cpm1[nn+1]-cp[nn]) + e1[nn]*(cm[nn]-cmm1[nn+1]);
      nn++;
    }
    
    nn=0;
    int lm2 = nl2-2;

    // odd indices
    for(int ii=1; ii<lm2; ii+=2){
      af[ii] = (e1[nn]+e3[nn])*(gama[nn+1]-1.);
      bf[ii] = (e2[nn]+e4[nn])*(gama[nn+1]-1.);
      cf[ii] = 2.*(1.-gama[nn+1]*gama[nn+1]);
      df[ii] = (gama[nn+1]-1.)*(cpm1[nn+1]-cp[nn]) + (1.-gama[nn+1])*(cm[nn]-cmm1[nn+1]);
      nn++;
    }

    af[nl2-1] = e1[nlayer-1] - rsf*e3[nlayer-1];
    bf[nl2-1] = e2[nlayer-1] - rsf*e4[nlayer-1];
    cf[nl2-1] = 0.;
    df[nl2-1] = bottom - cp[nlayer-1] + rsf*cm[nlayer-1]; // original says bsurf, but was called with bottom

    // DTRIDGL subroutine to compute the necessary xki array
    as[nl2-1] = af[nl2-1]/bf[nl2-1];
    ds[nl2-1] = df[nl2-1]/bf[nl2-1];
    for(int ii=2; ii<nl2; ii++){
      double xx = 1./(bf[nl2-ii] - cf[nl2-ii]*as[nl2-ii+1]);
      as[nl2-ii] = af[nl2-ii]*xx;
      ds[nl2-ii] = (df[nl2-ii] - cf[nl2-ii]*ds[nl2-ii+1])*xx;
    }
    xki[0] = ds[0];
    for(int ii=1; ii<nl2; ii++){
      xki[ii] = ds[ii] - as[ii]*xki[ii-1];
    }
    // End of DTRIDGL subroutine
    
    for(int nn=0; nn<nlayer; nn++){
      xk1[nn] = xki[2*nn] + xki[2*nn+1];
      xk2[nn] = xki[2*nn] - xki[2*nn+1];
      if(xk2[nn]!=0. && fabs(xk2[nn]/xk[2*nn] < 1.e-30)) xk2[nn] = 0.;
    }
    // End of DSolver subroutine.

    for(int ng=0; ng<nGauss; ng++){
      double ugauss = gaussPoints[ng];
      
      for(int j=0; j<nlayer; j++){
	if(w0[i][j]>=0.01){
	  double alphax = sqrt( (1.-w0[i][j])/(1.-w0[i][j]*cosbar[j]) );
	  gg[j] = 2.*pi*w0[i][j]*xk1[j]*(1.+cosbar[j]*alphax)/(1.+alphax);
	  hh[j] = 2.*pi*w0[i][j]*xk2[j]*(1.-cosbar[j]*alphax)/(1.+alphax);
	  xj[j] = 2.*pi*w0[i][j]*xk1[j]*(1.-cosbar[j]*alphax)/(1.+alphax);
	  xk[j] = 2.*pi*w0[i][j]*xk2[j]*(1.+cosbar[j]*alphax)/(1.+alphax);
	  alpha1[j] = 2.*pi*(b0[j] + bdiff[j]*(ubari*w0[i][j]*cosbar[j]/(1.-w0[i][j]*cosbar[j])));
	  alpha2[j] = 2.*pi*bdiff[j];
	  sigma1[j] = 2.*pi*(b0[j] - bdiff[j]*(ubari*w0[i][j]*cosbar[j]/(1.-w0[i][j]*cosbar[j])));
	  sigma2[j] = alpha2[j];
	}
	else{
	  gg[j] = 0.;
	  hh[j] = 0.;
	  xj[j] = 0.;
	  xk[j] = 0.;
	  alpha1[j] = 2.*pi*b0[j];
	  alpha2[j] = 2.*pi*bdiff[j];
	  sigma1[j] = alpha1[j];
	  sigma2[j] = alpha2[j];
	}
      }
      
      fpt[nlevel-1] = 2.*pi*(bsurf + bdiff[nlayer-1]*ugauss);
      //fmt[0] = 2.*pi*(1.-exp(-tautop/ugauss))*blackbodyL(tprof[0],wavelength);
      
      for(int j=0; j<nlayer; j++){
	em1[j] = exp(-lamda[j]*taulayer[i][j]);
	em2[j] = exp(-taulayer[i][j]/ugauss);
	em3[j] = em1[j]*em2[j];
	epp[j] = exp( min(lamda[j]*taulayer[i][j], 35.) );
	/*
	fmt[j+1] = fmt[j]*em2[j];
	fmt[j+1] += xj[j]/(lamda[j]*ugauss+1.)*(epp[j]-em2[j]);
	fmt[j+1] += xk[j]/(lamda[j]*ugauss-1.)*(em2[j]-em1[j]);
	fmt[j+1] += sigma1[j]*(1.-em2[j]);
	fmt[j+1] += sigma2[j]*(ugauss*em2[j]+taulayer[i][j]-ugauss);
	*/
      }
      
      for(int j=nlayer-1; j>=0; j--){
	fpt[j] = fpt[j+1]*em2[j];
	fpt[j] += gg[j]/(lamda[j]*ugauss-1.)*(epp[j]*em2[j]-1.);
	fpt[j] += hh[j]/(lamda[j]*ugauss+1.)*(1.-em3[j]);
	fpt[j] += alpha1[j]*(1.-em2[j]);
	fpt[j] += alpha2[j]*(ugauss - (taulayer[i][j]+ugauss)*em2[j]);
      }
      
      totalflux[i] += gaussWeights[ng]*fpt[0];
    } // end of ng for loop
  } // end of i for loop (wavens)
  
  return totalflux;
}

// computes the total spectral radiance of the planet in frequency space
vector<double> Planet::getFluxOneStream(vector<double> wavens)
{
  // 8 Gauss points
  vector<double> gaussPoints(8,0);
  vector<double> gaussWeights(8,0);
  
  // Cosines of angles to calculate the gauss points
  gaussPoints[0] = 0.0446339553;
  gaussPoints[1] = 0.1443662570;
  gaussPoints[2] = 0.2868247571;
  gaussPoints[3] = 0.4548133152;
  gaussPoints[4] = 0.6280678354;
  gaussPoints[5] = 0.7856915206;
  gaussPoints[6] = 0.9086763921;
  gaussPoints[7] = 0.9822200849;

  // Weights for the gauss points
  gaussWeights[0] = 0.0032951914;
  gaussWeights[1] = 0.0178429027;
  gaussWeights[2] = 0.0454393195;
  gaussWeights[3] = 0.0791995995;
  gaussWeights[4] = 0.1060473494;
  gaussWeights[5] = 0.1125057995;
  gaussWeights[6] = 0.0911190236;
  gaussWeights[7] = 0.0445508044;
  
  vector<double> totalflux(wavens.size(),0);
  vector<double> fracflux(wavens.size(),0);
  
  for(int j=0; j<8; j++){
    fracflux = getFluxes(wavens,gaussPoints[j],gaussWeights[j]);
    for(int i=0; i<wavens.size(); i++){
      totalflux[i] += 2.*pi*fracflux[i]*gaussWeights[j]*pi/2.; // erg/s/cm^2/Hz
      // extra factor of pi/2 from the mu substitution
    }
  }
  return totalflux;
}

// computes the spectral radiance of the planet emitted at a given angle
vector<double> Planet::getFluxes(vector<double> wavens, double cosmu, double delMu)
{
  vector<double> fluxes(wavens.size(),0);
  double sinmu = sqrt(1.-cosmu*cosmu);
  
  for(int i=0; i<wavens.size(); i++){
    double freq = wavens[i] * c;
    double wavelength = 1./wavens[i];
    double taucolumn = tauprof[i][nlayer-1];
    
    for(int j=nlayer-1; j>=0; j--){
      // Entire layer below cloud deck
      if(hprof[j]<hmin){
	fluxes[i] += 0.;
      }
      
      // Entire layer above cloud deck
      else if(hprof[j+1]>=hmin){
	double tmid = 0.5*(tprof[j+1] + tprof[j]);
	double freq = wavens[i] * c;
	double wavelength = 1./wavens[i];
	double I0 = blackbodyL(tmid,wavelength);             // erg/s/sr/cm^2/Hz
	double dflux = I0 * (1.-exp(-taulayer[i][j]/cosmu)); // source function
	if(j>0) dflux *= exp(-tauprof[i][j-1]/cosmu);        // absorption above layer
	dflux *= sinmu;                                      // multiply by area
	fluxes[i] += dflux;
      }
      
      // Layer overlaps cloud deck boundary
      else{
	double fvisible = (hprof[j]-hmin)/(hprof[j]-hprof[j+1]);
	double tcloud = fvisible*tprof[j+1] + (1.-fvisible)*tprof[j];
	double tmid = 0.5*fvisible*tprof[j+1] + (1.-0.5*fvisible)*tprof[j];
	double freq = wavens[i] * c;
	double wavelength = 1./wavens[i];
	
	double I0 = blackbodyL(tcloud,wavelength);           // erg/s/sr/cm^2/Hz
	double tauadjust = fvisible*taulayer[i][j];

	double dflux = I0;
	if(j>0) dflux *= exp(-(tauadjust+tauprof[i][j-1])/cosmu);
	else dflux *= exp(-tauadjust/cosmu);

	double I1 = blackbodyL(tmid,wavelength);
	double dflux1 = I1 * (1.-exp(-tauadjust/cosmu));      // source function
	if(j>0) dflux1 *= exp(-tauprof[i][j-1]/cosmu);        // absorption above layer
	dflux += dflux1;
	dflux *= sinmu;                                      // multiply by area
	fluxes[i] += dflux;
      }
    }
  }
  return fluxes;
}

// computes the flux ratio received at each wavelength IN TRANSIT
vector<double> Planet::transFlux(double rs, vector<double> wavens)
{
  //const int nthreads = max(atoi(getenv("OMP_NUM_THREADS")),4);
  //const int chunk = wavens.size()/nthreads;
  vector<double> fluxes(wavens.size(),0);
  /*
#pragma omp parallel for \
  default(shared) \
  private(i,j) \
  schedule(static,chunk) \
  num_threads(nthreads)
  */
  for(int i=0; i<wavens.size(); i++){
    double freq = wavens[i] * c;
    double wavelength = 1./wavens[i];
    double I0 = blackbodyL(tstar,wavelength);
    double dflux = pi*rp*rp;

    for(int j=nlayer-1; j>=0; j--){
      // Entire layer below cloud deck
      if(hprof[j]<hmin){
	dflux += (hprof[j]-hprof[j+1]) * 2*pi*(rp+hprof[j]);
      }
      
      // Entire layer above cloud deck
      else if(hprof[j+1]>=hmin){
	dflux += (1 - exp(-tauprof[i][j])) * (hprof[j]-hprof[j+1]) * 2*pi*(rp+hprof[j]);
      }
      
      // Layer overlaps cloud deck boundary
      else{
	double fvisible = (hprof[j]-hmin)/(hprof[j]-hprof[j+1]);
	dflux += (1 - exp(-tauprof[i][j])) * (hprof[j]-hprof[j+1]) * 2*pi*(rp+hprof[j]) * fvisible;
      }
    }

    fluxes[i] = dflux / (pi*rs*rs);
  }
  return fluxes;
}

// Computes a 90-layer T-P profile from the input profile
void Planet::getProfile(vector<double> tpprofile)
{
  tprof = tpprofile;
  prprof = vector<double>(nlevel);

  for(int i=0; i<nlevel; i++){
    prprof[i] = pow(10,minP + (maxP-minP)*i/(nlevel-1));
  }

  double dpdr;
  
  hprof.back() = 0.;
  
  for(int i=nlevel-2; i>=0; i--){
    dpdr = G*mp*mu/k/tprof[i]/(rp+hprof[i+1])/(rp+hprof[i+1]);
    hprof[i] = hprof[i+1] + log(prprof[i+1]/prprof[i])/dpdr;
  }
}
// end getProfile

// Computes the optical depth table from the opacity table.
// Length of tauprof is 1 less than tprof
// because it computes the optical depth between T-P points.
// Index 0 is at the top of the atmosphere.
void Planet::getTauProf(vector<double> wavens){
  tauprof = vector<vector<double> >(wavens.size(),vector<double>(nlayer,0));
  taulayer = vector<vector<double> >(wavens.size(),vector<double>(nlayer,0));
  w0 = vector<vector<double> >(wavens.size(),vector<double>(nlayer,0));
  double hazedepth = 0.0;
  for(int i=0; i<wavens.size(); i++){
    double wavel = 10000./wavens[i];
    
    for(int j=0; j<nlayer; j++){
      double dl = (hprof[j]-hprof[j+1]);
      double opac = opacprof[i][j];
      double pmid = sqrt(prprof[j]*prprof[j+1]);
      double tmid = 0.5*(tprof[j]+tprof[j+1]);
      
      double dtau = opac*pmid/k/tmid;
      if(j>0) tauprof[i][j] = tauprof[i][j-1] + dtau * dl;
      taulayer[i][j] = dtau * dl;
      
      double dlc = 0.; // thickness of cloud inside the layer
      double hazeabund = 0.;
      if(hazetype!=0){
	double hazexsec = atmosphere->getXsec(wavel,haze[1]);
	
	// Slab cloud model
	if(cloudmod==2){
	  hazedepth = hazexsec * haze[0];
	  // layer is strictly inside the cloud
	  if(prprof[j] >= haze[2] && prprof[j+1] <= haze[3]){
	    dlc = dl;
	  }
	  // layer overlaps bottom of cloud
	  else if(prprof[j] >= haze[2] && prprof[j] <= haze[3]){
	    dlc = dl * (log10(haze[3]/prprof[j]) / log10(haze[3]/haze[2]));
	  }
	  // layer overlaps top of cloud
	  else if(prprof[j+1] >= haze[2] && prprof[j+1] <= haze[3]){
	    dlc = dl * (log10(prprof[j+1]/haze[2]) / log10(haze[3]/haze[2]));
	  }
	  // layer contains entire cloud
	  else if(prprof[j] <= haze[2] && prprof[j+1] >= haze[3]){
	    dlc = dl * (log10(prprof[j+1]/haze[2]) / log10(haze[3]/haze[2]));
	  }
	}

	// Gaussian cloud model, all layers are nominally inside the cloud
	if(cloudmod==3){
	  hazeabund = haze[0]*gauss(log10(prprof[j]),haze[2],haze[3]);
	  hazedepth = hazexsec * hazeabund;
	  dlc = dl;
	}

      } // end if(hazetype!=0)
      else hazedepth = 0.;
	
      if(j==0){
	tauprof[i][j] = 0.;
      }
      else{
	tauprof[i][j] += hazedepth * dlc;
      }
      taulayer[i][j] += hazedepth * dlc;
      w0[i][j] = scatable[i][j]*pmid/k/tmid*0.999999/(dtau+hazedepth);
      // w0 is layer single scattering albedo with layer 0 on top
      // ratio of scattering to total opacity
    }
  }
  
  for(int j=0; j<nlayer; j++){
    double tross = 0.;
    double trosso = 0.;
    for(int i=0; i<wavens.size(); i++){
      tross += tauprof[i][j]/wavens.size();
    }
    if(tross>0.67){
      double dross = (tross-0.67)/(tross-trosso);
      teff = (1.-dross)*tprof[j] + dross*tprof[j-1];
      break;
    }
    else if(hprof[j]<hmin){
      double dross = (hmin-hprof[j])/(hprof[j+1]-hprof[j]);
      teff = (1.-dross)*tprof[j] + dross*tprof[j-1];
      break;
    }
    else{
      if(j==nlayer-1){
	teff = tprof[nlayer-1];
	break;
      }
      else{
	tross = trosso;
      }
    }
  }
}

// Computes the optical depth table IN TRANSIT from the opacity table
void Planet::transTauProf(vector<double> wavens){
  tauprof = vector<vector<double> >(wavens.size(),vector<double>(nlayer,0));
  double hazedepth = 0.0;

  vector<vector<double> > dlgrid(nlayer,vector<double>(nlayer,0));
  for(int i=1; i<nlayer; i++){
    for(int j=0; j<i; j++){
      double ltot = sqrt( (rp+hprof[j])*(rp+hprof[j]) - (rp+hprof[i])*(rp+hprof[i]) );
      double lnew = sqrt( (rp+hprof[j+1])*(rp+hprof[j+1]) - (rp+hprof[i])*(rp+hprof[i]) );
      dlgrid[i][j] = ltot - lnew;
    }
  }
  
  for(int ii=0; ii<wavens.size(); ii++){
    double wavel = 10000./wavens[ii];
    for(int i=1; i<nlayer; i++){
      for(int j=0; j<i; j++){
	double dtau = opacprof[ii][j] * prprof[j]/k/tprof[j];
	double dl = (hprof[j]-hprof[j+1]);
	double dlc = 0.; // thickness of cloud inside the layer
	tauprof[ii][i] += dtau * dlgrid[i][j] * (hprof[i]-hprof[i+1]);
	double hazeabund = 0.;
	
	if(hazetype!=0){
	  double hazexsec = atmosphere->getXsec(wavel,haze[1]);
	  
	  // Slab cloud model
	  if(cloudmod==2){
	    hazedepth = hazexsec * haze[0];
	    // layer is strictly inside the cloud
	    if(prprof[j] >= haze[2] && prprof[j+1] <= haze[3]){
	      dlc = dl;
	    }
	    // layer overlaps bottom of cloud
	    else if(prprof[j] >= haze[2] && prprof[j] <= haze[3]){
	      dlc = dl * (log10(haze[3]/prprof[j]) / log10(haze[3]/haze[2]));
	    }
	    // layer overlaps top of cloud
	    else if(prprof[j+1] >= haze[2] && prprof[j+1] <= haze[3]){
	      dlc = dl * (log10(prprof[j+1]/haze[2]) / log10(haze[3]/haze[2]));
	    }
	    // layer contains entire cloud
	    else if(prprof[j] <= haze[2] && prprof[j+1] >= haze[3]){
	      dlc = dl * (log10(prprof[j+1]/haze[2]) / log10(haze[3]/haze[2]));
	    }
	  }
	
	  // Gaussian cloud model, all layers are nominally inside the cloud
	  if(cloudmod==3){
	    hazeabund = haze[0]*gauss(log10(prprof[j]),haze[2],haze[3]);
	    hazedepth = hazexsec * hazeabund;
	    dlc = dl;
	  }
	    
	  tauprof[ii][i] += hazedepth * dlgrid[i][j] * dlc;
	  // dlgrid = path length, deltaH = layer height, but multiply by dlc/deltaH to account for haze layer thickness
	}
      }
      tauprof[ii][i] /= (hprof[0]-hprof[i]);
      tauprof[ii][i] *= 2.;
    }
  }
}

// Computes the opacity table for the computed T-P profile from the input opacities.
// This is also 1 shorter than the T-P profile.
void Planet::getOpacProf(double rxsec, vector<double> wavens, vector<double> abund){
  getSca(rxsec,wavens);
  opacprof = vector<vector<double> >(wavens.size(),vector<double>(nlayer,0));
  int degrade = round(10000.*log(wavens[0]/wavens[1]));
  int nmol = nspec;
  if(nmol==0) nmol = 1;
  
  double wmin = 10000./(0.6*exp(2.1204));
  double wmax = 10000./0.6;
  int numwave = wavens.size();
  double pmin = 0.0;            // original opacity files in mbar
  double pmax = 8.5;            // converted to cgs here
  int numpr   = 18;
  double tmin = log10(75.);
  double tmax = log10(75.)+34./20.3;
  int numtemp = 35;
  
  double wval = log(10000./wavens[0]/0.6)*10000./degrade;
  int wstart = (int)wval;

  // interpolation variables
  double opr1, opr2, opac;
  vector<double> xsec(4,0);

  // Compute the indexes for interpolation
  vector<double> logtemp(numtemp,0);
  for(int n=0; n<numtemp; n++){
    logtemp[n] = log10(75.*pow(10,n/20.3));
  }
  vector<double> logpr(numpr,0);
  for(int n=0; n<numpr; n++){
    logpr[n] = n*0.5;
  }
  // loop over wavelengths
  for(int m=0; m<numwave; m++){
    // loop over T-P profile
    for(int j=0; j<nlayer; j++){
      double tl = log10(0.5*(tprof[j]+tprof[j+1]));
      double deltat = (tl-tmin)/(tmax-tmin)*(numtemp-1.);
      int jt = (int)deltat;
      double dti = (tl-logtemp[jt])/(logtemp[jt+1]-logtemp[jt]);
      if(jt<0){
	jt=0.;
	dti=0.;
      }
      if(jt>numtemp-2){
	jt=numtemp-2;
	dti=0.999999;
      }
      double pl = 0.5*(log10(prprof[j]*prprof[j+1]));
      double deltap = (pl-pmin)/(pmax-pmin)*(numpr-1.);
      int jp = (int)deltap;
      double dpi = (pl-logpr[jp])/(logpr[jp+1]-logpr[jp]);
      if(jp<0){
	jp=0.;
	dpi=0.;
      }
      if(jp>numpr-2){
	jp=numpr-2;
	dpi=0.999999;
      }

      xsec[0] = 0.;
      xsec[1] = 0.;
      xsec[2] = 0.;
      xsec[3] = 0.;
      
      for(int iii=0; iii<nmol; iii++){
	xsec[0] += mastertable[jp][jt][m+wstart][iii]*abund[iii];
	xsec[1] += mastertable[jp][jt+1][m+wstart][iii]*abund[iii];
	xsec[2] += mastertable[jp+1][jt][m+wstart][iii]*abund[iii];
	xsec[3] += mastertable[jp+1][jt+1][m+wstart][iii]*abund[iii];
      }
      // interpolate opacity
      opr1 = xsec[0] + dpi*(xsec[2]-xsec[0]);
      opr2 = xsec[1] + dpi*(xsec[3]-xsec[1]);
      opac = opr1 + dti*(opr2-opr1);

      opacprof[m][j] = opac + scatable[m][j];
    }
  }
  return;
}
// end getOpacProf

// Computes the scattering opacity table
void Planet::getSca(double rxsec, vector<double> wavens){
  scatable = vector<vector<double> >(wavens.size(),vector<double>(nlayer,0));
  for(int m=0; m<wavens.size(); m++){
    for(int j=0; j<nlayer; j++){
      double nu = c*wavens[m];  
      scatable[m][j] = rxsec * pow(nu/5.0872638e14,4.0);
    }
  }
}
// end getSca