#ifndef PLANET_LAYER_H
#define PLANET_LAYER_H

#include "specincld.h"
#include "Atmosphere.h"

class Planet
{
 private:
  double Tbase;  // temperature at the base of the T-P profile
  
  double mp;     // planet mass
  double mu;     // mean molecular weight
  //double g;      // surface gravity
  double tstar;  // stellar temperature

  // implement as vectors
  vector<double> prprof;
  vector<double> tprof;
  vector<double> wavens;
  vector<vector<double > > scatable;
  vector<vector<double > > h2prof;
  vector<vector<double > > nakprof;

  int dsteps;    // integration steps in horizontal direction

 public:
  Planet(int modenum, vector<double> waves, int cloudnum, int hazenum, int wrange, vector<int> mollist, string opacname);

  int nlayer;
  int nlevel;
  int npress;
  int ntemp;
  int nwave;
  int nspec;
  int mode;
  int hazetype;
  int ntable;
  int degrade;
  int cloudmod;
  int range;
  double rp;     // planet radius
  double rs;
  double rxsec;
  double minP;
  double maxP;
  double sma;
  double grav;   // log(g)
  double hmin;   // minimum height treated in atmosphere
  double hmax;   // maximum height treated in atmosphere
  double wavel; // wavelength for optical depth calculation (in microns)
  double hazeopac;
  int bsteps;    // integration steps in vertial direction
  double lmin;
  double lmax;
  double res;
  vector<double> haze;
  vector<double> tpprof;
  vector<double> hprof;
  vector<vector<double > > opacprof;
  vector<vector<double > > tauprof;
  vector<vector<double > > taulayer;
  vector<vector<double > > w0;
  vector<vector<vector<double> > > opacities;
  vector<vector<vector<vector<double> > > > mastertable;
  string opacdir;
  Atmosphere* atmosphere;
  
  void setParams(vector<double> plparams, vector<double> abund, vector<double> tpprofile);
  vector<double> getSpectrum(int streams);
  void readopac(vector<int> mollist, vector<double> wavens, string opacdir);
  void setWave(int npoints, double rxsec, vector<double> wavens, vector<double> abund);
  double getP(double height);
  double getT(double height);
  double getH(double pr);
  vector<double> getFlux(vector<double> wavens);
  vector<double> getFluxOneStream(vector<double> wavens);
  vector<double> transFlux(double rs, vector<double> wavens);
  vector<double> getFluxes(vector<double> wavens, double cosmu, double delMu);
  //double chord(double b, int wave);
  //double Aefftrans(double b, int wave, double l, double z);
  void getTauProf(vector<double> wavens);
  void transTauProf(vector<double> wavens);
  void getOpacProf(double rxsec, vector<double> wavens, vector<double> abund);
  void getProfile(vector<double> tpprofile);
  void getSca(double rxsec, vector<double> wavens);
  double HminBoundFree(double t, double waven);
  double HminFreeFree(double t, double waven);
};

#endif // FASTPLANET_H
