#ifndef PLANET_H
#define PLANET_H

#include "specincld.h"
#include "Atmosphere.h"

class Planet
{
 private:
  double Tbase;  // temperature at the base of the T-P profile
  
  double mp;     // planet mass
  double mu;     // mean molecular weight
  double tstar;  // stellar temperature

  // implement as vectors
  vector<double> prprof;
  vector<double> rho;
  vector<double> taulist;
  vector<double> tprof;
  vector<double> wavens;
  vector<double> wavenslo;
  vector<vector<double > > scatable;
  vector<vector<double > > scatablelo;

 public:
  Planet(vector<int> switches, vector<double> waves, vector<double> waveslo, vector<int> mollist, string opacname, string hir, string lor);

  int nlayer;
  int nlevel;
  int npress;
  int ntemp;
  int nwave;
  int nwavelo;
  int nspec;
  int mode;
  int hazetype;
  int streams;
  int ntable;
  int degrade;
  int cloudmod;
  int tprofmode;
  double rp;     // planet radius
  double rs;
  double rxsec;
  double minP;
  double maxP;
  double sma;
  double grav;   // log(g)
  double hmin;   // minimum height treated in atmosphere
  double hmax;   // maximum height treated in atmosphere
  double wavel;  // wavelength for optical depth calculation (in microns)
  double pmin;
  double pmax;
  double tmin;
  double tmax;
  double lmin;
  double lmax;
  double wmin;
  double wmax;
  double res;
  double lminlo;
  double lmaxlo;
  double reslo;
  vector<double> haze;
  vector<double> tpprof;
  vector<double> hprof;
  vector<vector<double > > opacprof;
  vector<vector<double > > tauprof;
  vector<vector<double > > taulayer;
  vector<vector<double > > w0;
  vector<vector<double > > asym;
  vector<vector<double > > opacproflo;
  vector<vector<double > > tauproflo;
  vector<vector<double > > taulayerlo;
  vector<vector<double > > w0lo;
  vector<vector<double > > asymlo;
  vector<vector<vector<double> > > opacities;
  vector<vector<vector<vector<double> > > > mastertable;
  vector<vector<vector<vector<double> > > > lotable;
  string opacdir;
  bool doMie;
  string hires;
  string lores;
  Atmosphere* atmoshires;
  Atmosphere* atmoslores;
  
  void setParams(vector<double> plparams, vector<double> abund, vector<double> tpprofile);
  double getTeff();
  vector<double> getSpectrum();
  void readopac(vector<int> mollist, vector<double> wavens, string table, string opacdir);
  void setWave(int npoints, double rxsec, vector<double> wavelist, vector<double> abund);
  double getP(double height);
  double getT(double height);
  double getH(double pr);
  vector<double> getFlux(vector<double> wavens, string table);
  vector<double> getFluxOneStream(vector<double> wavens, string table);
  vector<double> getFluxes(vector<double> wavens, double cosmu, double delMu, string table);
  vector<double> transFlux(double rs, vector<double> wavens, string table);
  void getProfLayer(vector<double> tpprofile);
  void getProfParam(vector<double> tpprofile);
  void getTauProf(vector<double> wavens, string table);
  void transTauProf(vector<double> wavens, string table);
  void getOpacProf(double rxsec, vector<double> wavelist, vector<double> abund, string table);
  void getSca(double rxsec, vector<double> wavelist, string table);
  double HminBoundFree(double t, double waven);
  double HminFreeFree(double t, double waven);
};

#endif // PLANET_H
