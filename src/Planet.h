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
  // ada: Adding species-specific Rayleigh scattering arrays.
  vector<vector<vector<double > > > specscatable;
  vector<vector<vector<double > > > specscatablelo;

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
  double deltalogt;
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
  vector<vector<vector<double > > > specopacprof;
  vector<vector<double > > tauprof;
  vector<vector<double > > taulayer;
  // ada: Adding the individual optical depths from the gas and the cloud layers.
  vector<vector<double > > cloudtauprof;
  vector<vector<double > > cloudtaulayer;
  vector<vector<double > > gastauprof;
  vector<vector<double > > gastaulayer;
  // ada: Adding the blackbody spectrum per layer as a public attribute for use in calculating the contribution function. (blayer[i][j] -> b0[j])
  vector<vector<double > > blayer;
  vector<vector<double > > w0;
  vector<vector<double > > gasw0;
  vector<vector<double > > asym;
  vector<vector<double > > gasasym;
  vector<vector<double > > opacproflo;
  vector<vector<vector<double > > > specopacproflo;
  vector<vector<double > > tauproflo;
  vector<vector<double > > taulayerlo;
  vector<vector<double > > gastauproflo;
  vector<vector<double > > gastaulayerlo;
  vector<vector<double > > w0lo;
  vector<vector<double > > gasw0lo;
  vector<vector<double > > asymlo;
  vector<vector<double > > gasasymlo;
  vector<vector<vector<double> > > opacities;
  vector<vector<vector<vector<double> > > > mastertable;
  vector<vector<vector<vector<double> > > > lotable;
  string opacdir;
  bool doMie;
  string hires;
  string lores;
  Atmosphere* atmoshires;
  Atmosphere* atmoslores;
  
  void setParams(vector<double> plparams, vector<double> abund, vector<double>rxsecs, vector<double> tpprofile);
  //void setParams(vector<double> plparams, vector<vector<double> > abund, vector<double>rxsecs, vector<double> tpprofile);
  double getTeff();
  vector<double> getSpectrum();
  vector<double> getClearSpectrum();
  void readopac(vector<int> mollist, vector<double> wavens, string table, string opacdir);
  void setWave(int npoints, vector<double> rxsec, vector<double> wavelist, vector<double> abund);
  //void setWave(int npoints, vector<double> rxsec, vector<double> wavelist, vector<vector<double> > abund);
  vector<vector<double> > getContribution();
  vector<vector<double> > getCloudContribution();
  vector<vector<double> > getGasContribution();
  vector<vector<vector<double> > > getSpeciesContribution();
  double getP(double height);
  double getT(double height);
  double getH(double pr);
  vector<double> getFlux(vector<double> wavens, string table, vector<vector<double > > taulayer, vector<vector<double > > w0, vector<vector<double > > asym);
  vector<double> getFluxOneStream(vector<double> wavens, string table, vector<vector<double > > taulayer, vector<vector<double > > tauprof);
  vector<double> getFluxes(vector<double> wavens, double cosmu, double delMu, string table, vector<vector<double > > taulayer, vector<vector<double > > tauprof);
  vector<double> transFlux(double rs, vector<double> wavens, string table);
  void getProfLayer(vector<double> tpprofile);
  void getProfParam(vector<double> tpprofile);
  void getTauProf(vector<double> wavens, string table);
  void transTauProf(vector<double> wavens, string table);
  void getOpacProf(vector<double> rxsecs, vector<double> wavelist, vector<double> abund, string table);
  void getOpacProf(vector<double> rxsecs, vector<double> wavelist, vector<vector<double> > abund, string table);
  void getSca(vector<double> rxsecs, vector<double> wavelist, string table);
  double HminBoundFree(double t, double waven);
  double HminFreeFree(double t, double waven);
};

#endif // PLANET_H
