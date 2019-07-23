#ifndef ATMOSPHERE_H
#define ATMOSPHERE_H

#include "specincld.h"

class Atmosphere
{
 private:
  double s1;
  double s2;
  double pext;
  double qext;
  double s1r;
  double s1i;
  double s2r;
  double s2i;
  double* ar;
  double* ai;
  double* br;
  double* bi;
  double* rsr;
  double* rsi;
  double* rsx;
  double* px;
  int hazeNumVals;

  int smi1xs(double x, int nx, double xm, double ym);
  void sm5msx(double x, int i);
  void smi3sm(int i, double thd);
  void sm5pqs(double x, int nm);

 public:
  //Atmosphere(int metalint, double* mols, int hazeint, double* hazeparams);
  Atmosphere(int hazeint, vector<double> hazeparams);
  double qe;
  double qa;
  double qs;
  double qp;
  double asf;
  int metals;
  int hazetype;
  vector<double> molarray;
  vector<double> haze;
  double* einn;
  double* eink;
  double* elam;

  int nwave;
  double wmin;
  double wmax;
  int nsize;
  double smin;
  double smax;
  vector<double> logsize;
  vector<double> waves;
  vector<vector<double> > hazetab;

  double getXsec(double wave, double size);
  double value(double x, double* sigma, double* e);
  void MieLargeSzpara(double size, double wavelength, double xm, double ymm);


};

#endif // ATMOSPHERE_H
