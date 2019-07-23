/*
 * List of constants for the Conjunction Code
 * Alex R. Howe
 */
#ifndef CONJ_CONST_H
#define CONJ_CONST_H

// all in cgs units, so don't screw it up
// except wavelengths all in microns
namespace cons
{
  extern const double G;
  extern const double c;
  extern const double pi;
  extern const double h;
  extern const double hc;
  extern const double k;
  extern const double N_A;
  extern const double M_EARTH;
  extern const double R_EARTH;

  double integrator();
  double gauss(double x, double m, double s);
  double blackbodyN(double T, double freq);
  double blackbodyL(double T, double wave);
  double expint(int n, double x);
}
#endif
