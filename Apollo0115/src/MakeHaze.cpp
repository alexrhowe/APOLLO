#include "specincld.h"
#include "constants.h"
#include "Atmosphere.haze.h"
#include <iomanip>
#include <limits>
//#include "Planet.h"
//#include <valgrind/memcheck.h>

/* Required files:
   MakeHaze.cpp
   Atmosphere.haze.h
   Atmosphere.haze.cpp
   constants.h
   constants.cpp
   specincld.h
   Opacities directory
 */

using namespace cons;
using std::showpoint;

int main(int argc, char *argv[]){
  int hazetype = 4;
  vector<double> haze(4,0);
  vector<double> abund(7,0);
  double wave, rad, realn, imagn, xsec, axsec, sxsec, asfac;
  
  std::string opacdir = "../Opacities";
  std::string hazename = "corundum";
  std::string tabletype = "lores";

  double wmin = 0.6;
  double wmax = 30.0;
  double resolv = 200.;
  int nwave = (int)(ceil(log(wmax/wmin)*resolv))+1;
  
  if(argc>1){
    opacdir = argv[1];
  }
  if(argc>2){
    hazename = argv[2];
  }
  
  std::string fname = opacdir + "/indices/" + hazename + ".indices.dat";
  std::string faOut = opacdir + "/absorption/" + hazename + ".abs." + tabletype + ".dat";
  std::string fsOut = opacdir + "/scattering/" + hazename + ".sca." + tabletype + ".dat";
  std::string asfOut = opacdir + "/asymmetry/" + hazename + ".asf." + tabletype + ".dat";
  
  haze[0] = 1.0; haze[1] = 1.0; haze[2] = 0.1; haze[3] = 1.0;
  abund[0] = 1.0; abund[1] = 0.0; abund[2] = 0.0; abund[3] = 0.0; abund[4] = 0.0; abund[5] = 0.0; abund[6] = 0.0;
  
  Atmosphere atmosphere = Atmosphere(hazetype,haze,fname);

  std::ofstream aoutput;
  aoutput.open(faOut.c_str());
  std::ofstream soutput;
  soutput.open(fsOut.c_str());
  std::ofstream asfoutput;
  asfoutput.open(asfOut.c_str());
  
  int nsize = 161;
  double smin = -4.0; // particle radii in microns
  double smax = 4.0;
  
  aoutput << nwave << " " << wmin << " " << wmax << " " << nsize << " " << smin << " " << smax << "\n";
  soutput << nwave << " " << wmin << " " << wmax << " " << nsize << " " << smin << " " << smax << "\n";
  asfoutput << nwave << " " << wmin << " " << wmax << " " << nsize << " " << smin << " " << smax << "\n";
  
  for(int i=0; i<nwave; i++){
    wave = wmin*exp(i/resolv);
    aoutput << std::setprecision(10) << std::showpoint << wave;
    soutput << std::setprecision(10) << std::showpoint << wave;
    asfoutput << std::setprecision(10) << std::showpoint << wave;

    for(int j=0; j<161; j++){
      rad = 0.0001*pow(10,j/20.);
  
      realn = atmosphere.value(wave,atmosphere.einn,atmosphere.elam);
      imagn = atmosphere.value(wave,atmosphere.eink,atmosphere.elam);
      atmosphere.MieLargeSzpara(rad,wave,realn,imagn);
  
      xsec = atmosphere.qe*rad*rad*pi/4.e8;
      axsec = atmosphere.qa*rad*rad*pi/4.e8;
      sxsec = atmosphere.qs*rad*rad*pi/4.e8;
      //asfac = atmosphere.asf*rad*rad*pi/4.e8;
      asfac = atmosphere.asf;
      if(isnan(xsec)) xsec = 0.;
      if(isnan(axsec)) axsec = 0.;
      if(isnan(sxsec)) sxsec = 0.;
      if(isnan(asfac)) asfac = 0.;
      
      aoutput << " " << axsec;
      soutput << " " << sxsec;
      asfoutput << " " << asfac;
    }
    aoutput << "\n";
    soutput << "\n";
    asfoutput << "\n";
    if(i%100==0) printf("%d\n",i);
  }
}
