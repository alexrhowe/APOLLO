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
  int hazetype = 0;
  vector<double> haze(4,0);
  vector<double> abund(7,0);
  double wave, rad, realn, imagn, xsec;
  std::string opacdir = "../Opacities";
  std::string fname = "blank.r.dat";
  std::string fOut = "haze.table.dat";
  if(argc>1){
    opacdir = argv[1];
  }  
  if(argc>2){
    hazetype = 4;
    fname = opacdir + "/" + argv[2];
  }
  if(argc>3){
    fOut = opacdir + "/" + argv[3];
  }
  
  haze[0] = 1.0; haze[1] = 1.0; haze[2] = 0.1; haze[3] = 1.0;
  abund[0] = 1.0; abund[1] = 0.0; abund[2] = 0.0; abund[3] = 0.0; abund[4] = 0.0; abund[5] = 0.0; abund[6] = 0.0;
  
  Atmosphere atmosphere = Atmosphere(hazetype,haze,fname);

  std::ofstream output;
  output.open(fOut.c_str());

  double wmin = 0.6;
  double wmax = 30.0;
  double resolv = 10000.;
  int nwave = (int)(ceil(log(wmax/wmin)*resolv))+1;
  
  int nsize = 161;
  double smin = -4.0; // particle radii in microns
  double smax = 4.0;
  
  output << nwave << " " << wmin << " " << wmax << " " << nsize << " " << smin << " " << smax << "\n";
  
  for(int i=0; i<nwave; i++){
    wave = wmin*exp(i/resolv);
    output << std::setprecision(10) << std::showpoint << wave;

    for(int j=0; j<161; j++){
      rad = 0.0001*pow(10,j/20.);
  
      realn = atmosphere.value(wave,atmosphere.einn,atmosphere.elam);
      imagn = atmosphere.value(wave,atmosphere.eink,atmosphere.elam);
      atmosphere.MieLargeSzpara(rad,wave,realn,imagn);
  
      xsec = atmosphere.qe*rad*rad*pi/4.e8;
      if(isnan(xsec)) xsec = 0.;
      
      output << " " << xsec;
    }
    output << "\n";
    if(i%100==0) printf("%d\n",i);
  }
}
