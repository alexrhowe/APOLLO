// Computes the optical depth table IN TRANSIT from the opacity table
void Planet::transTauProf(vector<double> wavens, string table)
{
  if(table=="hires"){
    tauprof = vector<vector<double> >(wavens.size(),vector<double>(nlayer,0)); // Would be more accurately be called taulayer.
    double hazedepth = 0.0;
    doMie = false;

    vector<vector<double> > dlgrid(nlayer,vector<double>(nlayer,0));
    for(int i=1; i<nlayer; i++){
      for(int j=0; j<i; j++){
	double ltot = sqrt( (rp+hprof[j])*(rp+hprof[j]) - (rp+hprof[i])*(rp+hprof[i]) );
	double lnew = sqrt( (rp+hprof[j+1])*(rp+hprof[j+1]) - (rp+hprof[i])*(rp+hprof[i]) );
	dlgrid[i][j] = (ltot - lnew) / (hprof[j]-hprof[j+1]);
      }
    }
  
    for(int ii=0; ii<wavens.size(); ii++){
      double wavel = wavens[ii];
      for(int i=1; i<nlayer; i++){
	for(int j=0; j<i; j++){
	  double dl = (hprof[j]-hprof[j+1]);
	  double opac = opacprof[ii][j];
	  double pmid = sqrt(prprof[j]*prprof[j+1]);
	  double tmid = 0.5*(tprof[j]+tprof[j+1]);

	  double dtau = opac*pmid/k/tmid;
	  
	  tauprof[ii][i] += dtau * dlgrid[i][j] * dl;
	  
	  double dlc = 0.; // thickness of cloud inside the layer
	  double hazeabund = 0.;
	  
	  if((hazetype!=0 && cloudmod>1) || cloudmod==4){
	    double hazexsec;

	    if(cloudmod!=4){
	      double absorbxsec = atmoshires->getAbsXsec(wavel,haze[1]);
	      double scatterxsec = atmoshires->getScaXsec(wavel,haze[1]);
	      hazexsec = absorbxsec + scatterxsec;
	    }
	    // ada: For the power-law opacity cloud model, we use the hazexsec variable as the wavelength scaling (haze[1] is the power-law exponent). haze[0] is a constant optical depth per unit length, or linear attenuation coefficient, for the cloud.
	    if(cloudmod==4){
	      hazexsec = pow(wavel, haze[0]);
	    }
	    if(mode<=1 && streams==2 && cloudmod!=4){
	      doMie = true;
	    }
	    
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
	    if(cloudmod==4){
	      // ada: Since the optical depth profile is set directly in this model, rather than the haze density, we don't need the path length.
	      //hazeabund = 1.;
	      //double scale = (haze[1]*(haze[2]-1))/haze[2];
	      //dlc = exp((prprof[j+1]-haze[1])/scale) - exp((prprof[j]-haze[1])/scale);
	      //if(dlc > 100.){
	      //  dlc = 100.;
	      //}
	      hazeabund = haze[3];
	      double P_base = haze[1] * haze[2];
	      dlc = (prprof[j+1]*prprof[j+1]/(prprof[j]*prprof[j])) / ((P_base*P_base)-(haze[1]*haze[1]));
	      hazedepth = hazexsec * hazeabund;
	    }
	    
	    tauprof[ii][i] += hazedepth * dlgrid[i][j] * dlc;
	    // dlgrid = path length, deltaH = layer height, but multiply by dlc/deltaH to account for haze layer thickness
	  }
	}
	tauprof[ii][i] *= 2.; // For the sunward and antisunward halves of the limb of the planet.
      }
    }
  }
  
  if(table=="lores"){
    tauproflo = vector<vector<double> >(wavenslo.size(),vector<double>(nlayer,0)); // Would be more accurately be called taulayer.
    double hazedepth = 0.0;
    doMie = false;

    vector<vector<double> > dlgrid(nlayer,vector<double>(nlayer,0));
    for(int i=1; i<nlayer; i++){
      for(int j=0; j<i; j++){
	double ltot = sqrt( (rp+hprof[j])*(rp+hprof[j]) - (rp+hprof[i])*(rp+hprof[i]) );
	double lnew = sqrt( (rp+hprof[j+1])*(rp+hprof[j+1]) - (rp+hprof[i])*(rp+hprof[i]) );
	dlgrid[i][j] = (ltot - lnew) / (hprof[j]-hprof[j+1]);
      }
    }
  
    for(int ii=0; ii<wavenslo.size(); ii++){
      double wavel = wavenslo[ii];
      for(int i=1; i<nlayer; i++){
	for(int j=0; j<i; j++){
	  double dl = (hprof[j]-hprof[j+1]);
	  double opac = opacprof[ii][j];
	  double pmid = sqrt(prprof[j]*prprof[j+1]);
	  double tmid = 0.5*(tprof[j]+tprof[j+1]);

	  double dtau = opac*pmid/k/tmid;
	  
	  tauproflo[ii][i] += dtau * dlgrid[i][j] * dl;
	  
	  double dlc = 0.; // thickness of cloud inside the layer
	  double hazeabund = 0.;
	  
	  if((hazetype!=0 && cloudmod>1) || cloudmod==4){
	    double hazexsec;

	    if(cloudmod!=4){
	      double absorbxsec = atmoshires->getAbsXsec(wavel,haze[1]);
	      double scatterxsec = atmoshires->getScaXsec(wavel,haze[1]);
	      hazexsec = absorbxsec + scatterxsec;
	    }
	    // ada: For the power-law opacity cloud model, we use the hazexsec variable as the wavelength scaling (haze[1] is the power-law exponent). haze[0] is a constant optical depth per unit length, or linear attenuation coefficient, for the cloud.
	    if(cloudmod==4){
	      hazexsec = pow(wavel, haze[0]);
	    }
	    if(mode<=1 && streams==2 && cloudmod!=4){
	      doMie = true;
	    }
	    
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
	    if(cloudmod==4){
	      // ada: Since the optical depth profile is set directly in this model, rather than the haze density, we don't need the path length.
	      dlc = 1.;
	      double scale = (haze[1]*(haze[2]-1))/haze[2];
	      hazeabund = exp((prprof[j+1]-haze[1])/scale) - exp((prprof[j]-haze[1])/scale);
	      if(hazeabund > 100.){
		hazeabund = 100.;
	      }
	      hazedepth = hazexsec * hazeabund;
	    }
	    
	    tauproflo[ii][i] += hazedepth * dlgrid[i][j] * dlc;
	    // dlgrid = path length, deltaH = layer height, but multiply by dlc/deltaH to account for haze layer thickness
	  }
	}
	tauproflo[ii][i] *= 2.; // For the sunward and antisunward halves of the limb of the planet.
      }
    }
  }
}
// end transTauProf
