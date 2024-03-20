void Planet::getTauProf(vector<double> wavens, string table)
{
  if(table=="hires"){
    tauprof = vector<vector<double> >(wavens.size(),vector<double>(nlayer,0));
    taulayer = vector<vector<double> >(wavens.size(),vector<double>(nlayer,0));
    w0 = vector<vector<double> >(wavens.size(),vector<double>(nlayer,0));
    gasw0 = vector<vector<double> >(wavens.size(),vector<double>(nlayer,0));
    asym = vector<vector<double> >(wavens.size(),vector<double>(nlayer,0));
    gasasym = vector<vector<double> >(wavens.size(),vector<double>(nlayer,0));

    // ada: Adding the individual optical depths from the gas and the cloud layers.
    cloudtauprof = vector<vector<double> >(wavens.size(),vector<double>(nlayer,0));
    cloudtaulayer = vector<vector<double> >(wavens.size(),vector<double>(nlayer,0));
    gastauprof = vector<vector<double> >(wavens.size(),vector<double>(nlayer,0));
    gastaulayer = vector<vector<double> >(wavens.size(),vector<double>(nlayer,0));
    
    double forwardfrac = 0.;
    double backwardfrac = 0.;
    
    for(int i=0; i<wavens.size(); i++){
      double wavel = wavens[i];
      
      for(int j=0; j<nlayer; j++){
	double dl = (hprof[j]-hprof[j+1]);
	double opac = opacprof[i][j];
	double sca  = scatable[i][j];
	double abs  = opac - sca;
	double pmid = sqrt(prprof[j]*prprof[j+1]);
	double tmid = 0.5*(tprof[j]+tprof[j+1]);
	double nden = pmid/k/tmid; // Number density of gas molecules.
	
	double dtau = opac * nden;
	if(j>0){
	  gastauprof[i][j] = gastauprof[i][j-1] + dtau * dl;
	  tauprof[i][j] = tauprof[i][j-1] + dtau * dl;
	}
	gastaulayer[i][j] = dtau * dl;
	taulayer[i][j] = dtau * dl;

	gasw0[i][j] = sca / opac;
	
	double dlc = 0.; // thickness of cloud inside the layer
	double hazeabund = 0.;
	double hazedepth = 0.;
	double absdepth = 0.;
	double scadepth = 0.;

	double absorbxsec = 0.;
	double scatterxsec = 0.;
	double hazexsec = 0.;
	
	doMie = false;
	
	if((hazetype!=0 && cloudmod>1) || cloudmod==4){
	  if(cloudmod!=4){
	    absorbxsec = atmoshires->getAbsXsec(wavel,haze[1]);
	    scatterxsec = atmoshires->getScaXsec(wavel,haze[1]);
	    hazexsec = absorbxsec + scatterxsec;
	    //if(j==0 && wavel<1.20){
	    // printf("The hazexsec at wavel %f is %.10e.\n", wavel, hazexsec);
	    //}

	    if(mode<=1 && streams==2){
	      doMie = true;
	    }
	  }
	  // ada: For the power-law opacity cloud model, we use the hazexsec variable as the wavelength scaling (haze[1] is the power-law exponent). haze[0] is a constant optical depth per unit length, or linear attenuation coefficient, for the cloud.
	  else{
	    hazexsec = pow(wavel, haze[0]);
	  }
	  
	  // Slab cloud model
	  if(cloudmod==2){
	    hazeabund = haze[0];
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
	    double P_top = haze[1];
	    double P_base = P_top * haze[2];
	    hazeabund = haze[3] * ((prprof[j+1]*prprof[j+1])-(prprof[j]*prprof[j])) / ((P_base*P_base)-(P_top*P_top));
	    //if (i == wavens.size()-1) {
	    //  printf("dl in layer %d (logp is %f) is %f.\n", j, log10(prprof[j]), dl);
	    //}
	    // layer is strictly inside the cloud
	    if(prprof[j] >= P_top && prprof[j+1] <= P_base){
	      dlc = 1.;
	    }
	    // layer overlaps bottom of cloud
	    else if(prprof[j] >= P_top && prprof[j] <= P_base){
	      dlc = log10(P_base/prprof[j]) / log10(P_base/P_top);
	    }
	    // layer overlaps top of cloud
	    else if(prprof[j+1] >= P_top && prprof[j+1] <= P_base){
	      dlc = log10(prprof[j+1]/P_top) / log10(P_base/P_top);
	    }
	    // layer contains entire cloud
	    else if(prprof[j] <= P_top && prprof[j+1] >= P_base){
	      dlc = log10(prprof[j+1]/P_top) / log10(P_base/P_top);
	    }
	  }
	  
	  hazedepth = hazexsec * hazeabund;
	  absdepth = absorbxsec * hazeabund;
	  scadepth = scatterxsec * hazeabund;
	  
	} // end if(hazetype!=0)
	
	if(j==0){
	  cloudtauprof[i][j] = 0.;
	  tauprof[i][j] = 0.;
	}
	else{
	  cloudtauprof[i][j] += hazedepth * dlc;
	  tauprof[i][j] += hazedepth * dlc;
	}
	cloudtaulayer[i][j] += hazedepth * dlc;
	taulayer[i][j] += hazedepth * dlc;
        
	if(doMie==true){
	  // w0 is the single scattering albedo: the ratio of the scattering optical depth to the total optical depth
	  w0[i][j] = (sca*nden + scadepth) / (dtau + hazedepth);
	  
	  if(j==0){
	    backwardfrac = atmoshires->getAsym(wavel,haze[1]);
	    forwardfrac = 1.-backwardfrac;
	  }
	  // Hemispheric approximation to the asymmetry parameter integral.
	  asym[i][j] = (forwardfrac-backwardfrac)*scatterxsec / (sca + (forwardfrac+backwardfrac)*scatterxsec);
	}
	// ada: For the power-law opacity cloud model, the single-scattering albedo is a free parameter (constant in wavelength).
	else if(cloudmod==4){
	  w0[i][j] = (sca*nden + haze[4]*hazedepth) / (dtau + hazedepth);
	}
	else w0[i][j] = sca*nden / (dtau + hazedepth);
	// w0 is layer single scattering albedo with layer 0 on top
	// ratio of scattering to total opacity
      } // end for(j)
    } // end for(i)
  } // end if(table=="hires")
  
  if(table=="lores"){
    tauproflo = vector<vector<double> >(wavenslo.size(),vector<double>(nlayer,0));
    taulayerlo = vector<vector<double> >(wavenslo.size(),vector<double>(nlayer,0));
    w0lo = vector<vector<double> >(wavenslo.size(),vector<double>(nlayer,0));
    gasw0lo = vector<vector<double> >(wavenslo.size(),vector<double>(nlayer,0));
    asymlo = vector<vector<double> >(wavenslo.size(),vector<double>(nlayer,0));
    gasasymlo = vector<vector<double> >(wavenslo.size(),vector<double>(nlayer,0));

    // ada: Adding the individual optical depths from the gas.
    gastauproflo = vector<vector<double> >(wavens.size(),vector<double>(nlayer,0));
    gastaulayerlo = vector<vector<double> >(wavens.size(),vector<double>(nlayer,0));
    
    double forwardfrac = 0.;
    double backwardfrac = 0.;
    
    for(int i=0; i<wavenslo.size(); i++){
      double wavel = wavenslo[i];
      
      for(int j=0; j<nlayer; j++){
	double dl = (hprof[j]-hprof[j+1]);
	double opac = opacproflo[i][j];
	double sca  = scatablelo[i][j];
	double abs  = opac - sca;
	double pmid = sqrt(prprof[j]*prprof[j+1]);
	double tmid = 0.5*(tprof[j]+tprof[j+1]);
	double nden = pmid/k/tmid; // Number density of gas molecules.
	
	double dtau = opac*nden;
	if(j>0){
	  gastauproflo[i][j] = gastauproflo[i][j-1] + dtau * dl;
	  tauproflo[i][j] = tauproflo[i][j-1] + dtau * dl;
	}
	gastaulayerlo[i][j] = dtau * dl;
	taulayerlo[i][j] = dtau * dl;
	
	gasw0lo[i][j] = sca / opac;
	
	double dlc = 0.; // thickness of cloud inside the layer
	double hazeabund = 0.;
	double hazedepth = 0.;
	double absdepth = 0.;
	double scadepth = 0.;

	double absorbxsec = 0.;
	double scatterxsec = 0.;
	double hazexsec = 0.;

	doMie = false;
	
	if((hazetype!=0 && cloudmod>1) || cloudmod==4){
	  if(cloudmod!=4){
	    absorbxsec = atmoslores->getAbsXsec(wavel,haze[1]);
	    scatterxsec = atmoslores->getScaXsec(wavel,haze[1]);
	    hazexsec = absorbxsec + scatterxsec;

	    if(mode<=1 && streams==2){
	      doMie = true;
	    }
	  }
	  // ada: For the power-law opacity cloud model, we use the hazexsec variable as the wavelength scaling (haze[1] is the power-law exponent). haze[0] is a constant optical depth per unit length, or linear attenuation coefficient, for the cloud.
	  else{
	    hazexsec = pow(wavel, haze[0]);
	  }
	  
	  // Slab cloud model
	  if(cloudmod==2){
	    hazeabund = haze[0];
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
	    double P_top = haze[1];
	    double P_base = P_top * haze[2];
	    hazeabund = haze[3] * ((prprof[j+1]*prprof[j+1])-(prprof[j]*prprof[j])) / ((P_base*P_base)-(P_top*P_top));
	    // layer is strictly inside the cloud
	    if(prprof[j] >= P_top && prprof[j+1] <= P_base){
	      dlc = dl;
	    }
	    // layer overlaps bottom of cloud
	    else if(prprof[j] >= P_top && prprof[j] <= P_base){
	      dlc = (log10(P_base/prprof[j]) / log10(P_base/P_top));
	    }
	    // layer overlaps top of cloud
	    else if(prprof[j+1] >= P_top && prprof[j+1] <= P_base){
	      dlc = (log10(prprof[j+1]/P_top) / log10(P_base/P_top));
	    }
	    // layer contains entire cloud
	    else if(prprof[j] <= P_top && prprof[j+1] >= P_base){
	      dlc = (log10(prprof[j+1]/P_top) / log10(P_base/P_top));
	    }
	  }
	  
	  hazedepth = hazexsec * hazeabund;
	  absdepth = absorbxsec * hazeabund;
	  scadepth = scatterxsec * hazeabund;
	  
	} // end if(hazetype!=0)
        
	if(j==0){
	  tauproflo[i][j] = 0.;
	}
	else{
	  tauproflo[i][j] += hazedepth * dlc;
	}
	taulayerlo[i][j] += hazedepth * dlc;

	if(doMie==true){
	  // w0 is the single scattering albedo: the ratio of the scattering optical depth to the total optical depth
	  w0lo[i][j] = (sca*nden + scadepth) / (dtau + hazedepth);
	  
	  if(j==0){
	    backwardfrac = atmoslores->getAsym(wavel,haze[1]);
	    forwardfrac = 1.-backwardfrac;
	  }
	  // Hemispheric approximation to the asymmetry parameter integral.
	  asymlo[i][j] = (forwardfrac-backwardfrac)*scatterxsec / (sca + (forwardfrac+backwardfrac)*scatterxsec);
	}
	// ada: For the power-law opacity cloud model, the single-scattering albedo is a free parameter (constant in wavelength).
	else if(cloudmod==4){
	  w0lo[i][j] = haze[4];
	}
	else w0lo[i][j] = sca*nden / (dtau + hazedepth);
	// w0 is layer single scattering albedo with layer 0 on top
	// ratio of scattering to total opacity
      } // end for(j)
    } // end for(i)
  } // end if(table=="lores")
}
// end getTauProf
