from __future__ import print_function
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.special import expn
import corner

# Am I getting the pressure wrong here?

def getProfile(rad,grav,tpparams):
    Tint = tpparams[0]
    kIR = 10**tpparams[1]
    gamma1 = tpparams[2]
    gamma2 = tpparams[3]
    alpha = tpparams[4]
    
    taulist = np.zeros(101)
    for i in range(101):
        taulist[i] = 0.001 * 10**(10.*i/100.)
        
    Tirr = 0.
    prprof = np.zeros(101)
    hprof = np.zeros(101)
    rho = np.zeros(101)
    tprof = np.zeros(101)
    
    hprof[0] = 0.
    prprof[0] = grav*taulist[0]/kIR
    
    for i in range(0,101):
        E21 = expn(2,gamma1*taulist[i])
        E22 = expn(2,gamma2*taulist[i])
        xi1 = 2./3. + 2./(3.*gamma1) * (1.+(gamma1*taulist[i]/2. - 1.)*np.exp(-gamma1*taulist[i])) + 2.*gamma1/3.*(1.-taulist[i]*taulist[i]/2.)*E21
        xi2 = 2./3. + 2./(3.*gamma2) * (1.+(gamma2*taulist[i]/2. - 1.)*np.exp(-gamma2*taulist[i])) + 2.*gamma2/3.*(1.-taulist[i]*taulist[i]/2.)*E22
        tprof[i] = (3.*pow(Tint,4)/4.*(2./3.+taulist[i])) + (3.*pow(Tirr,4)/4.*(1-alpha)*xi1) + (3.*pow(Tirr,4)/4.*alpha*xi2)
        tprof[i] = pow(tprof[i],0.25)
        
        if tprof[i]<75.: tprof[i]=75.
        if tprof[i]>4000.: tprof[i]=4000.
        rho[i] = prprof[i]/tprof[i]*2.21/6.022e23/1.38e-16
        deltaH = 0
        
        if i<100:
            deltaH = (taulist[i+1]-taulist[i])/(kIR*rho[i])
            hprof[i+1] = hprof[i]+deltaH
            prprof[i+1] = prprof[i] + grav*rad*rad/(rad+hprof[i+1])**2 * rho[i]*deltaH;
    print(prprof)
    return tprof,prprof
            
if len(sys.argv)>1:
    fdata = open(sys.argv[1],'r')
else:
    print('Parameters: Input File, Object Name, [Reference Spectrum] OR [Confidence Interval]')
    sys.exit()
    
outfile = ''
if len(sys.argv)>2:
    outfile = outfile + sys.argv[2] + '.'

extra = ''
if len(sys.argv)>3:
    extra = sys.argv[3]

header = fdata.readline().split()
datatype = header[0]
    
if datatype == 'Spectrum':
    if extra != '':
        fref = open(extra,'r')  # Don't need this for samples.
    else:
        print('Error: Reference spectrum not specified.')
        sys.exit()
    
    # Processing the Apollo spectrum
    dlines = fdata.readlines()
    dlen   = len(dlines)
    dcalhi = np.zeros(dlen)
    dcallo = np.zeros(dlen)
    dflux  = np.zeros(dlen)
    dnoise = np.zeros(dlen)

    for i in range(0,dlen):
        dcalhi[i] = (float)(dlines[i].split()[0])
        dcallo[i] = (float)(dlines[i].split()[1])
        dflux[i]  = (float)(dlines[i].split()[2])
        dnoise[i] = (float)(dlines[i].split()[3])

    dcalmid = (dcalhi + dcallo)/2.
    dbandwidth = np.fabs(np.gradient(dcalmid)) * 2.99792458e10

    # Processing the reference spectrum
    
    rlines = fref.readlines()
    rlen   = len(rlines)
    rcalhi = np.zeros(rlen)
    rcallo = np.zeros(rlen)
    rflux  = np.zeros(rlen)
    rnoise = np.zeros(rlen)

    for i in range(0,rlen):
        rcalhi[i] = (float)(rlines[i].split()[0])
        rcallo[i] = (float)(rlines[i].split()[1])
        rflux[i]  = (float)(rlines[i].split()[2])
        rnoise[i] = (float)(rlines[i].split()[3])

    rcalmid = (rcalhi + rcallo)/2. + 0.0000 # Is that the wavelength offset?
    rbandwidth = np.fabs(np.gradient(rcalmid)) * 2.99792458e10

    # End input processing

    fr = interp1d(dcalmid,dflux,kind='linear')  # Cubic?
    bflux = fr(rcalmid)

    fd = interp1d(dcalmid,dnoise,kind='linear') # Cubic?
    nfac = float(rlen)/float(dlen)
    bnoise = fd(rcalmid) * np.sqrt(nfac)  # Note rcalmid

    # Output binned spectrum

    fout = open('modelspectra/binned.spectrum.dat','w') # Need a filename format. Easier as part of makespectrum.
    for i in range(0,rlen):
        fout.write('{0:8.2f} {1:8.2f} {2:8.5e} {3:8.5e} {3:8.5e} {4:8.5e}\n'.format(rcalhi[i],rcallo[i],bflux[i],bnoise[i],bflux[i]))

    # Plot spectra together

    xmin = max(min(dcalmid),min(rcalmid))
    xmax = min(max(dcalmid),max(rcalmid))
    xmin = xmin - 0.05*(xmax-xmin)
    xmax = xmax + 0.05*(xmax-xmin)

    dlist = [i for i in range(0,len(dcalmid)) if xmin < dcalmid[i] < xmax]
    rlist = [i for i in range(0,len(rcalmid)) if xmin < rcalmid[i] < xmax]

    yref = max(max(dflux[dlist]),max(rflux[rlist]))

    ymin = -0.20 * yref
    ymax =  1.05 * yref

    # Plot with binned spectrum

    fig = plt.figure(figsize=(10,7))
    ax = fig.add_subplot(111)
    plt.axis((xmin,xmax,ymin,ymax))
    
    plt.xlabel('$\lambda$ ($\mu$m)',fontsize=14)
    plt.ylabel('Flux (cgs)',fontsize=14)
    plt.tick_params(axis='both',which='major',labelsize=12)

    # May need a way to set the labels based on the input files

    ax.plot(rcalmid,bflux,'-',linewidth=1,label='Retrieved Spectrum (binned)',c='b')
    #ax.plot(rcalmid,rflux,'o',linewidth=1,label='Observations',c='k')
    ax.errorbar(rcalmid,rflux,rnoise,capsize=3,marker='o',linestyle='',linewidth=1,label='Observations',c='k')

    residuals = bflux-rflux
    ax.plot(rcalmid,residuals+ymin/2.,'-',linewidth=1,label='Residuals (binned and offset)',c='r')
    ax.plot([xmin,xmax],[0.,0.],'-',c='k')
    ax.plot([xmin,xmax],[ymin/2.,ymin/2.],'--',c='k')

    # Add average error bar symbol

    plt.legend(fontsize=12)
    figname = 'plots/' + outfile + 'binned.fit.png'
    plt.savefig(figname) # Need to procedurally generate filename

    # Plot with full-res spectrum
    
    fig = plt.figure(figsize=(10,7))
    ax = fig.add_subplot(111)
    plt.axis((xmin,xmax,ymin,ymax))
    
    plt.xlabel('$\lambda$ ($\mu$m)',fontsize=14)
    plt.ylabel('Flux (cgs)',fontsize=14)
    plt.tick_params(axis='both',which='major',labelsize=12)

    # May need a way to set the labels based on the input files

    ax.plot(dcalmid,dflux,'-',linewidth=1,label='Forward Model',c='b')
    #ax.plot(rcalmid,rflux,'o',linewidth=1,label='Observations',c='k')
    ax.errorbar(rcalmid,rflux,rnoise,capsize=3,marker='o',linestyle='',linewidth=1,label='Observations',c='k')

    residuals = bflux-rflux
    ax.plot(rcalmid,residuals+ymin/2.,'-',linewidth=1,label='Residuals (binned and offset)',c='r')
    ax.plot([xmin,xmax],[0.,0.],'-',c='k')
    ax.plot([xmin,xmax],[ymin/2.,ymin/2.],'--',c='k')

    # Add average error bar symbol

    plt.legend(fontsize=12)
    figname = 'plots/' + outfile + 'fullres.fit.png'
    #plt.savefig(figname) # Need to procedurally generate filename

    # Option to display plot
    plt.show()

elif datatype == 'Samples':
    confidence = 0.99
    
    if len(sys.argv)>4:
        if extra != '':
            confidence = float(extra) # Optional confidence interval parameter.
        
    # Processing the sample file    
    dim = [(int)(header[1]),(int)(header[2]),(int)(header[3])]
    profile = header[4]

    if profile == 'Layered':
        npress = (int)(header[5])
        minP = (float)(header[6])
        maxP = (float)(header[7])
    pnames = fdata.readline().split()
    
    baselist = ['Rad','Log(g)','Cloud_Base','P_cl','Mass','C/O','[Fe/H]','Teff']
    bnamelist = ['Radius (R$_J$)','log(g)','Base Pressure (bar)','Base Pressure (bar)','Mass (M$_J$)','C/O','[M/H]','T$_{eff}$ (K)']
    
    gaslist = ['h2','h2only','he','h-','h2o','ch4','co','co2','nh3','h2s','Burrows_alk','Lupu_alk','crh','feh','tio','vo','hcn','n2','ph3']
    gnamelist = ['H$_2$+He','H$_2$','He','[H-]','[H$_2$O]','[CH$_4$]','[CO]','[CO$_2$]','[NH$_3$]','[H$_2$S]','[Na,K]','[Na,K]','[CrH]','[FeH]','[TiO]','[VO]','[HCN]','[N2]','[PH3]']
    basic = []
    bnames = []
    bnames2 = []
    gases = [] 
    gnames = []

    tplist = ['Tint','log(kappaIR)','gamma1','gamma2','alpha']
    tpnamelist = ['T$_{int}$','log($\\kappa_{IR}$)','$\\gamma_1$','$\\gamma_2$','$\\alpha$']
    tp = []
    tpnames = []

    statlist = ['S0','S1','S2','S3','deltaL','deltaL0','deltaL1','deltaL2','deltaL3','logf','logf0','logf1','logf2','logf3','Likelihood']
    statnamelist = ['Band 1 Scale','Band 2 Scale','Band 3 Scale','Band 4 Scale','$\Delta\lambda$','$\Delta\lambda_1$','$\Delta\lambda_2$','$\Delta\lambda_3$','$\Delta\lambda_4$','log(f)','log($f_1$)','log($f_2$)','log($f_2$)','log($f_3$)','log(Likelihood)']
    stats = []
    statnames = []

    temps = []

    for i in range(0,len(pnames)):
        if pnames[i] in baselist:
            basic.append(i)
            j = baselist.index(pnames[i])
            bnames.append(bnamelist[j])
            bnames2.append(pnames[i])
                
        if pnames[i] in gaslist:
            gases.append(i)
            j = gaslist.index(pnames[i])
            gnames.append(gnamelist[j])
            
        if pnames[i] in tplist:
            tp.append(i)
            j = tplist.index(pnames[i])
            tpnames.append(tpnamelist[j])
            
        if pnames[i] in statlist:
            stats.append(i)
            j = statlist.index(pnames[i])
            statnames.append(statnamelist[j])
            
        if pnames[i][0]=='T' and not pnames[i]=='Teff': temps.append(i)

    samples = np.zeros((dim[0]*dim[1],dim[2]))
    for i in range(0,len(samples)):
        line = fdata.readline().split()
        for j in range(0,dim[2]):
            samples[i,j] = (float)(line[j])
            if np.isnan(samples[i,j]): samples[i,j]=0.
            if j==0: samples[i,j] = samples[i,j]/11.2
            
    bsamples = samples[:,basic]
    gsamples = samples[:,gases]
    tsamples = samples[:,temps]
    tpsamples = samples[:,tp]
    statsamples = samples[:,stats]
    
    # Used for plotting radius corrections based on photometry
    '''
    samples[:,0] = samples[:,0]*11.2*0.99083
    bsamples[:,0] = bsamples[:,0]*11.2*0.99083
    bsamples[:,2] = bsamples[:,2]*0.98175
    '''
    
    grange = np.zeros(len(gnames))
    for i in range(0,len(gnames)): grange[i]=confidence
    brange = np.zeros(len(bnames))
    for i in range(0,len(bnames)): brange[i]=confidence
    tprange = np.zeros(len(tpnames))
    for i in range(0,len(tpnames)): tprange[i]=confidence
    statrange = np.zeros(len(statnames))
    for i in range(0,len(statnames)): statrange[i]=confidence
    
    fig = corner.corner(gsamples,labels=gnames,range=grange,plot_datapoints=False,show_titles=False,title_kwargs={'fontsize': 14},label_kwargs={'fontsize': 18})
    fig.subplots_adjust(left=0.10,bottom=0.10,wspace=0,hspace=0)
    fig1name = 'plots/' + outfile + 'gases.png'
    fig.savefig(fig1name)

    fig2 = corner.corner(bsamples,labels=bnames,range=brange,plot_datapoints=False,show_titles=False,title_kwargs={'fontsize': 14},label_kwargs={'fontsize': 18})
    fig2.subplots_adjust(left=0.10,bottom=0.10,wspace=0,hspace=0)
    fig2name = 'plots/' + outfile + 'basic.png'
    fig2.savefig(fig2name)

    if profile=='Parametric':
        fig3 = corner.corner(tpsamples,labels=tpnames,range=tprange,plot_datapoints=False,label_kwargs={'fontsize': 18})
        fig3.subplots_adjust(left=0.10,bottom=0.10,wspace=0,hspace=0)
        fig3name = 'plots/' + outfile + 'tpparams.eps'
        fig3.savefig(fig3name)
    
    #statsamples2 = [s for s in statsamples if not np.isinf(s[-1])]
    fig4 = corner.corner(statsamples,labels=statnames,range=statrange,plot_datapoints=False,labelsize=24)
    fig4.subplots_adjust(left=0.10,bottom=0.10,wspace=0,hspace=0)
    fig4name = 'plots/' + outfile + 'stats.eps'
    fig4.savefig(fig4name)
    
    if profile=='Layered':
        plist = np.zeros(len(temps))
        for i in range(0,len(temps)): plist[i] = 10**(maxP + (minP-maxP)*i/(len(temps)-1))
        tlist = np.percentile(tsamples,[16,50,84],axis=0)
        
        # Plot the T-P profile
        fit = np.percentile(samples,50,axis=0)
        
        fig5 = plt.figure(figsize=(8,6))
        ax = fig5.add_subplot(111)
        plt.axis((0,3000,10**(maxP),10**(minP)))
        ax.set_yscale('log')
        plt.xlabel('T (K)',fontsize=14)
        plt.ylabel('P (bar)', fontsize=14)
        
        ax.fill_betweenx(plist,tlist[0],tlist[2],facecolor='#ff8080')
        ax.plot(tlist[0],plist,c='r')
        ax.plot(tlist[1],plist,c='k')
        ax.plot(tlist[2],plist,c='r')
        
        fig5name = 'plots/' + outfile + 'TP.png'
        fig5.savefig(fig5name)

    elif profile == 'Parametric':
        tpplist = np.percentile(tpsamples,[16,50,84],axis=0)
        blist = np.percentile(bsamples,[16,50,84],axis=0)
    
        fig5 = plt.figure(figsize=(8,6))
        ax = fig5.add_subplot(111)
        ax.set_yscale('log')
        plt.xlabel('T (K)',fontsize=14)
        plt.ylabel('P (bar)', fontsize=14)

        print(blist)
        #print(plist)
        for i in range(0,20):
            j = np.random.randint(0,len(tpsamples))
            tpparams = tpsamples[j]
            #if 'Rad' in bnames:
            if 'Rad' in pnames:
                jr = pnames.index('Rad')
                rad = 6.371e8*samples[j][jr]
            else: rad = 6.371e8*11.2
            #if 'logg' in bnames:
            if 'Log(g)' in pnames:
                jg = pnames.index('Log(g)')
                grav = 10**samples[j][jg]
            else: grav = 10**4.87
            tprof,prprof = getProfile(rad,grav,tpparams)
            print(rad,grav,tpparams)
            ax.plot(tprof,prprof/1.e6,c='#a0a0a0')
            if i==19: ax.plot(tprof,prprof/1.e6,c='k')
            #print(prprof)
            
        ### WARNING: This is specific to GJ 229B!
        
        tpparams = tpplist[1]

        tpparams = [868.76,-1.54,0.89,1.14,0.46]
        tpparams = [869.18,-1.5443,0.902,1.090,0.493]
        
        tlayer = [[2641.52931,2712.96167,2784.73001],
		  [2260.05527,2297.38351,2333.51886],
		  [1867.77311,1884.83269,1900.34831],
		  [1446.35641,1452.16609,1459.30833],
		  [1349.21063,1354.03474,1358.71897],
		  [1019.16816,1028.61313,1037.69971],
		  [853.39312, 859.72377, 866.04291],
		  [712.28971, 724.35058, 735.58452],
		  [579.57611, 599.35291, 618.54968],
		  [482.02401, 515.34428, 542.83721],
		  [400.01529, 448.30752, 491.61250],
		  [330.63486, 394.27780, 458.57395],
		  [269.87095, 346.78162, 432.19712],
		  [214.19864, 304.03192, 418.90273],
		  [153.04883, 265.81270, 407.33757]]

        minP = -3.
        maxP = 2.5
        plist = np.zeros(len(tlayer))
        for i in range(0,len(tlayer)):
            plist[i] = 10**(maxP + (minP-maxP)*i/(len(tlayer)-1))
            
        tmins = []
        tmeds = []
        tmaxs = []
        for i in range(0,15):
            tmins.append(tlayer[i][0])
            tmeds.append(tlayer[i][1])
            tmaxs.append(tlayer[i][2])
            
        print(blist[0])
        print(blist[1])
        #rad = 6.371e8*blist[1][0]*11.2
        #grav = 10**blist[1][1]
        rad = 6.371e8*11.57
        grav = 10**4.514

        print(plist)
        print(tlayer[:][0])
        
        ax.fill_betweenx(plist,tmins,tmaxs,facecolor='#ff8080')
        ax.plot(tmins,plist,c='r')
        ax.plot(tmeds,plist,c='k')
        ax.plot(tmaxs,plist,c='r')
        
        tmid,pmid = getProfile(rad,grav,tpparams)
        print(pmid)
        print(rad,grav,tpparams)
                            
        plt.axis((0,3000,min(10**2.5,prprof[-1]/1.e6),min(10**-3,prprof[0]/1.e6)))
        #plt.axis((0,4000,10**2.5,10**-3))
        #ax.plot(tmid,pmid/1.e6,c='k')
    
        fig5name = 'plots/' + outfile + 'TP.png'
        fig5.savefig(fig5name)
        
    print('1-sigma bounds on basic parameters')
    blist = np.percentile(samples,[16,50,84],axis=0)
    for i in range(0,len(blist[0])):
        if len(pnames[i])<=6: print('{0:s}:\t\t   {1:6.5f}\t {2:6.5f}\t {3:6.5f}'.format(pnames[i],blist[0][i],blist[1][i],blist[2][i]))
        else: print('{0:s}:\t\t   {1:6.5f}\t {2:6.5f}\t {3:6.5f}'.format(pnames[i],blist[0][i],blist[1][i],blist[2][i]))