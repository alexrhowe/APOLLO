import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import corner

if(sys.argv)>1:
    fdata = open(sys.argv[1],'r')
else:
    print 'Parameters: Input File, Object Name, Data Type, [Reference Spectrum]'
    sys.exit()

outfile = ''
if len(sys.argv)>2:
    outfile = outfile + sys.argv[2] + '.'
    
datatype = ''
if len(sys.argv)>3:
    datatype = sys.argv[3]

if datatype == '':
    print 'Error: Input type not specified. Spectrum or Samples.'
    sys.exit()

if datatype == 'Spectrum':
    if len(sys.argv)>4:
        fref = open(sys.argv[4],'r')  # Don't need this for samples.
    else:
        print 'Error: Reference spectrum not specified.'
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

    dcalmid = 10000./((dcalhi + dcallo)/2.)
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

    rcalmid = 10000./((rcalhi + rcallo)/2.) + 0.0000 # Is that the wavelength offset?
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

    ax.plot(dcalmid,dflux,'-',linewidth=1,label='Retrieved Spectrum',c='b')
    #ax.plot(rcalmid,rflux,'o',linewidth=1,label='Observations',c='k')
    ax.errorbar(rcalmid,rflux,rnoise,capsize=3,marker='o',linestyle='',linewidth=1,label='Observations',c='k')

    residuals = bflux-rflux
    ax.plot(rcalmid,residuals+ymin/2.,'-',linewidth=1,label='Residuals (binned and offset)',c='r')
    ax.plot([xmin,xmax],[0.,0.],'-',c='k')
    ax.plot([xmin,xmax],[ymin/2.,ymin/2.],'--',c='k')

    # Add average error bar symbol

    plt.legend(fontsize=12)
    figname = 'plots/' + outfile + 'fullres.fit.png'
    plt.savefig(figname) # Need to procedurally generate filename

    # Option to display plot
    #plt.show()

elif datatype == 'Samples':
    confidence = 0.99
    
    if len(sys.argv)>4:
        confidence = (float)(sys.argv[4])  # Optional confidence interval parameter.
        
    # Processing the sample file    
    line = fdata.readline().split()
    dim = [(int)(line[0]),(int)(line[1]),(int)(line[2])]
    if len(line)>3:
        npress = (int)(line[3])
        minP = (float)(line[4])
        maxP = (float)(line[5])
    pnames = fdata.readline().split()
        
    baselist = ['Rad','RtoD','RtoD2U','Log(g)','Cloud_Base','P_cl','Mass','C/O','[Fe/H]','Teff']
    bnamelist = ['Radius (R$_J$)','Radius (R$_J$)','Radius (R$_J$)','log(g)','Base Pressure (bar)','Base Pressure (bar)','Mass (M$_J$)','C/O','[Fe/H]','T$_{eff}$ (K)']
    
    gaslist = ['h2','h2only','he','h-','h2o','ch4','co','co2','nh3','h2s','Burrows_alk','Lupu_alk','crh','feh','tio','vo','hcn','n2','ph3']
    gnamelist = ['H$_2$+He','H$_2$','He','[H-]','[H$_2$O]','[CH$_4$]','[CO]','[CO$_2$]','[NH$_3$]','[H$_2$S]','[Na,K]','[Na,K]','[CrH]','[FeH]','[TiO]','[VO]','[HCN]','[N2]','[PH3]']
    basic = []
    bnames = []
    bnames2 = []
    gases = [] 
    gnames = []
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
        if pnames[i][0]=='T' and not pnames[i]=='Teff': temps.append(i)

    samples = np.zeros((dim[0]*dim[1],dim[2]))
    for i in range(0,len(samples)):
        line = fdata.readline().split()
        for j in range(0,dim[2]):
            samples[i,j] = (float)(line[j])
            
    bsamples = samples[:,basic]
    gsamples = samples[:,gases]
    tsamples = samples[:,temps]
    
    grange = np.zeros(len(gnames))
    for i in range(0,len(gnames)): grange[i]=confidence
    brange = np.zeros(len(bnames))
    for i in range(0,len(bnames)): brange[i]=confidence
    
    fig = corner.corner(gsamples,labels=gnames,range=grange,plot_datapoints=False,labelsize=24)
    fig.subplots_adjust(left=0.10,bottom=0.10,wspace=0,hspace=0)
    fig1name = 'plots/' + outfile + 'gases.png'
    fig.savefig(fig1name)

    fig2 = corner.corner(bsamples,labels=bnames,range=brange,plot_datapoints=False,labelsize=24)
    fig2.subplots_adjust(left=0.10,bottom=0.10,wspace=0,hspace=0)
    fig2name = 'plots/' + outfile + 'basic.png'
    fig2.savefig(fig2name)

    plist = np.zeros(len(temps))
    for i in range(0,len(temps)): plist[i] = 10**(maxP + (minP-maxP)*i/(len(temps)-1))
    tlist = np.percentile(tsamples,[16,50,84],axis=0)
    
    # Plot the T-P profile
    fit = np.percentile(samples,50,axis=0)
    
    fig3 = plt.figure(figsize=(8,6))
    ax = fig3.add_subplot(111)
    plt.axis((0,3000,10**(maxP),10**(minP)))
    ax.set_yscale('log')
    plt.xlabel('T (K)',fontsize=14)
    plt.ylabel('P (bar)', fontsize=14)
    
    ax.fill_betweenx(plist,tlist[0],tlist[2],facecolor='#ff8080')
    ax.plot(tlist[0],plist,c='r')
    ax.plot(tlist[1],plist,c='k')
    ax.plot(tlist[2],plist,c='r')
    
    fig3name = 'plots/' + outfile + 'TP.png'
    fig3.savefig(fig3name)

    print '1-sigma bounds on basic parameters'
    blist = np.percentile(bsamples,[16,50,84],axis=0)
    for i in range(0,len(blist[0])):
        if len(bnames2[i])<=6: print '{0:s}:\t\t {1:10.5f} {2:10.5f} {3:10.5f}'.format(bnames2[i],blist[0][i],blist[1][i],blist[2][i])
        else: print '{0:s}:\t {1:10.5f} {2:10.5f} {3:10.5f}'.format(bnames2[i],blist[0][i],blist[1][i],blist[2][i])
    
else:
    print 'Error: input type wrongly specified. Spectrum or Samples.'
    sys.exit()
        
