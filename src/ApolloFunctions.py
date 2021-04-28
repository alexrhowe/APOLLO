import numpy as np
from astropy.convolution import convolve

def FindBands(wavehi,wavelo,flux,err):
    whibands = [[wavehi[0]]]
    wlobands = [[wavelo[0]]]
    fbands = [[flux[0]]]
    ebands = [[err[0]]]
    bindex = [[0]]
    length = len(wavehi)
    bi = 0
    for i in range(1,length):
        if wavehi[i]==wavelo[i-1]:
            whibands[bi].append(wavehi[i])
            wlobands[bi].append(wavelo[i])
            fbands[bi].append(flux[i])
            ebands[bi].append(err[i])
        else:
            bindex[bi].append(i-1)
            bindex.append([i])
            bi = bi+1
            whibands.append([wavehi[i]])
            wlobands.append([wavelo[i]])
            fbands.append([flux[i]])
            ebands.append([err[i]])
    bindex[bi].append(length-1)

    for i in range(0,bi+1):
        whibands[i] = np.asarray(whibands[i])
        wlobands[i] = np.asarray(wlobands[i])
        fbands[i] = np.asarray(fbands[i])
        ebands[i] = np.asarray(ebands[i])
    
    return bindex,whibands,wlobands,fbands,ebands

def ConvBands(bandflux,banderr,dataconv):
    if dataconv <= 1:
        convflux = bandflux
        converr = banderr
    else:
        convflux = []
        converr = []
        for i in range(0,len(bandflux)):
            convflux.append(ConvSpec(bandflux[i],dataconv))
            converr.append(ConvSpec(banderr[i],dataconv))
        
    return convflux,converr

def BinBands(bandhi,bandlo,convflux,converr,databin):
    if databin <= 1:
        binflux = convflux
        binerr = converr
        binhi = bandhi
        binlo = bandlo
    if databin > 1:
        binflux = convflux[0]
        binerr = converr[0]
        binhi = bandhi[0]
        binlo = bandlo[0]
        if len(convflux)>1:
            for i in range(1,len(convflux)):
                binflux = np.r_[binflux,convflux[i]]
                binerr = np.r_[binerr,converr[i]]
                binhi = np.r_[binhi,bandhi[i]]
                binlo = np.r_[binlo,bandlo[i]]
        else:
            binflux,binerr,binhi,binlo = BinSpec(convflux[0],converr[0],bandhi[0],bandlo[0],databin)
            if len(convflux)>1:
                for i in range(1,len(convflux)):
                    bf,be,bh,bl = BinSpec(convflux[i],converr[i],bandhi[i],bandlo[i],databin)
                    binflux = np.r_[binflux,bf]
                    binerr = np.r_[binerr,be]
                    binhi = np.r_[binhi,bh]
                    binlo = np.r_[binlo,bl]

    return binhi,binlo,binflux,binerr

def SliceModel(bandhi,bandlo,modwave):
    bindex = []
    modindex = [[0]]
    for i in range(0,len(bandhi)):
        bstart = bandhi[i][0]
        bend = bandlo[i][-1]
        js = np.where(modwave<bstart)[0][0]-1
        je = np.where(modwave<bend)[0][0]+1
        bindex.append([js,je])
        
    heads = []
    for i in range(0,len(bindex)):
        heads.append(bindex[i][0])
    hsort = np.argsort(heads)

    slwave = []
    for i in range(0,len(hsort)):
        js = bindex[hsort[i]][0]
        je = bindex[hsort[i]][1]
        slwave = np.r_[slwave,modwave[js:je]]
        modindex[i].append(modindex[i][0]+je-js)
        if i<len(hsort)-1:
            modindex.append([modindex[i][0]+je-js])
            
    return bindex,modindex,slwave

def NormSpec(wave,flux,startsnorm,endsnorm):    
    normwave = np.zeros(len(startsnorm))
    normpoints = np.zeros(len(startsnorm))

    for i in range(0,len(normpoints)):
        stemp = (int)(startsnorm[i])
        etemp = (int)(endsnorm[i])
        normwave[i] = (wave[stemp] + wave[etemp-1])/2.
        normpoints[i] = np.mean(flux[stemp:etemp])

    fit = np.polyfit(normwave,normpoints,len(normwave))
    poly = np.poly1d(fit)
    
    return flux/poly(wave)

def ConvSpec(flux,binw):
    
    binw6 = binw * 6.0
    sigmab = binw/2.35
    kwid = (int)(binw6)
    kernel = np.zeros(kwid)
    
    for i in range(0,kwid):
        kernel[i] = np.exp(-0.5*(i-binw6/2.)*(i-binw6/2.)/sigmab/sigmab)
    kernel = kernel/np.sum(kernel)

    convflux = np.convolve(flux,kernel,mode='same')
    #convflux = convolve(flux,kernel,boundary='extend')
    
    return convflux

def BinSpec(flux,err,wavehi,wavelo,binw):

    blen = (int)(len(flux)/binw)
    binflux = np.zeros(blen)
    binerr = np.zeros(blen)
    ibinhi = np.zeros(blen)
    ibinlo = np.zeros(blen)
    fbinhi = np.zeros(blen)
    fbinlo = np.zeros(blen)
    binhi = np.zeros(blen)
    binlo = np.zeros(blen)

    for i in range(0,blen):
        ibinhi[i] = i*binw
        ibinlo[i] = (i+1)*binw
        fbinhi[i] = np.modf(ibinhi[i])[0]
        fbinlo[i] = np.modf(ibinlo[i])[0]
        if fbinhi[i]==0.:
            ibinhi[i] = ibinhi[i] + 0.000001
            fbinhi[i] = 0.000001
        if fbinlo[i]==0.:
            ibinlo[i] = ibinlo[i] + 0.000001
            fbinlo[i] = 0.000001

    for i in range(0,len(ibinhi)):

        binflux[i] = np.sum(flux[(int)(np.ceil(ibinhi[i])):(int)(np.floor(ibinlo[i]))])
        binerr[i] = np.sum(err[(int)(np.ceil(ibinhi[i])):(int)(np.floor(ibinlo[i]))])
        
        binflux[i] = binflux[i] + (1.-fbinhi[i])*flux[(int)(np.floor(ibinhi[i]))]
        binerr[i] = binerr[i] + (1.-fbinhi[i])*err[(int)(np.floor(ibinhi[i]))]

        if (int)(np.floor(ibinhi[i]))==(int)(len(flux)-1):
            binhi[i] = (1.-fbinhi[i])*wavehi[(int)(np.floor(ibinhi[i]))]
        else:
            binhi[i] = (1.-fbinhi[i])*wavehi[(int)(np.floor(ibinhi[i]))] + fbinhi[i]*wavehi[(int)(np.floor(ibinhi[i]))+1]

        if (int)(np.floor(ibinlo[i]))==len(flux):
            binlo[i] = (1.-fbinlo[i])*wavelo[(int)(np.floor(ibinlo[i]))-1]
        else:
            binlo[i] = (1.-fbinlo[i])*wavelo[(int)(np.floor(ibinlo[i]))-1] + fbinlo[i]*wavelo[(int)(np.floor(ibinlo[i]))]

        if (int)(np.floor(ibinlo[i]))>=len(flux):
            binflux[i] = binflux[i]
            binerr[i] = binerr[i]
        else:
            binflux[i] = binflux[i] + fbinlo[i]*flux[(int)(np.floor(ibinlo[i]))]
            binerr[i] = binerr[i] + fbinlo[i]*err[(int)(np.floor(ibinlo[i]))]
        
        binflux[i] = binflux[i]/binw
        binerr[i] = binerr[i]/binw
        binerr[i] = binerr[i]/np.sqrt(binw-1.)

    return binflux,binerr,binhi,binlo
    
def BinModel(flux,binhi,binlo):
    binflux = np.zeros(len(binhi))
    fbinhi = (1.-np.modf(binhi)[0])-0.5
    fbinlo = np.modf(binlo)[0]-0.5
    binw = binlo-binhi
    for i in range(0,len(binhi)):
        
        binflux[i] = np.sum(flux[(int)(np.ceil(binhi[i])):(int)(np.ceil(binlo[i]))])
        binflux[i] = binflux[i] + fbinhi[i]*flux[(int)(np.floor(binhi[i]))]
        if (int)(np.ceil(binlo[i]))==len(flux):
            binflux[i] = binflux[i] + fbinlo[i]*flux[(int)(np.floor(binlo[i]))]
        else:
            binflux[i] = binflux[i] + fbinlo[i]*flux[(int)(np.ceil(binlo[i]))]

        if i==0 and binw[i]==0: binflux[i] = binflux[i]/binw[i+1]
        elif i==len(binhi)-1 and binw[i]==0: binflux[i] = binflux[i]/binw[i-1]
        else: binflux[i] = binflux[i]/binw[i]

    return binflux

def GetBins(specwave,obshi,obslo):
    binshi = np.zeros(len(obshi))
    binslo = np.zeros(len(obshi))
    # Not sure if I still might need this.
    '''
    if(obshi[0]>specwave[0] or obslo[-1]<specwave[-1]):
        print "Wavelength out of range."
        return [0.,0.]
    '''
    for i in range(0,len(obshi)):
        for j in range(0,len(specwave)):
            if(obshi[i]>specwave[j]):
                binshi[i] = float(j) + (obshi[i]-specwave[j-1])/(specwave[j]-specwave[j-1])
                break
        for j in range(0,len(specwave)):
            if(obslo[i]>=specwave[j]):
                binslo[i] = float(j) + (obslo[i]-specwave[j-1])/(specwave[j]-specwave[j-1])
                break
    return [binshi,binslo]

def GetScaOpac(gases,abunds):
    filler = 1. - np.sum(10**abunds)
    gaslist = ["h2he","h2","he","h-","h2o","ch4","co","co2","nh3","h2s","Burrows_alk","Lupu_alk","crh","feh","tio","vo","hcn","n2","ph3"]
    mmwlist = [2.28, 2.00, 4.00, 1.00, 18.0, 16.0, 28.0, 44.0, 17.0, 34.1, 24.1, 24.1, 53.0, 56.8, 63.9, 66.9, 27.0, 28.0, 34.0]
    scalist = [0.672e-27, 0.605e-27, 0.047e-27, 19.36e-27, 2.454e-27, 6.50e-27, 4.14e-27, 6.82e-27, 4.80e-27, 14.36e-27, 718.9e-27, 718.9e-27, 84.0e-27, 84.0e-27, 183.3e-27, 131.3e-27, 7.32e-27, 3.18e-27, 19.55e-27]
    '''
    Notes on scattering coefficients:
    Estimated h2 (=h2+he) from the mixing ratio.
    H-, CrH, TiO, and VO estimated based on theoretical models.
    FeH was not available; set equal to CrH.
    '''
    # Default filler is h2.
    if gases[0] in gaslist:
        i = gaslist.index(gases[0])
        mmw = mmwlist[i] * filler
        scaopac = scalist[i] * filler
    else:
        mmw = mmwlist[0] * filler
        scaopac = scalist[0] * filler
        
    for n in range(1,len(gases)):
        if gases[n] in gaslist:
            i = gaslist.index(gases[n])
            mmw = mmw + mmwlist[i] * 10**(abunds[n-1])
            scaopac = scaopac + scalist[i] * 10**(abunds[n-1])
        n = n+1
    return mmw,scaopac

def GetMollist(gases):
    mollist = np.zeros(len(gases))
    gaslist = ["h2he","h2","he","h-","h2o","ch4","co","co2","nh3","h2s","Burrows_alk","Lupu_alk","crh","feh","tio","vo","hcn","n2","ph3"]
    for i in range(0,len(gases)):
        if gases[i] in gaslist: mollist[i] = gaslist.index(gases[i])
        else: mollist[i] = 0
    return mollist
