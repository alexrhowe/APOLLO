import numpy as np
import matplotlib.pyplot as plt
import sys

fin = open('modelspectra/example.Resolved.27params30k.ensemble.dat')
lines = fin.readlines()

params = []
pvals = []

nlines = len(lines)
l0 = 0

normmin = -1
normmax = -1
normtomax = False

print(len(sys.argv))

if len(sys.argv)==2:
    if float(sys.argv[1])==0: normtomax = True
elif len(sys.argv)>2:
    print(sys.argv)
    normmin = float(sys.argv[1])
    normmax = float(sys.argv[2])
    
for i in range(0,nlines):
    if lines[i][0]==' ':
        params.append(lines[i].split()[0])
        pvals.append(lines[i].split()[1:])        
    else:
        nlines = nlines-i
        l0 = i
        break

lines2 = lines[l0:]

print(params)
print(pvals[0][0:16])
print(pvals[1][0:16])
print(pvals[2][0:16])
print(pvals[3][0:16])

wave = np.zeros(nlines)
spectra = np.zeros((len(pvals[0]),nlines))

colors = []
for i in range(0,len(pvals[0])):
    r=0
    g=0
    b=0
    h2o = i%5
    ch4 = int(i/5)%5
    co = int(i/25)%5
    co2 = int(i/125)%5
    if h2o==0:
        r = 1.0
    if h2o==1:
        r = 1.0
        g = 1.0
    if h2o==2:
        g = 1.0
    if h2o==3:
        g = 1.0
        b = 1.0
    if h2o==4:
        b = 1.0
    color = [r*(ch4+1.)/5.,g*(ch4+1.)/5.,b*(ch4+1.)/5.]
    colors.append(color)

testlist = np.zeros(len(pvals[0]))
for i in range(0,nlines):
    line = lines2[i].split()
    wave[i] = 10000./float(line[0])
    for j in range(0,len(pvals[0])):
        spectra[j][i] = float(line[j+2])
        if i==0: testlist[j] = line[j+2]

if wave[0]<normmin<wave[-1] and wave[0]<normmax<wave[-1] and normmin<normmax:
    print('Normalizing to range ({0:f},{1:f})'.format(normmin,normmax))
    print(normmin,normmax,wave)
    wmin = np.where(np.logical_and(normmin<wave,wave<normmax))
    for i in range(0,len(pvals[0])):
        spectra[i] = spectra[i]/np.mean(spectra[i][wmin])
elif normtomax:
    print('Normalizing to maximum values.')
    for i in range(0,len(pvals[0])):
        spectra[i] = spectra[i]/np.max(spectra[i])
elif normmin==-1 and normmax==-1:
    print('No normalization. Original values used.')
else:
    print('Invalid normalization.')

fig = plt.figure()
ax = fig.add_subplot(111)

width = (wave[-1]-wave[0])/(len(pvals[0])-1)
nx = (wave[-1]-wave[0])/nlines

for i in range(0,len(pvals[0])):
    #if i%25==18:
    ax.plot(wave,spectra[i],linewidth=0.2,c=colors[i])
    #x = i*width
    #ind = int(i*nlines/len(pvals[0]))
    #ax.text(wave[ind],spectra[i][ind],'{0:d}'.format(i))

# NEED TO ALSO NORMALIZE OBSERVATIONS

fin2 = open('../Apollo0113/GJ570D.YJHK.dat','r')
#fin2 = open('examples/example.obs.dat')
lines2 = fin2.readlines()
lendat = len(lines2)
wdat = np.zeros(lendat)
fdat = np.zeros(lendat)
for i in range(0,lendat):
    wdat[i] = 10000./float(lines2[i].split()[0])
    fdat[i] = float(lines2[i].split()[2])
#fdat = fdat/np.mean(fdat[16:30])
fdat = fdat/np.mean(fdat[854:1328])

ax.plot(wdat,fdat,'o--',c='m')
    
plt.show()
