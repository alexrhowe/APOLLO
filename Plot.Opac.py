from __future__ import print_function
import sys
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv)<=3:
    print('Error: arguments not specified. File1, File2, Pressure, Temperature')
    sys.exit()

file1 = sys.argv[1]
file2 = sys.argv[2]
press = np.log10((float)(sys.argv[3]))
temp = np.log10((float)(sys.argv[4]))
    
table1 = open(file1,'r')
header1 = table1.readline().split()

np1   =   (int)(header1[0])
minP1 = (float)(header1[1])
maxP1 = (float)(header1[2])
nt1   =   (int)(header1[3])
minT1 = (float)(header1[4])
maxT1 = (float)(header1[5])
len1  =   (int)(header1[6])
min1  = (float)(header1[7])
max1  = (float)(header1[8])

if press<minP1 or press>maxP1:
    print('Error: pressure out of range.')
    sys.exit()
if temp<minT1 or temp>maxT1:
    print('Error: temperature out of range.')
    sys.exit()

ip1 = (int)(np.floor((press-minP1)/(maxP1-minP1)*(np1-1)))
jt1 = (int)(np.floor((temp-minT1)/(maxT1-minT1)*(nt1-1)))

x1 = np.exp(np.linspace(np.log(min1),np.log(max1),len1))
y1 = np.zeros(len1)

for i in range(0,np1):
    for j in range(0,nt1):
        line = table1.readline()
        if i==ip1 and j==jt1:
            line = line.split()
            for k in range(0,len1):
                y1[k] = (float)(line[k])
table1.close()

table2 = open(file2,'r')
header2 = table2.readline().split()

np2   =   (int)(header2[0])
minP2 = (float)(header2[1])
maxP2 = (float)(header2[2])
nt2   =   (int)(header2[3])
minT2 = (float)(header2[4])
maxT2 = (float)(header2[5])
len2  =   (int)(header2[6])
min2  = (float)(header2[7])
max2  = (float)(header2[8])

if press<minP2 or press>maxP2:
    print('Error: pressure out of range.')
    sys.exit()
if temp<minT2 or temp>maxT2:
    print('Error: temperature out of range.')
    sys.exit()
ip2 = (int)(np.floor((press-minP2)/(maxP2-minP2)*(np2-1)))
jt2 = (int)(np.floor((temp -minT2)/(maxT2-minT2)*(nt2-1)))

x2 = np.exp(np.linspace(np.log(min2),np.log(max2),len2))
y2 = np.zeros(len2)

for i in range(0,np2):
    for j in range(0,nt2):
        line = table2.readline()
        if i==ip2 and j==jt2:
            line = line.split()
            for k in range(0,len2):
                y2[k] = (float)(line[k])
table2.close()
                
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_yscale('log')
plt.xlabel('Wavelength (microns)')
plt.ylabel('Cross Section (cm$^2$)')
#plt.axis((0.6,2.5,1e-31,1e-10))

ax.plot(x1,y1,'-',linewidth=1,c='k')
ax.plot(x2,y2,'-',linewidth=1,c='r')

plt.show()
