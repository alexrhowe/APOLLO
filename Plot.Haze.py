from __future__ import print_function
import sys
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv)<=3:
    print('Error: arguments not specified. File1, File2, Particle Size in microns.')
    sys.exit()

file1 = sys.argv[1]
file2 = sys.argv[2]
size = np.log10((float)(sys.argv[3]))
    
table1 = open(file1,'r')
header1 = table1.readline().split()

len1  =   (int)(header1[0])
min1  = (float)(header1[1])
max1  = (float)(header1[2])
ns1   =   (int)(header1[3])
minS1 = (float)(header1[4])
maxS1 = (float)(header1[5])

if size<minS1 or size>maxS1:
    print('Error: particle size out of range.')

is1 = (int)(np.floor((size-minS1)/(maxS1-minS1)*(ns1-1)))

x1 = np.exp(np.linspace(np.log(min1),np.log(max1),len1))
y1 = np.zeros(len1)

for i in range(0,len1):
    line = table1.readline().split()
    y1[i] = (float)(line[is1+1])
table1.close()

table2 = open(file2,'r')
header2 = table2.readline().split()

len2  =   (int)(header2[0])
min2  = (float)(header2[1])
max2  = (float)(header2[2])
ns2   =   (int)(header2[3])
minS2 = (float)(header2[4])
maxS2 = (float)(header2[5])

if size<minS2 or size>maxS2:
    print('Error: particle size out of range.')
    sys.exit()
    
is2 = (int)(np.floor((size-minS2)/(maxS2-minS2)*(ns2-1)))

x2 = np.exp(np.linspace(np.log(min2),np.log(max2),len2))
y2 = np.zeros(len2)

for i in range(0,len2):
    line = table2.readline().split()
    y2[i] = (float)(line[is2+1])
'''
for i in range(0,ns2):
    line = table2.readline()
    if i==is2:
        line = line.split()
        for k in range(0,len2):
            y2[k] = (float)(line[k])
'''
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
