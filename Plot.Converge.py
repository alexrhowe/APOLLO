from __future__ import print_function
import sys
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) > 1:
    fin = open(sys.argv[1], "r")
else:
    "Input file not specified."
    sys.exit()

line = fin.readline().split()
nwalkers = int(line[0])
nsteps = int(line[1])
ndim = int(line[2])

pnames = fin.readline().split()

lines = fin.readlines()
print(nwalkers, nsteps, ndim, len(lines))
samples = np.zeros((nsteps, nwalkers, ndim))
nsamp = np.linspace(1, nsteps, nsteps)

for i in range(0, nsteps):
    if i % 100 == 0:
        print(i)
    for j in range(0, nwalkers):
        for k in range(0, ndim):
            # print(len(lines[i*nwalkers+j]))
            samples[i, j, k] = lines[2 * i * nwalkers + 2 * j].split()[k]

for n in range(0, ndim):
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    plt.xlabel("Steps", fontsize=14)
    plt.ylabel(pnames[n], fontsize=14)
    for i in range(0, ndim):
        ax1.plot(nsamp, samples[:, i, n])
    plt.show()
