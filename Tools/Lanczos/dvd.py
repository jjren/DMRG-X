#!/usr/bin/python
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

eAtodebye = 4.8032038
stateidx = 1

# read energy
f = open("dvde.out")
nstate = int(f.readline())
e = np.zeros((nstate))
for istate in range(0, nstate):
    e[istate] = float(f.readline())
f.close()

f1 = open("dvd_dyn.out")
dip = np.zeros((3,nstate))
for istate in range(0,nstate):
    line = f1.readline()
    for i in range(0,3):
        dip[i,istate] = float(line.split()[i]) / eAtodebye
f1.close()

#print e
#print dip

omega = np.arange(0.0, 40, 0.01)
npoints = omega.shape[0]
spectrum = np.zeros((3,npoints))

# set the figure
Fig = plt.figure()
ax = Fig.add_subplot(111)
ax.set_xlabel('eV')
ax.set_ylabel(r'$G(\omega)$')
ax.set_xlim(0.0, 40.0)  # set the xlim to xmin, xmax
ax.set_ylim(0.0, 4.5)   
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

color = ('b', 'g', 'r')
markers = ('o', 'v', '^')
label = ('x', 'y', 'z')


eta = 1.0e-1

for i in range(0,3):
    for j in range(0, nstate):
        spectrum[i,:] += eta*dip[i,j]**2 / ((omega+e[stateidx-1]-e[j])**2 + eta**2) / np.pi
    ax.plot(omega, spectrum[i,:], marker = markers[i],
                markersize=0, c=color[i],  label=label[i])
    ax.vlines(e[:]-e[stateidx-1], [0], dip[i,:]**2, colors=color[i], linestyle='solid',
            linewidth=1)
ax.legend(loc='upper center', shadow=True)
#plt.show()
Fig.savefig("state2-dvd50-256.eps", format='eps', dpi=400)
