#!/usr/bin/python
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

# independent electron dynamic correlation function calculation

nocc = 5
nvir = 5

f = open("dipint.out")
norbs = int(f.readline())
if nocc + nvir != norbs :
    sys.exit("Error message")
dipint = np.zeros((norbs,norbs,3))
line = f.readline()
while line:
    for i in range(0,3):
        dipint[int(line.split()[0])-1, int(line.split()[1])-1, i] = float(line.split()[i+2])
        dipint[int(line.split()[1])-1, int(line.split()[0])-1, i] = float(line.split()[i+2])
    line = f.readline()

#print dipint
f.close()

f2 = open("MO.out")
emo = np.zeros((norbs))
line = f2.readline()
while line:
    if len(line.split()) == 2 :
       emo[int(line.split()[0])-1] = float(line.split()[1]) 
    line = f2.readline()

#print emo
f2.close()

eex = np.zeros((nocc*nvir))
dip2 = np.zeros((nocc*nvir,3))

# absorption line

i = 0
for iocc in range(0, nocc):
    for ivir in range(nocc, norbs):
        eex[i] = (emo[ivir] - emo[iocc])
        dip2[i,:] = 2.0 * dipint[iocc,ivir,:]**2
        i += 1

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
    for j in range(0, nocc*nvir):
        spectrum[i,:] += eta*dip2[j,i] / ((omega-eex[j])**2 + eta**2) / np.pi
    plt.plot(omega, spectrum[i,:], marker = markers[i],
                markersize=0, c=color[i],  label=label[i])
    plt.vlines(eex, [0], dip2[:,i], colors=color[i], linestyle='solid',
            linewidth=1)
ax.legend(loc='upper center', shadow=True)
#plt.show()
Fig.savefig("exact.eps", format='eps', dpi=400)
