#!/usr/bin/python
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

def getdata(dim, f):
    data = np.zeros((dim))
    count = 0
    while count < dim:
        line = f.readline()
        for i in range(0,len(line.split())):
            data[count+i] = float(line.split()[i])
        count += len(line.split())
    return data


def dyn(dim,npoints,omega,e,ev,diag,offdiag, dyn_initstate):
    # Lorentzian broaden
    eta = 1.0e-1

    absorb = np.zeros((npoints))

    # calculate the dynamic correlation function
    for ipoint in range(0,npoints):
        for ilanc in range(0,dim):
            absorb[ipoint] += ev[ilanc]*ev[ilanc]*eta / ((omega[ipoint]+e[dyn_initstate-1]-e[ilanc])**2+eta*eta)
    
    return absorb
    


def dyn_corr(dim, norm, figname, omega, dyn_initstate):   
 
    f = open("lanczos.out")
    # e: energy; ev: eigenvector coeff of initial \psi
    # diag: the diagonal element; offdiag: the off-diagonal elementl 
    e = np.zeros((3,dim))
    ev = np.zeros((3,dim))
    diag = np.zeros((3,dim))
    offdiag = np.zeros((3,dim-1))
    
    # get the e; ev; diag; offdiag from lanczos.out
    # can only treat the lanczos space of x,y,z component has equal size 
    section = -1
    line = f.readline()
    while line:
        if line.split()[0] == "=========================":
            section += 1 
        elif line.split()[0] == "eigenvalue":
            e[section,:] = getdata(dim,f)
        elif line.split()[0] == "eigenvector":
            ev[section,:] = getdata(dim,f)
        elif line.split()[0] == "diagonal":
            diag[section,:] = getdata(dim,f)
        elif line.split()[0] == "off-diagonal":
            offdiag[section,:] = getdata(dim-1,f)
        line = f.readline()
    f.close()
    
    # calculate the correlation function
    npoints = omega.shape[0]
    # initialize the correlation function
    absorb = np.zeros((3,npoints))
    
    # set the figure
    Fig = plt.figure()
    ax = Fig.add_subplot(111)
    ax.set_xlabel('eV')
    ax.set_ylabel(r'$G(\omega)$')
    ax.set_xlim(omega[0], omega[npoints-1])  # set the xlim to xmin, xmax
    ax.set_ylim(0.0, 4.5)   
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    markers = ('o', 'v', '^')
    colors = ('b', 'g', 'r')
    label = ('x', 'y', 'z')
    for i in range(0,3):
        # calculate the correlation function
        absorb[i,:] = dyn(dim,npoints,omega,e[i,:],ev[i,:],diag[i,:],offdiag[i,:],
                dyn_initstate) *  norm[i] / np.pi
        ax.plot(omega, absorb[i,:], marker = markers[i], markersize=0,
                c=colors[i], label=label[i])
    
    ax.legend(loc='upper center', shadow=True)
    
    #plt.show()
    Fig.savefig(figname+".eps", format='eps', dpi=400)
    
    return None


if __name__ == "__main__":
    # input 
    # the dimension of lanczos space
    dim = 100
    
    # set the omega range
    omega = np.arange(-0.1, 40, 0.005)
    
    # normalization of initial lanczos space (norm = <\psi|\psi>)
    norm = np.array([4.59992887193602, 2.88692686189203, 0.0])
    
    # set the init state index
    dyn_initstate = 1

    dyn_corr(dim, norm, "default", omega, dyn_initstate)
