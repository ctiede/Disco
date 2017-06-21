import sys
import math
import h5py as h5
import matplotlib.pyplot as plt
import numpy as np

def loadCheckpoint(filename):

    f = h5.File(filename, "r")

    piph = f['Data']['Cells'][:,-1][...]
    prim = f['Data']['Cells'][:,:-1][...]
    index = f['Grid']['Index'][...]
    idPhi0 = f['Grid']['Id_phi0'][...]
    nphi = f['Grid']['Np'][...]
    t = f['Grid']['T'][0]
    riph = f['Grid']['r_jph'][...]
    ziph = f['Grid']['z_kph'][...]

    r = np.zeros(piph.shape)
    z = np.zeros(piph.shape)
    phi = np.zeros(piph.shape)
    R = 0.5*(riph[1:] + riph[:-1])
    Z = 0.5*(ziph[1:] + ziph[:-1])
    primPhi0 = np.zeros((index.shape[0], index.shape[1], prim.shape[1]))
    for k in range(index.shape[0]):
        for j in range(index.shape[1]):
            ind0 = index[k,j]
            ind1 = ind0 + nphi[k,j]
            r[ind0:ind1] = R[j]
            z[ind0:ind1] = Z[k]
            piph_strip = piph[ind0:ind1]
            pimh = np.roll(piph_strip, 1)
            pimh[pimh>piph_strip] -= 2*np.pi
            phi[ind0:ind1] = 0.5*(pimh+piph_strip)
            primPhi0[k,j,:] = prim[idPhi0[k,j],:]

    return t, r, phi, z, prim, (riph, ziph, primPhi0, piph)

def loadDiagRZ(filename):

    f = h5.File(filename, "r")

    t = f['Grid']['T'][0]
    rjph = f['Grid']['r_jph'][...]
    zkph = f['Grid']['z_kph'][...]
    diag = f['Data']['Poloidal_Diagnostics'][...]

    f.close()

    Nr = rjph.shape[0]-1
    Nz = zkph.shape[0]-1
    Nq = diag.shape[-1]
    #diag = np.resize(diag, (Nz,Nr,Nq))

    R = 0.5*(rjph[1:]+rjph[:-1])
    Z = 0.5*(zkph[1:]+zkph[:-1])

    r = np.empty((Nz,Nr))
    z = np.empty((Nz,Nr))

    r[:,:] = R[None,:]
    z[:,:] = Z[:,None]

    return t, r, z, diag, rjph, zkph

def plotAx(ax, x, y, xscale, yscale, xlabel, ylabel, *args, **kwargs):
    ax.plot(x, y, *args, **kwargs)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    ax.set_xscale(xscale)
    if (y>0).any():
        ax.set_yscale(yscale)
    else:
        ax.set_yscale("linear")
