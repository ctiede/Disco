import sys
import h5py as h5
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import discoUtil as du

def plotCheckpoint(file, vars=None, noGhost=False):
    
    print("Loading {0:s}...".format(file))

    t, r, phi, z, prim, dat = du.loadCheckpoint(file)
    rjph = dat[0]
    zkph = dat[1]
    primPhi0 = dat[2]
    piph = dat[3]

    nq = prim.shape[1]

    Zs = np.unique(z)
    z_eq = Zs[len(Zs)/2]
    eq_ind = (z==z_eq)
    title = "DISCO t = {0:.1f}".format(t)
    name = file.split('/')[-1].split('.')[0].split('_')[-1]

    if vars is None:
        vars = range(nq)

    for q in vars:
        print("   Plotting...")
        fig, ax = plt.subplots(1,1, figsize=(10,12),subplot_kw={'projection':'polar'})
        #ax.pcolormesh(rjph, thkph, primPhi0[:,:,q], cmap=plt.cm.inferno)
        vmin = prim[:,q].min()
        vmax = prim[:,q].max()
        if q in [5,6,7]:
            vmin, vmax = 0.0, 1.0
        Rs = np.unique(r)
        if noGhost:
            Rs = Rs[:-2]
        for i, R in enumerate(Rs):
            ind = r==R
            phif = np.empty(len(piph[ind])+1)
            phif[1:] = piph[ind]
            phif[0] = piph[ind][-1]
            C = ax.pcolormesh(phif, rjph[i:i+2], prim[None,ind,q], 
                    cmap=plt.cm.inferno, vmin=vmin, vmax=vmax)
        #ax.plot(np.linspace(0, 2*np.pi, 100), 2*np.ones(100), color='grey', lw=6, ls='-')
        #ax.plot(np.linspace(0, 2*np.pi, 100), 6*np.ones(100), color='grey', lw=6, ls='--')
        ax.set_theta_direction(-1)
        ax.set_theta_offset(0.5*np.pi)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        #ax.tick_params(axis='x', which='both', bottom='off', top='off', 
        #                labelbottom='off')
        #ax.tick_params(axis='y', which='both', bottom='off', top='off', 
        #                labelbottom='off')
        if q not in [5,6]:
            fig.colorbar(C)
        fig.suptitle(title, fontsize=24)
        fig.tight_layout()
        plotname = "plot_eq_{0:s}_lin_{1:d}.png".format(name, q)
        print("   Saving {0:s}...".format(plotname))
        fig.savefig(plotname)
        plt.close(fig)
        
        if q in [0,1]:
            fig, ax = plt.subplots(1,1, subplot_kw={'projection':'polar'})
            #ax.pcolormesh(rjph, thkph, primPhi0[:,:,q], cmap=plt.cm.inferno)
            for i, R in enumerate(Rs):
                ind = r==R
                phif = np.empty(len(piph[ind])+1)
                phif[1:] = piph[ind]
                phif[0] = piph[ind][-1]
                C = ax.pcolormesh(phif, rjph[i:i+2], prim[None,ind,q], 
                        cmap=plt.cm.inferno, vmin=vmin, vmax=vmax)
            #ax.set_theta_direction(-1)
            #ax.set_theta_offset(0.5*np.pi)
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            #ax.tick_params(axis='x', which='both', bottom='off', top='off', 
            #                labelbottom='off')
            #ax.tick_params(axis='y', which='both', bottom='off', top='off', 
            #            labelbottom='off')
            fig.colorbar(C)
            fig.suptitle(title, fontsize=18)
            fig.tight_layout()
            plotname = "plot_eq_{0:s}_log_{1:d}.png".format(name, q)
            print("   Saving {0:s}...".format(plotname))
            fig.savefig(plotname)
            plt.close(fig)

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Makes plots of Disco prims as a function of r.")
        print("usage: python plotDiscoR.py <checkpoint.h5 ...>")
        sys.exit()

    files = sys.argv[1:]
    for f in files:
        plotCheckpoint(f, [0,1,2,3], noGhost=False)
