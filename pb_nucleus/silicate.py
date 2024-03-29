""" This code plots spectra extracted from the crntre of the nucleus and the north region of M31
in two panels.
You should have .tbl files obtained from CUBISM for both regions to run this code."""

import numpy as np                     # Allows Numpy functions to be called directly
import matplotlib.pyplot as plt         # Graphing routines
import sys
import scipy 
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import matplotlib


def getspectrum(sl1,sl2,ll2):
    # Read spectra from 3 modules and stitch them.
    # Returns the final stitched spectra

    W_SL2,F_SL2,SL2_unc  = sl2[:,0] , sl2[:,1],sl2[:,2]
    W_SL1,F_SL1,SL1_unc  = sl1[:,0] , sl1[:,1], sl1[:,2]
    W_LL2,F_LL2,LL2_unc  = ll2[:,0] , ll2[:,1], ll2[:,2]
    
    Y = []
    X = []
    Yerr = []


    # The below For loops cut the spectra at specific wavelengths and merge them into a
    # one array. Cutting points were determined by observing the unstitched spectrum.
    # It adds wavelengths in to X, Flux into Y and uncertainty into Yerr.
    # Note : We ignored SL3 spectrum

    for (i,aval) in enumerate(W_SL2):
        if aval < 7.57 :
            Y.append(F_SL2[i])
            X.append(aval)
            Yerr.append(SL2_unc[i])
      
    for (i,aval) in enumerate(W_SL1):
        if (aval > 7.57) & (aval < 14.50):
            Y.append(F_SL1[i])
            X.append(aval)
            Yerr.append(SL1_unc[i])

    for (i,aval) in enumerate(W_LL2):
        if aval > 14.50 :
            Y.append(F_LL2[i])
            X.append(aval)
            Yerr.append(LL2_unc[i])

    return X,Y,Yerr



def plotting(X,Y,Yerr,panel,axes,col=True):

    # Plot spectra into panels.
    if col:
        axes[panel].errorbar(X,Y,Yerr,0,'-')
    else:
        axes[panel].errorbar(X,Y,Yerr,0,'-', color='k')
    axes[panel].tick_params(axis='both', which='major', labelsize=24)
    axes[panel].tick_params(axis='both', which='minor', labelsize=24)

def doplot(color=True):
    matplotlib.rcdefaults()
    fig,axes = plt.subplots(2,1,sharex=True)
    minorLocator = AutoMinorLocator(5)
    for ax in axes:
        ax.xaxis.set_minor_locator(minorLocator)
        ax.tick_params(which='both', width=2)
        ax.tick_params(which='major', length=10)
        ax.tick_params(which='minor', length=7, color='k')
        ax.tick_params(which='both',labelsize=24)
        ax.tick_params(axis='both', which='major', labelsize=24)
        ax.tick_params(axis='both', which='minor', labelsize=24)
    fig.text(0.03, 0.5, 'Intensity (MJy/sr)', ha='center', va='center', rotation='vertical',fontsize=30)
    
    # Reading .tbl files got from CUBISM 
    sl2 = np.loadtxt("m31nuc_sl2_nucCentre.tbl", skiprows = 15 )
    sl1 = np.loadtxt("m31nuc_sl1_nucCentre.tbl", skiprows = 15 )
    ll2 = np.loadtxt("m31nuc_ll2_nucCentre.tbl", skiprows = 15 )
    Wavelength,Flux,FluxUnc = getspectrum(sl1,sl2,ll2)
    
    plotting(Wavelength,Flux,FluxUnc,1,axes,col=color)
    axes[1].set_xlabel("Wavelength ($\mu m$)",fontsize=24)
    axes[1].annotate('Silicate', xy=(9.7, 50),  xycoords='data',xytext=None,size=30, textcoords='offset points',arrowprops=dict(arrowstyle="->", linewidth = 4))
    axes[1].set_ylim(10,88)
    axes[1].yaxis.set_major_locator(MultipleLocator(20))    

#    sl2 = np.loadtxt("m31nuc_sl2_nucUP.tbl", skiprows = 15 )
#    sl1 = np.loadtxt("m31nuc_sl1_nucUP.tbl", skiprows = 15 )
#    ll2 = np.loadtxt("m31nuc_ll2_nucUP.tbl", skiprows = 15 )
    
    Wavelength,Flux,FluxUnc = np.loadtxt('m31nuc_nucUP_correct.dat',usecols=(0,1,2),unpack=True)
    
    plotting(Wavelength,Flux,FluxUnc,0,axes,col=color)
    axes[0].set_ylim(9,45)
    axes[0].yaxis.set_major_locator(MultipleLocator(10))

    fig.subplots_adjust(left=0.15, bottom=0.16, hspace=0.0)
    plt.show()
    return
