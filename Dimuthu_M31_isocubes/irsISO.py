""" Plots the ISO cam spectra from Nucleus, Bulge and the Region 9 with its IRS
spectra in three panels. You must have iso_nuc.py, iso_bulge.py and iso_isocvf.py codes in the same
directory with their data files to run this code.

Modified by PB, 2015/08/17:
- make figure greyscale rather than colour
- change references to pyfits to astrpy.io.fits
- remove tabs in indentation
(-resisted urge to refactor code) 
"""

import numpy as np                     
import matplotlib.pyplot as plt         
import pylab as py
from iso_nuc import plotvaluesNuc
from iso_bulge import plotvaluesBul
from iso_isocvf import plotvaluesIso



def plotting():

    Xiso,Yiso,Xirs,Yirs = plotvaluesNuc()
    axes[0].plot(Xiso,Yiso,'-',linewidth=3,color='0.6', label = "ISOCAM")
    axes[0].plot(Xirs,Yirs,'-',linewidth=3,color='k', label = 'IRS')
    axes[0].legend( loc='upper right', prop={'size':15} )
    axes[0].annotate('Nucleus', xy=(18, 30),  xycoords='data',xytext=None,size=20)
    axes[0].tick_params(axis='y', labelsize=20)
    axes[0].set_ylim(0,70)    

    Xiso,Yiso,Xirs,Yirs = plotvaluesBul()
    axes[1].plot(Xiso,Yiso,'-',linewidth=3, color='0.6')
    axes[1].plot(Xirs,Yirs,'-',linewidth=3, color='k')
    axes[1].annotate('Bulge', xy=(18, 0),  xycoords='data',xytext=None,size=20)
    axes[1].tick_params(axis='y', labelsize=20)
    axes[1].set_ylim(-3,9)
    
    Xiso,Yiso,Xirs,Yirs = plotvaluesIso()
    axes[2].plot(Xiso,Yiso,'-',linewidth=3, color='0.6')
    axes[2].plot(Xirs,Yirs,'-',linewidth=3, color='k')
    axes[2].annotate('Region 9', xy=(18, -4),  xycoords='data',xytext=None,size=20)
    axes[2].set_ylim(-5, 6.2)
    
    axes[2].set_xlabel("Wavelength ($\mu m$)",fontsize=24)
    plt.xticks([4,6,8,10,12,14,16,18,20,22])
    axes[2].tick_params(axis='x', labelsize=20)
    axes[2].tick_params(axis='y', labelsize=20)

    fig.text(0.055, 0.5, 'Intensity (MJy/sr)', ha='center', va='center', rotation='vertical',fontsize=24)
    
    plt.show()
    
fig,axes = plt.subplots(3,1,sharex=True, figsize=(4.5,5.75))
fig.subplots_adjust(hspace=0)

plotting()
