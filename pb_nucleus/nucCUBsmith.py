"""
Purpose : Plots the spectum from the nucleus of M31 with that of 6 other galaxies.
          Also this insert a subplot inside this plot. This  sub plot shows the the
          plot obtained from silicate.py

Inputs  : 1) .tbl files of the nucleus obtained from CUBISM
          2) spectra data for 6 galaxis (Smith et al. 2007)
          3) plot generated by the silicate.py as a png file

Notes  : The 1st section is similar to nucCUBISM.py

"""


import numpy as np                     # Allows Numpy functions to be called directly
import matplotlib.pyplot as plt         # Graphing routines
import matplotlib as mpl
import sys
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
from matplotlib.ticker import FormatStrFormatter,ScalarFormatter
from matplotlib.ticker import LogFormatterExponent,LogLocator
from matplotlib import ticker
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, \
     AnnotationBbox
from matplotlib._png import read_png
from mpl_toolkits.axes_grid.inset_locator import inset_axes

# original values
#spec_mult_offsets = {'M31':1.0, 'NGC1316':12.0, 'NGC4594':14.0, 'NGC1404':14.0, 'NGC2841':10.0, 'NGC4552':10.0, 'NGC4125':8.0}
# worked out by trial and error
spec_mult_offsets = {'M31':0.25, 'NGC1316':3.3, 'NGC4594':3.5, 'NGC1404':3.5, 'NGC2841':2.0, 'NGC4552':1.5, 'NGC4125':1.0}

axlabelsize=20

def doplot(color=True):
    fig, ax = plt.subplots(figsize=(6,3))
    mpl.rcParams['font.weight']='normal'        
    mpl.rcParams['axes.labelweight']='normal'    
    
    ax = plt.subplot(111)       
    
    ax.set_xlabel("Wavelength ($\mu$m)", fontsize=axlabelsize)
    ax.set_ylabel(r"Scaled intensity $I_{\nu}$ (arbitrary units)", fontsize=axlabelsize)
    
    plt.tick_params(which='both', width=2)
    plt.tick_params(which='major', length=10)
    plt.tick_params(which='minor', length=7, color='k')
    plt.tick_params(which='both', labelsize=axlabelsize*0.8)
#    plt.show()
    
    # I changed this code from others to match it with Smith's spectra
    
    X, Y = np.loadtxt("nucFLUX", usecols=(0,1), unpack=True)
    Yerr = np.loadtxt("nucUNC")
#    Yerr = Yerr[:,1]
    
    Y = (np.array(Y))*3/(np.array(X))  # changed from Intensity to FLux. There should be a 10^-6 term.
                                       # This doesn't matter because smith's spectra also has that term out.
    
    #####################################################################
    #Plotting with colors.

    if color:
    
        plt.plot(X,Y*spec_mult_offsets['M31'],'-', label = 'M31', linewidth = 4)
        
        spec = np.loadtxt("ngc1316.dat", skiprows = 2 )
        Wave,flux,unc  = spec[:,0] , spec[:,1], spec[:,2]
        unc = unc/Wave
        plt.plot(Wave, ((spec_mult_offsets['NGC1316']*flux/Wave)),'-', label = 'NGC 1316', linewidth = 1.5) # 12 is wrong correct no is 3
        
        spec = np.loadtxt("ngc4594.dat", skiprows = 2 )
        Wave,flux,unc  = spec[:,0] , spec[:,1], spec[:,2]
        unc = unc/Wave
        plt.plot(Wave, ((spec_mult_offsets['NGC4594']*flux/Wave)),'-', label = 'NGC 4594', linewidth = 1.5) # 12 is wrong correct no is 3
        
        
        spec = np.loadtxt("ngc1404.dat", skiprows = 2 )
        Wave,flux,unc  = spec[:,0] , spec[:,1], spec[:,2]
        unc = unc/Wave
        plt.plot(Wave, ((spec_mult_offsets['NGC1404']*flux/Wave)),'-', label = 'NGC 1404', linewidth = 1.5) # 12 is wrong correct no is 3
        
        spec = np.loadtxt("ngc2841.dat", skiprows = 2 )
        Wave,flux,unc  = spec[:,0] , spec[:,1], spec[:,2]
        unc = unc/Wave
        plt.plot(Wave, ((spec_mult_offsets['NGC2841']*flux/Wave)),'-', label = 'NGC 2841', linewidth = 1.5) # 12 is wrong correct no is 3
        
        spec = np.loadtxt("ngc4552.dat", skiprows = 2 )
        Wave,flux,unc  = spec[:,0] , spec[:,1], spec[:,2]
        unc = unc/Wave
        plt.plot(Wave, ((spec_mult_offsets['NGC4552']*flux/Wave)),'-', label = 'NGC 4552', linewidth = 1.5) # 12 is wrong correct no is 3
        
        
        spec = np.loadtxt("ngc4125.dat", skiprows = 2 )
        Wave,flux,unc  = spec[:,0] , spec[:,1], spec[:,2]
        unc = unc/Wave
        plt.plot(Wave, ((spec_mult_offsets['NGC4125']*flux/Wave)),'-', label = 'NGC 4125', linewidth = 1.5) # 12 is wrong correct no is 3
    
    else:
    ######################################################################
    #Black and white plotting

        plt.plot(X,Y*spec_mult_offsets['M31'],'-', label = 'M31', linewidth = 4, color='k')
        
        spec = np.loadtxt("ngc1316.dat", skiprows = 2 )
        Wave,flux,unc  = spec[:,0] , spec[:,1], spec[:,2]
        unc = unc/Wave
        plt.plot(Wave, ((spec_mult_offsets['NGC1316']*flux/Wave)),'-', label = 'NGC 1316',color = '0.10', linewidth = 1.5) 
        
        spec = np.loadtxt("ngc4594.dat", skiprows = 2 )
        Wave,flux,unc  = spec[:,0] , spec[:,1], spec[:,2]
        unc = unc/Wave
        plt.plot(Wave, ((spec_mult_offsets['NGC4594']*flux/Wave)),'-', label = 'NGC 4594',color = '0.20', ls='dotted') 
        
        
        spec = np.loadtxt("ngc1404.dat", skiprows = 2 )
        Wave,flux,unc  = spec[:,0] , spec[:,1], spec[:,2]
        unc = unc/Wave
        plt.plot(Wave, ((spec_mult_offsets['NGC1404']*flux/Wave)),'-', label = 'NGC 1404',color = '0.40') 
        
        spec = np.loadtxt("ngc2841.dat", skiprows = 2 )
        Wave,flux,unc  = spec[:,0] , spec[:,1], spec[:,2]
        unc = unc/Wave
        plt.plot(Wave, ((spec_mult_offsets['NGC2841']*flux/Wave)),'-', label = 'NGC 2841',color = '0.50', ls='dashed') 
        
        spec = np.loadtxt("ngc4552.dat", skiprows = 2 )
        Wave,flux,unc  = spec[:,0] , spec[:,1], spec[:,2]
        unc = unc/Wave
        plt.plot(Wave, ((spec_mult_offsets['NGC4552']*flux/Wave)),'-', label = 'NGC 4552',color = '0.70') 
        
        spec = np.loadtxt("ngc4125.dat", skiprows = 2 )
        Wave,flux,unc  = spec[:,0] , spec[:,1], spec[:,2]
        unc = unc/Wave
        plt.plot(Wave, ((spec_mult_offsets['NGC4125']*flux/Wave)),'-', label = 'NGC 4125',color = '0.05') 
    
    
    ########################################################################################
    # This section inserts the north and nucleus spectra into the main plot.

    # Reading .tbl files got from CUBISM 
    sl2 = np.loadtxt("m31nuc_sl2_nucCentre.tbl", skiprows = 15 )
    sl1 = np.loadtxt("m31nuc_sl1_nucCentre.tbl", skiprows = 15 )
    ll2 = np.loadtxt("m31nuc_ll2_nucCentre.tbl", skiprows = 15 )
    Wavelength,Flux,FluxUnc = getspectrum(sl1,sl2,ll2)
    
    inax = inset_axes(ax, width = "33%" , height = "20%", loc=(10,4.5))
    plotting(Wavelength,Flux,FluxUnc,inax,col=color)
    inax.set_xlabel("Wavelength ($\mu$m)", fontsize=axlabelsize)
    inax.annotate('Silicate', xy=(9.7, 50), xycoords='data',xytext=(0.6,0.7),textcoords='axes fraction',size=axlabelsize*0.8 ,arrowprops=dict(arrowstyle="->", linewidth = 3))
    inax.annotate('Nucleus', xy=(5,20),xycoords='data',size=axlabelsize*0.8)
    inax.set_ylim(10,88)
    inax.set_xlim(4,22)
    inax.yaxis.set_major_locator(MultipleLocator(20))    
    inax.xaxis.set_major_locator(MultipleLocator(5))

    Wavelength,Flux,FluxUnc = np.loadtxt('m31nuc_nucUP_correct.dat',usecols=(0,1,2),unpack=True)
    inax2 = inset_axes(ax, width = "33%" , height = "20%", loc=9)
    
    plotting(Wavelength,Flux,FluxUnc,inax2,col=color)
    inax2.set_ylim(9,45)
    inax2.set_xlim(4,22)
    inax2.annotate('North', xy=(5,12),xycoords='data',size=axlabelsize*0.8)
    inax2.yaxis.set_major_locator(MultipleLocator(10))
    inax2.xaxis.set_major_locator(MultipleLocator(5))

#    minorLocator = AutoMinorLocator(5)
 #   for ax in [inax,inax2]:
 #       ax.xaxis.set_minor_locator(minorLocator)
 #       ax.tick_params(which='both', width=2)
 #       ax.tick_params(which='major', length=10)
 #       ax.tick_params(which='minor', length=7, color='k')
 #       ax.tick_params(which='both',labelsize=axlabelsize*0.7)
 #       ax.tick_params(axis='both', which='major', labelsize=axlabelsize*0.5)
 #       ax.tick_params(axis='both', which='minor', labelsize=axlabelsize*0.5)
#    inax.text(0.03, 0.5, 'Intensity (MJy/sr)', ha='center', va='center', rotation='vertical',fontsize=axlabelsize)
    inax.tick_params(labelbottom='off')
    
#   legend, formatting
    
    ax.set_xlim(5.1, 37)
    ax.set_ylim(0, 9)
    ax.set_xscale('log',subsx=[6,7,8,9,10,11,12,15,20,25,30])
    ax.xaxis.set_minor_formatter(FormatStrFormatter('%1.0f'))
    ax.xaxis.set_major_formatter(ScalarFormatter())
    
    ax.tick_params(\
        axis='x',          # changes apply to the x-axis
        which='major',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='on',         # ticks along the top edge are on
        labelbottom='off') # labels along the bottom edge are off

    ax.legend(loc='upper right',prop={'size':13})
    
#    plt.draw()
    plt.show()
    return

""" This code plots spectra extracted from the crntre of the nucleus and the north region of M31
in two panels.
You should have .tbl files obtained from CUBISM for both regions to run this code."""

def stitch_nucleus():
    sl2 = np.loadtxt("m31nuc_sl_sl2_cube__fits__SL2_sp.tbl", skiprows = 15 )
    W_SL2,F_SL2  = sl2[:,0] , sl2[:,1]
    
    sl3 = np.loadtxt("m31nuc_sl_sl3_cube__fits__SL2_sp.tbl", skiprows = 15 )
    W_SL3,F_SL3  = sl3[:,0] , sl3[:,1]
    
    sl1 = np.loadtxt("m31nuc_sl_sl1_cube__fits__SL1_sp.tbl", skiprows = 15 )
    W_SL1,F_SL1  = sl1[:,0] , sl1[:,1]
    
    ll2 = np.loadtxt("m31nuc_ll_ll2_cube__fits__LL2_sp.tbl", skiprows = 15 )
    W_LL2,F_LL2  = ll2[:,0] , ll2[:,1]
    
    
    Xval = W_SL1[(W_SL1 >= 14.33) & (W_SL1 <= 14.68)]
    XLL2 = W_LL2[(W_LL2 >= 14.33) & (W_LL2 <= 14.68)]
    YLL2 = F_LL2[(W_LL2 >= 14.33) & (W_LL2 <= 14.68)]
    YSL1 = F_SL1[(W_SL1 >= 14.33) & (W_SL1 <= 14.68)]
    Yval = np.interp(Xval,XLL2, YLL2)
    
    fractions3 = []
    for YSL1,Yval in zip(YSL1,Yval):
        ratio = Yval - YSL1
        fractions3.append(ratio)
    
    factor3 = np.mean(fractions3)
    
    Y = []
    X = []
    
    
    for (i,aval) in enumerate(W_SL2):
        if aval < 7.57 :
            Y.append(F_SL2[i])
            X.append(aval)
    
    #for (i,aval) in enumerate(W_SL3):
       # if (aval > 7.55) & (aval < 7.67) :
          #  Y.append(F_SL3[i])
          #  X.append(aval)
    
    for (i,aval) in enumerate(W_SL1):
        if (aval > 7.57) & (aval < 14.50):
            Y.append(F_SL1[i])
            X.append(aval)
    
    for (i,aval) in enumerate(W_LL2):
        if aval > 14.50 :
            Y.append(F_LL2[i])
            X.append(aval)
    
    with open('nucFLUX', 'w') as f:
        for x,y in zip(X,Y):
            f.write(str(x))
            f.write('\t'+str(y)+'\n')
    f.close()
    return

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


def plotting(X,Y,Yerr,ax,col=True):

    # Plot spectra into panels.
    if col:
        ax.errorbar(X,Y,Yerr,0,'-')
    else:
        ax.errorbar(X,Y,Yerr,0,'-', color='k')
    ax.tick_params(axis='both', which='major', labelsize=24)
    ax.tick_params(axis='both', which='minor', labelsize=24)
    return

