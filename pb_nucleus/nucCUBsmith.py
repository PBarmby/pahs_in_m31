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
import sys
import scipy 
from matplotlib.ticker import FormatStrFormatter,ScalarFormatter
from matplotlib.ticker import LogFormatterExponent,LogLocator
from matplotlib import ticker
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, \
     AnnotationBbox
from matplotlib.cbook import get_sample_data
from matplotlib._png import read_png

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

def doplot(color=True):
    fig, ax = plt.subplots()
    
    ax = plt.subplot(111)       
    
    ax.set_xlabel(" Wavelength ($\mu m$)", fontsize=20)
    ax.set_ylabel("Flux (Arbitrary Units)", fontsize=20)
    
    plt.tick_params(which='both', width=2)
    plt.tick_params(which='major', length=10)
    plt.tick_params(which='minor', length=7, color='k')
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
    
        plt.plot(X,Y,'-', label = 'M31', linewidth = 4)
        
        spec = np.loadtxt("ngc1316.dat", skiprows = 2 )
        Wave,flux,unc  = spec[:,0] , spec[:,1], spec[:,2]
        
        unc = unc/Wave
        plt.plot(Wave, ((12*flux/Wave)),'-', label = 'NGC1316', linewidth = 1.5) # 12 is wrong correct no is 3
        
        spec = np.loadtxt("ngc4594.dat", skiprows = 2 )
        Wave,flux,unc  = spec[:,0] , spec[:,1], spec[:,2]
        
        unc = unc/Wave
        plt.plot(Wave, ((14*flux/Wave)),'-', label = 'NGC4594', linewidth = 1.5) # 12 is wrong correct no is 3
        
        
        spec = np.loadtxt("ngc1404.dat", skiprows = 2 )
        Wave,flux,unc  = spec[:,0] , spec[:,1], spec[:,2]
        
        unc = unc/Wave
        plt.plot(Wave, ((14*flux/Wave)),'-', label = 'NGC1404', linewidth = 1.5) # 12 is wrong correct no is 3
        
        spec = np.loadtxt("ngc2841.dat", skiprows = 2 )
        Wave,flux,unc  = spec[:,0] , spec[:,1], spec[:,2]
        
        unc = unc/Wave
        plt.plot(Wave, ((10*flux/Wave)),'-', label = 'NGC2841', linewidth = 1.5) # 12 is wrong correct no is 3
        
        spec = np.loadtxt("ngc4552.dat", skiprows = 2 )
        Wave,flux,unc  = spec[:,0] , spec[:,1], spec[:,2]
        
        unc = unc/Wave
        plt.plot(Wave, ((10*flux/Wave)),'-', label = 'NGC4552', linewidth = 1.5) # 12 is wrong correct no is 3
        
        
        spec = np.loadtxt("ngc4125.dat", skiprows = 2 )
        Wave,flux,unc  = spec[:,0] , spec[:,1], spec[:,2]
        
        unc = unc/Wave
        plt.plot(Wave, ((8*flux/Wave)),'-', label = 'NGC4125', linewidth = 1.5) # 12 is wrong correct no is 3
    
    else:
    ######################################################################
    #Black and white plotting

        plt.plot(X,Y,'-', label = 'M31', linewidth = 4)
        
        spec = np.loadtxt("ngc1316.dat", skiprows = 2 )
        Wave,flux,unc  = spec[:,0] , spec[:,1], spec[:,2]
        
        unc = unc/Wave
        plt.plot(Wave, ((12*flux/Wave)),'-', label = 'NGC1316',color = '0.10', linewidth = 1.5) # 12 is wrong correct no is 3
        
        spec = np.loadtxt("ngc4594.dat", skiprows = 2 )
        Wave,flux,unc  = spec[:,0] , spec[:,1], spec[:,2]
        
        unc = unc/Wave
        plt.plot(Wave, ((14*flux/Wave)),'-', label = 'NGC4594',color = '0.20') # 12 is wrong correct no is 3
        
        
        spec = np.loadtxt("ngc1404.dat", skiprows = 2 )
        Wave,flux,unc  = spec[:,0] , spec[:,1], spec[:,2]
        
        unc = unc/Wave
        plt.plot(Wave, ((14*flux/Wave)),'-', label = 'NGC1404',color = '0.40') # 12 is wrong correct no is 3
        
        spec = np.loadtxt("ngc2841.dat", skiprows = 2 )
        Wave,flux,unc  = spec[:,0] , spec[:,1], spec[:,2]
        
        unc = unc/Wave
        plt.plot(Wave, ((10*flux/Wave)),'-', label = 'NGC2841',color = '0.50') # 12 is wrong correct no is 3
        
        spec = np.loadtxt("ngc4552.dat", skiprows = 2 )
        Wave,flux,unc  = spec[:,0] , spec[:,1], spec[:,2]
        
        unc = unc/Wave
        plt.plot(Wave, ((10*flux/Wave)),'-', label = 'NGC4552',color = '0.70') # 12 is wrong correct no is 3
        
        
        spec = np.loadtxt("ngc4125.dat", skiprows = 2 )
        Wave,flux,unc  = spec[:,0] , spec[:,1], spec[:,2]
        
        unc = unc/Wave
        plt.plot(Wave, ((8*flux/Wave)),'-', label = 'NGC4125',color = '0.80') # 12 is wrong correct no is 3
    
    
    ########################################################################################
    # This section inserts the subplot from silicate.py in to the main plot.
#    fn = get_sample_dat("./silicate.png", asfileobj=False)
    arr_lena = read_png("./silicate.png")
    
    imagebox = OffsetImage(arr_lena, zoom=0.25)
    
    ab = AnnotationBbox(imagebox, xy=(17, 20),
                            xycoords='data',
                            boxcoords="offset points")
                            
    ax.add_artist(ab)
    
#   legend, formatting
    
    ax.set_xlim(5.1, 37)
    ax.set_ylim(0, 31)
    ax.xaxis.set_major_formatter(ScalarFormatter())
    plt.xscale('log',subsx=[6,7,8,9,10,11,12,15,20,25,30])
    ax.xaxis.set_minor_formatter(FormatStrFormatter('%1.0f'))
    ax.xaxis.set_major_formatter(ScalarFormatter())
    
    plt.tick_params(\
        axis='x',          # changes apply to the x-axis
        which='major',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='on',         # ticks along the top edge are off
        labelbottom='off') # labels along the bottom edge are off

    ax.legend(loc='upper right',prop={'size':20})
    
#    plt.draw()
    plt.show()
    return

