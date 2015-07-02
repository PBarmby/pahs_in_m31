"""
Purpose : Plots line intensities 7.7/11.3 vs 6.2/11.3 PAH features for M31 and SINGS survey. This also

Inputs  : 1) Two data files that has Line intensities and their uncertainties for 11 PAH features.
          2) Table 4 from Smith at al. 2007

Outputs : 1) Plot of 7.7/11.3 vs 6.2/11.3 PAH features for M31 and SINGS survey

Notes  : 

"""

import numpy as np                     
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
from astropy.table import Table
from matplotlib.ticker import AutoMinorLocator

## section: defaults for plotting
ms = 15 # markersize
label_font_size = 15
legend_size = 12

def make_plot():

    mpl.rcParams['font.weight']='normal'        
    mpl.rcParams['axes.labelweight']='normal'        

    Sings_PAH= np.loadtxt("SINGS_PAH_val.txt")
    feature_plots(Sings_PAH[0,:],Sings_PAH[2,:],Sings_PAH[8,:],Sings_PAH[1,:],Sings_PAH[3,:],Sings_PAH[9,:],2)
    
    m31_dat = Table.read('m31_alldat.fits')
    feature_plots(m31_dat['PAH6.2flx'],m31_dat['PAH7.7flx'],m31_dat['PAH11.3flx'],\
       m31_dat['PAH6.2flx_unc'],m31_dat['PAH7.7flx_unc'],m31_dat['PAH11.3flx_unc'],1)    
    plt.show()
    return

def feature_plots(PAH_val6_2,PAH_val7_7,PAH_val11_3,PAH_unc6_2,PAH_unc7_7,PAH_unc11_3,data) :
    #Plotting  7.7/11.3 vs 6.2/11.3. key word "data" is an integer which defines which data set you want to plot. 1 = M31, 2 = SINGS
    
    xaxis = []
    yaxis = []
    x_err = []
    y_err = []
    for i in range(np.shape(PAH_val6_2)[0]):
        if PAH_val6_2[i] != 0 and PAH_val7_7[i] != 0 and PAH_val11_3[i] != 0 :
            
            xaxis.append(PAH_val6_2[i]/ PAH_val11_3[i])
            yaxis.append(PAH_val7_7[i]/ PAH_val11_3[i])
            y_err.append(math.sqrt(((PAH_unc7_7[i]/PAH_val11_3[i])**2) + ((PAH_unc11_3[i]/PAH_val7_7[i])**2)))
            x_err.append(math.sqrt(((PAH_unc6_2[i]/PAH_val11_3[i])**2) + ((PAH_unc11_3[i]/PAH_val6_2[i])**2)))
           
    
    y_err = (np.array(y_err)/np.array(yaxis))*0.434
    x_err = (np.array(x_err)/np.array(xaxis))*0.434
    W = 1/(y_err)

    X = [np.log10(a) for a in xaxis]  #use this for logscale
    Y = [np.log10(a) for a in yaxis]  #use this for logscale

    if data == 1:
     plt.errorbar(X,Y,y_err,x_err,'bo', markersize = ms, mfc = 'k', label = 'M31', color = '0.75')

    else:
     plt.errorbar(X,Y,y_err,x_err,'s', color = '0.75', markersize = ms, mfc = 'white', label = 'Smith at al. 2007')


    plt.ylabel("Log(PAH_7.7/PAH_11.2)", fontsize = label_font_size)
    plt.xlabel("Log(PAH_6.2/PAH_11.2)", fontsize = label_font_size)
    plt.legend(loc='lower right',prop={'size':legend_size})
    minorLocator   = AutoMinorLocator(5)
    ax = plt.subplot(111)
    ax.xaxis.set_minor_locator(minorLocator)
    #ax.yaxis.set_minor_locator(minorLocator)

    plt.tick_params(which='both', width=2)
    plt.tick_params(which='major', length=10)
    plt.tick_params(which='minor', length=7, color='k')

    return
