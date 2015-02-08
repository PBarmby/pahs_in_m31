import numpy as np                     
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import math

##### SECTION: calculate stuff

#TODO: finish (figure out correct way to do this..)
#TODO: apply (add to data tables)
def compute_rhi(atomic_lines, atomic_line_unc):
    """ Calculate single RHI value
    Inputs: atomic_lines: tuple with SIV,SIII,NeIII, NeII fluxes
            atomic_line_unc: tuple with SIV,SIII,NeIII, NeII flux uncertainties
            convention: val=num, unc=NaN means that val is an upper limit,
                        val=NaN, unc=NaN means 'nodata'.
                        val=0, unc=NaN also means 'nodata' 

    Output : tuple of (RHI value, uncertainty)."""

    if not any(math.isnan(element) for element in atomic_lines):

        RHI = ((np.log10(NeIII/NeII)) + ( 0.71 + (1.58*(np.log10(SIV/SIII)))))/2
        RHI_unc = delII

    # uncertainy propagation
    delNe = math.sqrt(((NeIIIunc/NeIII)**2) + ((NeIIunc/NeII)**2))*(NeIII/NeII)
    dellogNe = (delNe/(NeIII/NeII))*0.434
    delS = math.sqrt(((SIVunc/SIV)**2) + ((SIIIunc/SIII)**2))*(SIV/SIII)
    dellogS = (delS/(SIV/SIII))*0.434
    RHI_unc = math.sqrt((dellogS**2)+(dellogNe**2))
    
    if SIV == 0.000 :
        RHI = np.log10(NeIII/NeII)
        RHI_unc = dellogNe
        
    if NeII == 0.00000 or NeIII == 0.00000 :

        RHI = ( 0.71 + (1.58*(np.log10(SIV/SIII))))
        RHI_unc = dellogS
        
    if SIV != 0.00 and NeII != 000 and NeIII != 0.00 :

        #II = ((np.log10(NeIII/NeII))/(dellogNe**2) + ( 0.71 + (1.58*(np.log10(SIV/SIII))))/(dellogS**2))/(dellogS**(-2)+ dellogNe**(-2))

   
    return(RHI,RHI_unc)


# SECTION: make plots

def make_fig_10_plot(engel_tab, m31_tab):
    # plot Engelbracht data
    Y= engel_tab('PAH8')
    Yerr = engel_tab('PAH8_unc')    
    Yerr = (Yerr/Y)*0.434
    Y = np.log10(Y)
    X = engel_tab('RHI')
    Xerr = engel_tab('RHI_unc')
    plt.errorbar(X,Y,Yerr,Xerr,'ks',linewidth=2.0)
    plt.plot(X,Y,'ks',label = 'Engelbracht et al. 2008', markersize=15 , mfc = 'white')
       
    # plot M31 data
    X1,X1err,Y1,Y1err = m31tab['RHI'], m31tab['RHI_unc'], m31tab['PAH8'], m31tab['PAH8_unc'] # NB: PAH8 = 7.7+8.3+8.6, may need to add
    plt.plot(X1[1], Y1[1], 'b<', markersize=20,linewidth=2.0) #Upper Limit
    plt.errorbar(X1[3], Y1[3] ,Y1err[3],X1err[3], 'bo', markersize=15,linewidth=2.0)
    plt.plot(X1[4], Y1[4] , 'b>', markersize=20,linewidth=2.0)
    plt.errorbar(X1[5], Y1[5] ,Y1err[5],X1err[5], 'bo', markersize=15,linewidth=2.0)
    plt.errorbar(X1[6], Y1[6] ,Y1err[6],X1err[6], 'bo', markersize=15,linewidth=2.0)
    plt.plot(X1[7], Y1[7] , 'b<', markersize=20,linewidth=2.0) # Lowe limit
    plt.errorbar(X1[8], Y1[8] ,Y1err[8],X1err[8], 'bo', markersize=15,linewidth=2.0)
    plt.plot(X1[8], Y1[8] , 'bo', markersize=15,linewidth=2.0, label = 'M31')
    
    # plot formatting
    plt.xlabel("Log(RHI)", fontsize = 30)
    plt.ylabel("Log(EQW_8($\mu m$))", fontsize = 30)
    plt.legend( loc='lower left',prop={'size':20} )
    plt.yticks([-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0,1.2])
    ax = plt.subplot(111)
    minorLocator   = AutoMinorLocator(5)
    ax.xaxis.set_minor_locator( minorLocator)
    plt.tick_params(which='both', width=2)
    plt.tick_params(which='major', length=10)
    plt.tick_params(which='minor', length=7, color='k')
    plt.show()
    return

def make_figure_11(m31dat, gord_dat):    
    fig,axes = plt.subplots(2,1,sharex=True)
    
    plotting(m31dat, gord_dat, 'PAH7.7',0,"Log(EQW_7.7($\mu m$))")
    plotting(m31dat, gord_dat, 'PAH11.3',1,"Log(EQW_11.3($\mu m$))")
    return


def plotting(m31_tab, gordon_tab, feature,panel,Ylabel):
    """ Plotting EQWs vs II for M31 regions and also over plotting it with Gordon et al. 2008.

    X values : II values
    Y values : EQWs    """ 

    X1,X1err,Y1,Y1err = m31_tab['RHI'], m31_tab['RHI_unc'], m31_tab[feature], m31_tab[feature+'_unc']  # M31 data
    X2,X2err,Y2,Y2err = gordon_tab['RHI'], gordon_tab['RHI_unc'], gordon_tab[feature], gordon_tab[feature+'_unc']  # Gordon's data

    # this is a mess, clean it up
    axes[panel].errorbar(X2,Y2,Y2err,X2err,'ks', markersize=15,linewidth=2.0)
    axes[panel].plot(X2,Y2,'gs', markersize=15,linewidth=2.0, mfc = "white", label = 'Gordon et al. 2008')

    axes[panel].plot(X1[1], Y1[1], 'b<', markersize=20,linewidth=2.0) #Upper Limit
    axes[panel].errorbar(X1[3], Y1[3] ,Y1err[3],X1err[3], 'bo', markersize=15,linewidth=2.0)
    axes[panel].plot(X1[4], Y1[4], 'b>', markersize=20,linewidth=2.0)
    axes[panel].errorbar(X1[5], Y1[5] ,Y1err[5],X1err[5], 'bo', markersize=15,linewidth=2.0)
    axes[panel].errorbar(X1[6], Y1[6] ,Y1err[6],X1err[6], 'bo', markersize=15,linewidth=2.0)
    axes[panel].plot(X1[7], Y1[7] ,'b<', markersize=20,linewidth=2.0)
    axes[panel].errorbar(X1[8], Y1[8] ,Y1err[8],X1err[8], 'bo', markersize=15,linewidth=2.0)
    axes[panel].plot(X1[8], Y1[8] , 'bo', markersize=15,linewidth=2.0, label = 'M31')

    axes[panel].set_ylabel(Ylabel,fontsize=24)
    axes[1].set_xlabel("Log(RHI)",fontsize=24)
    axes[panel].tick_params(axis='both', which='major', labelsize=30)
    axes[panel].tick_params(axis='both', which='minor', labelsize=30)
    minorLocator   = AutoMinorLocator(5)
    axes[1].xaxis.set_minor_locator(minorLocator)
    axes[0].xaxis.set_minor_locator(minorLocator)
    axes[1].tick_params(which='both', width=2)
    axes[0].tick_params(which='both', width=2)
    axes[1].tick_params(which='major', length=10)
    axes[0].tick_params(which='major', length=10)
    axes[1].tick_params(which='minor', length=7, color='k')
    axes[0].tick_params(which='minor', length=7, color='k')
    

    axes[panel].legend( loc='lower left',prop={'size':20} )

    plt.show()
    return

  
def make_figure_12(engel_tab, m31_tab, feature_list):
#   plot formatting
    fig, ax = plt.subplots()
    minorLocator   = MultipleLocator(5)
    ax.xaxis.set_minor_locator(MultipleLocator(4))
    ax.yaxis.set_minor_locator(minorLocator)
    
    plt.tick_params(which='both', width=2)
    plt.tick_params(which='major', length=10)
    plt.tick_params(which='minor', length=7, color='k')
    ax.xaxis.set_minor_locator(MultipleLocator(5))

    # get data
    Oxy = engel_tab['12plogOH']
    Oxyunc = engel_tab['12plogOH_unc']
    
    # Removing IRC3 data. (TODO: still need to do this)
    X = m31_tab['12plogOH']
    Xerr = m31_tab['12plogOH_unc']

    # loop ovr features to be plotted
    for feat in feature_list:
        Y = m31_tab[feat]
        Yerr = m31_tab[feat+_'unc']
        plt.errorbar(X,Xerr, Yerr, Xerr,symbol,color = col, linewidth=2.0)
        plt.plot(X-0.35, Xerr,symbol,color = col, markersize=15,label = featureL)

        EQW= engel_tab[feat]
        EQW_unc = enget_tab[feat+'_unc']
        plt.errorbar(Oxy,EQW,EQWunc,Oxyunc,'o',color = '0.75', linewidth=2.0)
        plt.plot(Oxy,EQW,'o',mfc = 'white', markersize=15)

    plt.xlabel("12+ log[O/H] " ,fontsize=28)
    plt.ylabel("EQW ($\mu m$)" ,fontsize=28)
    plt.legend( loc='lower left' ,prop={'size':20} )
    plt.yscale('log')
    plt.show()

    return
	
