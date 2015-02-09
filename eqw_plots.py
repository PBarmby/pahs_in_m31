import numpy as np                     
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import math, string
from astropy.table import Table, Column, join
import make_tab

## section: defaults for plotting
ms = 20 # markersize
label_font_size = 20
legend_size = 20

##### SECTION: calculate stuff

#TODO: finish figuring out dealing with limits
def compute_rhi(atomic_lines, atomic_lines_unc):
    """ Calculate single RHI value
    Inputs: atomic_lines: tuple with SIV,SIII,NeIII, NeII fluxes
            atomic_lines_unc: tuple with SIV,SIII,NeIII, NeII flux uncertainties
            convention: val=num, unc=NaN means that val is an upper limit,
                        val=NaN, unc=NaN means 'nodata'.
                        val=0, unc=NaN also means 'nodata' 

    Output : tuple of (RHI value, uncertainty)."""

    SIV, SIII, NeIII, NeII = atomic_lines
    SIV_unc, SIII_unc, NeIII_unc, NeII_unc = atomic_lines_unc

    if math.isnan(NeIII_unc) or math.isnan(NeII_unc):
        term1 =  0.71 + 1.58*(np.log10(SIV/SIII))
        RHI_unc = float('NaN') # possibly not completely correct
    else:
        term1 = np.log10(NeIII/NeII)

    if math.isnan(SIV_unc) or math.isnan(SIII_unc):
        term2 =  (np.log10(NeIII/NeII) - 0.71)/1.58
        RHI_unc = float('NaN')
    else:
        term2 = 0.71 + 1.58*(np.log10(SIV/SIII))
    
    RHI = term1+term2

    if not any(math.isnan(element) for element in atomic_lines_unc):
        RHI = ((np.log10(NeIII/NeII)) + ( 0.71 + (1.58*(np.log10(SIV/SIII)))))/2
	# uncertainy propagation
	delNe = math.sqrt(((NeIII_unc/NeIII)**2) + ((NeII_unc/NeII)**2))*(NeIII/NeII)
	dellogNe = (delNe/(NeIII/NeII))*0.434
	delS = math.sqrt(((SIV_unc/SIV)**2) + ((SIII_unc/SIII)**2))*(SIV/SIII)
	dellogS = (delS/(SIV/SIII))*0.434
	RHI_unc = math.sqrt((dellogS**2)+(dellogNe**2))
	    
    return(RHI,RHI_unc)

def add_rhi(tab_file, tabformat='ascii.commented_header'):
    tab = Table.read(tab_file,format=tabformat)
    if 'RHI' not in tab.colnames:
        nrows = len(tab)
        rhi = np.zeros(nrows)
        rhi_unc = np.zeros(nrows)
        for row in range(0,nrows):
            atlines = (tab['SIV'][row],tab['SIII'][row],tab['NeIII'][row],tab['NeII'][row])
            atlines_unc = (tab['SIV_unc'][row],tab['SIII_unc'][row],tab['NeIII_unc'][row],tab['NeII_unc'][row])
            rhi[row],rhi_unc[row] = compute_rhi(atlines,atlines_unc)
#            print atlines, atlines_unc, rhi[row], rhi_unc[row]
        tab['RHI'] = rhi
        tab['RHI_unc'] = rhi_unc
    else:
        print 'Table already has RHI'
    return(tab)

def process_comparison_data():
    # add RHI to Gordon data
    grd_rhi = add_rhi('gordon_atomic.dat')
    grd_eqw = Table.read('GordonEQW',format='ascii.commented_header')
    # join atomic and EQW files
    gordon_full = join(grd_eqw, grd_rhi,keys='ID')
    gordon_full.write('gordon_m101.dat', format='ascii.commented_header')

    # add RHI to Engelbracht data
    eng_in = add_rhi('englbrt.dat')
    eng_met = Table.read('englbrt_eqw_oxy',format='ascii.commented_header')
    eng_out = join(eng_in, eng_met, keys=['ID','PAH8eqw','PAH8eqw_unc'])
    eng_out.write('englbrt_sb.dat',format='ascii.commented_header')
    return

def process_m31_data():
    outtab = add_rhi('m31mega.fits',tabformat='fits')
    outtab['PAH8eqw'],outtab['PAH8eqw_unc'] = make_tab.compute_complex(outtab,['PAH7.7eqw','PAH8.3eqw','PAH8.6eqw'])
    outtab.write('m31_alldat.fits',format='fits')
    return


# SECTION: make plots
def make_fig_10_plot(engel_tab, m31_tab):
    fix, ax = plt.subplots()
    # plot Engelbracht data
    Y= engel_tab['PAH8eqw']
    Yerr = engel_tab['PAH8eqw_unc']    
    Yerr = (Yerr/Y)*0.434
    Y = np.log10(Y)
    X = engel_tab['RHI']
    Xerr = engel_tab['RHI_unc']
    ax.errorbar(X,Y,Yerr,Xerr,'ks',linewidth=2.0)
    ax.plot(X,Y,'ks',label = 'E08 starburst', markersize=ms*0.75 , mfc = 'white')
       
    # plot M31 data
    X1,X1err,Y1,Y1err = m31_tab['RHI'], m31_tab['RHI_unc'], np.log10(m31_tab['PAH8eqw']), 0.434*(m31_tab['PAH8eqw_unc']/m31_tab['PAH8eqw'])
    ax.plot(X1,Y1, marker='o',markersize=ms*0.75, color='b', label='M31')
#    ax.errorbar(X1,Y1,Y1err,X1err, color='b', marker=None,markersize=ms*0.75,linewidth=2.0)
#    ax.plot(X1[1], Y1[1], 'b<', markersize=ms,linewidth=2.0) #Upper Limit
#    ax.errorbar(X1[3], Y1[3] ,Y1err[3],X1err[3], 'bo', markersize=ms*0.75,linewidth=2.0)
#    ax.plot(X1[4], Y1[4] , 'b>', markersize=ms,linewidth=2.0)
#    ax.errorbar(X1[5], Y1[5] ,Y1err[5],X1err[5], 'bo', markersize=ms*0.75,linewidth=2.0)
#    ax.errorbar(X1[6], Y1[6] ,Y1err[6],X1err[6], 'bo', markersize=ms*0.75,linewidth=2.0)
#    ax.plot(X1[7], Y1[7] , 'b<', markersize=ms,linewidth=2.0) # Lowe limit
#    ax.errorbar(X1[8], Y1[8] ,Y1err[8],X1err[8], 'bo', markersize=ms*0.75,linewidth=2.0)
#    ax.plot(X1[8], Y1[8] , 'bo', markersize=ms*0.75,linewidth=2.0, label = 'M31')
    
    # plot formatting
    ax.set_xlabel("log(RHI)", fontsize = label_font_size)
    ax.set_ylabel(r'log(8$\mathregular{\mu m}$ EQW)', fontsize = label_font_size)
    ax.legend(loc='best',prop={'size':legend_size} )
#    ax.yticks([-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0,1.2])
#    minorLocator   = AutoMinorLocator(5)
#    ax.xaxis.set_minor_locator( minorLocator)
#    ax.tick_params(which='both', width=2)
#    ax.tick_params(which='major', length=10)
#    ax.tick_params(which='minor', length=7, color='k')
    plt.draw()
    return

def make_figure_11(m31dat, gord_dat):    
    fig,axes = plt.subplots(2,1,sharex=True)
    
    plotting(m31dat, gord_dat, 'PAH7.7eqw',axes[0],r'log(7.7$\mathregular{\mu m}$ EQW)')
    plotting(m31dat, gord_dat, 'PAH11.3eqw',axes[1],r'log(11.3$\mathregular{\mu m}$ EQW)')
    axes[1].set_xlabel("Log(RHI)",fontsize = label_font_size)
    return


def plotting(m31_tab, gordon_tab, feature, ax, Ylabel):
    """ Plotting EQWs vs II for M31 regions and also over plotting it with Gordon et al. 2008.

    X values : II values
    Y values : EQWs    """ 

    X1,X1err,Y1,Y1err = m31_tab['RHI'], m31_tab['RHI_unc'], np.log10(m31_tab[feature]), 0.434*(m31_tab[feature+'_unc']/m31_tab[feature])  # M31 data
    X2,X2err,Y2,Y2err = gordon_tab['RHI'], gordon_tab['RHI_unc'], np.log10(gordon_tab[feature]),0.434*(gordon_tab[feature+'_unc']/gordon_tab[feature]) # Gordon's data

    ax.errorbar(X2,Y2,Y2err,X2err,'ks', markersize=ms*0.75,linewidth=2.0)
    ax.plot(X2,Y2,'gs', markersize=ms*0.75,linewidth=2.0, mfc = "white", label = 'G08 M101')

    ax.plot(X1,Y1, marker='o',markersize=ms*0.75, color='b', label='M31')
    ax.errorbar(X1, Y1,yerr=Y1err, marker=None, linewidth=2.0)

#    ax.plot(X1[1], Y1[1], 'b<', markersize=ms,linewidth=2.0) #Upper Limit
#    ax.errorbar(X1[3], Y1[3] ,Y1err[3],X1err[3], 'bo', markersize=ms*0.75,linewidth=2.0)
#    ax.plot(X1[4], Y1[4], 'b>', markersize=ms,linewidth=2.0)
#    ax.errorbar(X1[5], Y1[5] ,Y1err[5],X1err[5], 'bo', markersize=ms*0.75,linewidth=2.0)
#    ax.errorbar(X1[6], Y1[6] ,Y1err[6],X1err[6], 'bo', markersize=ms*0.75,linewidth=2.0)
#    ax.plot(X1[7], Y1[7] ,'b<', markersize=ms,linewidth=2.0)
#    ax.errorbar(X1[8], Y1[8] ,Y1err[8],X1err[8], 'bo', markersize=ms*0.75,linewidth=2.0)
#    ax.plot(X1[8], Y1[8] , 'bo', markersize=ms*0.75,linewidth=2.0, label = 'M31')

    ax.set_ylabel(Ylabel,fontsize = label_font_size)
    minorLocator   = AutoMinorLocator(5)
    ax.xaxis.set_minor_locator(minorLocator)
    ax.tick_params(axis='both', which='major', labelsize=20, length=10)
    ax.tick_params(axis='both', which='minor', labelsize=20, length=7, color='k')
    ax.legend( loc='best',prop={'size': legend_size} )

    plt.draw()
    return

  
def make_figure_12(engel_tab, m31_tab, feature_list, ax=None, xquant = '12plogOH', xlab = "12+ log[O/H]"):
    symlist = ['o','s']
    collist = ['b','r']
    # get data
    Oxy = engel_tab[xquant]
    Oxy_unc = engel_tab[xquant+'_unc']
    EQW= engel_tab['PAH8eqw']
    EQW_unc = engel_tab['PAH8eqw_unc']
    
    # m31 data, except for IRC3
    X = m31_tab[xquant][m31_tab['ID']!='irc3'] - 0.35
    Xerr = m31_tab[xquant+'_unc'][m31_tab['ID']!='irc3']

    if ax == None:
        fig, ax = plt.subplots()
    plt.errorbar(Oxy,EQW,EQW_unc,Oxy_unc,'o',color = '0.75', linewidth=2.0)
    plt.plot(Oxy,EQW,'o',mfc = 'white', markersize=ms*0.75, label='E08: PAH8')

    # loop ovr features to be plotted
    for i,feat in enumerate(feature_list):
        feature_lab = 'M31 '+ feat[:string.find(feat,'eqw')]
        Y = m31_tab[feat][m31_tab['ID']!='irc3']
        Yerr = m31_tab[feat+'_unc'][m31_tab['ID']!='irc3']
        
        ax.errorbar(X,Y, Yerr, Xerr,symlist[i],color = collist[i], linewidth=2.0)
        ax.plot(X,Y,symlist[i],color = collist[i], markersize=ms*0.75,label = feature_lab)

#   plot formatting
    minorLocator   = MultipleLocator(5)
    ax.xaxis.set_minor_locator(MultipleLocator(4))
    ax.yaxis.set_minor_locator(minorLocator)
    
    ax.tick_params(which='both', width=2)
    ax.tick_params(which='major', length=10)
    ax.tick_params(which='minor', length=7, color='k')
    ax.xaxis.set_minor_locator(MultipleLocator(5))

    ax.set_xlabel(xlab ,fontsize = label_font_size)
    ax.set_ylabel(r'EQW ($\mathregular{\mu m}$)' ,fontsize = label_font_size)
    ax.legend(loc='best' ,prop={'size': legend_size} )
    ax.set_yscale('log')
    ax.set_xlim(7.8,9.1)
    plt.draw()

    return
	
def combined_fig(m31dat, gord_dat, eng_dat):    
    fig,ax = plt.subplots(1,2,sharey=True)

    symlist = ['o','s']
    collist = ['b','r']

    ax[1].errorbar(eng_dat['12plogOH'],eng_dat['PAH8eqw'],eng_dat['PAH8eqw_unc'],eng_dat['12plogOH_unc'],'o',color = '0.75', linewidth=2.0)
    ax[1].plot(eng_dat['12plogOH'],eng_dat['PAH8eqw'],'p',mfc = 'white', markersize=ms*0.75, label='E08: PAH8')

    
    # m31 data, except for IRC3
    X = m31dat['12plogOH'][m31dat['ID']!='irc3'] - 0.35
    Xerr = m31dat['12plogOH_unc'][m31dat['ID']!='irc3']
    Y = m31dat['PAH8eqw'][m31dat['ID']!='irc3']
    Yerr = m31dat['PAH8eqw_unc'][m31dat['ID']!='irc3']
    ax[1].errorbar(X,Y, Yerr, Xerr,'o',color = 'w', linewidth=2.0)
    ax[1].plot(X,Y,'o',color = 'w', markersize=ms*0.75,label = 'M31')

    # loop ovr features to be plotted
    for i,feat in enumerate(['PAH7.7eqw','PAH11.3eqw']):
        feature_lab = 'M31 '+ feat[:string.find(feat,'eqw')]
        Y = m31dat[feat][m31dat['ID']!='irc3']
        Yerr = m31dat[feat+'_unc'][m31dat['ID']!='irc3']
        
        ax[1].errorbar(X,Y, Yerr, Xerr,'o',color = collist[i], linewidth=2.0)
        ax[1].plot(X,Y,'o',color = collist[i], markersize=ms*0.75,label = feature_lab)

        ax[1].errorbar(gord_dat['12plogOH'],gord_dat[feat],gord_dat[feat+'_unc'],gord_dat['12plogOH_unc'],'o',color = '0.75', linewidth=2.0)
        ax[1].plot(gord_dat['12plogOH'],gord_dat[feat],'s',mfc = collist[i], markersize=ms*0.75, label='G08')


    ax[1].set_xlabel("12+ log[O/H]" ,fontsize = label_font_size)
#    ax[1].set_ylabel(r'EQW ($\mathregular{\mu m}$)' ,fontsize = label_font_size)
    ax[1].legend(loc='best' ,prop={'size': legend_size} )
    ax[1].set_yscale('log')
    ax[1].set_xlim(7.8,9.1)
    plt.draw()
    

    return
