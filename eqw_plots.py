import numpy as np                     
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import math

##### SECTION: calculate stuff


def getIIforEngel():
    Line= np.loadtxt("englbrt II")
    SIV_SIII = Line[:,8]/Line[:,6]
    NeIII_NeII = Line[:,4]/Line[:,2]
    II = (np.log10((np.array([NeIII_NeII]))) + ( 0.71 + (1.58*(np.log10(np.array(SIV_SIII))))))/2
    meanII = mean(II)
    return (II)


def getIIerrorforEngel(NeII,NeIII,SIII,SIV,NeIIval,NeIIIval,SIIIval,SIVval,IIval):
    """ Calculates uncertainties for II, NeIII/NeII and SIV/SIII
    Inputs : first four parameters are the uncertainties of atomic lines and next four are their values. IIval has II values
             All are 1D arrays.

    Output : Uncertainties of II values."""
    
    IIerr = []
    for i in range(0,(np.size(SIV))):

        delS = math.sqrt(((SIV[i]/SIVval[i])**2) + ((SIII[i]/SIIIval[i])**2))*(SIVval[i]/SIIIval[i])
        delNe = math.sqrt(((NeIII[i]/NeIIIval[i])**2) + ((NeII[i]/NeIIval[i])**2))*(NeIIIval[i]/NeIIval[i])
        delII = math.sqrt((delS**2)+(delNe**2))
        #print delII
        IIerr.append((delII/ (10**(IIval[i])))*0.434 )
    return IIerr

def getIIforGordon():
    #Calculating II and returns II and itz uncertainties.
    Line= np.loadtxt("gor_atom", skiprows=2)
    delNeII,delSII,IIerror = getIIerrorforGordon(Line[:,1],Line[:,3], Line[:,5],Line[:,7],Line[:,0],Line[:,2], Line[:,4],Line[:,6])
    
    SIV_SIII = Line[:,0]/Line[:,6]
    NeIII_NeII = Line[:,4]/Line[:,2]

    II = ((np.log10((np.array([NeIII_NeII])))/(delNeII**2) + ( 0.71 + (1.58*(np.log10(np.array(SIV_SIII)))))/(delSII**2))/(delSII**(-2)+ delNeII**(-2)))
    # There are two missing emission lines in 1st and second rows. 1st is SIV, Second is NeIII. I used the alternative method to get II.
    II[0,0] = np.log10(Line[0,4]/Line[0,2])
    II[0,1] = ( 0.71 + (1.58*(np.log10(np.array(Line[1,0]/Line[1,6])))))    
    return II,IIerror
   
def getIIerrorforGordon(NeII,NeIII,SIII,SIV,NeIIval,NeIIIval,SIIIval,SIVval):
    """ Calculates uncertainties for II, NeIII/NeII and SIV/SIII
    Inputs : first four parameters are the uncertainties of atomic lines and next four are their values.
             All are 1D arrays.

    Output : II values and their errors.

    Notes  : More details about calculating RHI (II) values can be found in the paper.
             We are using different methods to calculate II according to the line they are missing
             SIII is not missing in any region. Generating NAN values will not matter beacause they will not be plotted"""
    
    
    IIerr = []
    delNeII = []
    delSII = []
    delII = []
    for i in range(0,(np.size(SIV))):

        delS = math.sqrt(((SIV[i]/SIVval[i])**2) + ((SIII[i]/SIIIval[i])**2))*(SIVval[i]/SIIIval[i])
        dellogS = (delS/(SIVval[i]/SIIIval[i]))*0.434
        delSII.append(dellogS)
        delNe = math.sqrt(((NeIII[i]/NeIIIval[i])**2) + ((NeII[i]/NeIIval[i])**2))*(NeIIIval[i]/NeIIval[i])
        dellogNe = (delNe/(NeIIIval[i]/NeIIval[i]))*0.434
        delNeII.append(dellogNe)
        
        delII.append(math.sqrt((dellogS**2)+(dellogNe**2)))
        
    #1st two regions of Gordon et al. 2008 was missing some atomic lines. so we set the II to that calculated by Ne and S lines
    delII[0] = delNeII[0] 
    delII[1] = delSII[1]
    #print delII
    return (np.array(delNeII)),(np.array(delSII)),delII  


def norm_gordon():
    # Finding the normalization factor for EQWs and normalize both EQW and Uncertainty
    
    Norm_Fact_Arr= [0.87,0.87,3.01,3.01,1.04,1.04,1.32,1.32,0.56,0.56] # These are given in Gordon et al. 2008
    Norm_Fact_Arr = 1/ np.array(Norm_Fact_Arr)
    eqw = np.loadtxt("gordonEQW", skiprows = 2)
    
    i = 0
    # Thi sgoes through all 5 PAH features. It skips the uncertainty by saying i = i + 2
    while i <10 :
        eqw[:,i] = (eqw[0:,i])*(Norm_Fact_Arr[i])
        i = i + 2
    np.savetxt('Norm_Gord_EQW.txt', eqw,fmt='%.2f')
    return

def compute_rhi(atomic_lines, atomic_line_unc):
    """ Calculate single RHI value
    Inputs: atomic_lines: tuple with SIV,SIII,NeIII, NeII fluxes
            atomic_line_unc: tuple with SIV,SIII,NeIII, NeII flux uncertainties
            (note convention for upper limits here)

    Output : RHI value and uncertainty."""

    SIV, SIII,  NeIII, NeII = atomic_lines
    SIVunc, SIIIunc, NeIIIunc, NeIIunc = atomic_line_unc
    
    # TODO: fix calculation
    delNe = math.sqrt(((NeIIIunc/NeIII)**2) + ((NeIIunc/NeII)**2))*(NeIII/NeII)
    dellogNe = (delNe/(NeIII/NeII))*0.434

    #Calculating uncertainty for log(SIV/SIII)
    delS = math.sqrt(((SIVunc/SIV)**2) + ((SIIIunc/SIII)**2))*(SIV/SIII)
    dellogS = (delS/(SIV/SIII))*0.434

    #Calculating uncertainty for II
    delII = (math.sqrt((dellogS**2)+(dellogNe**2)))
    #delII = (dellogS + dellogNe)/2
    
    if SIV == 0.000 :
    # If any Ne line is missing then II value will be infinite. But they will not be plotted. Don't care about division by zero of Log(0).
    # Same for the other if conditions.
        II = np.log10(NeIII/NeII)
        IIerror = dellogNe
        
        
    if NeII == 0.00000 or NeIII == 0.00000 :

        RHI = ( 0.71 + (1.58*(np.log10(SIV/SIII))))
        RHI_unc = dellogS
        
    if SIV != 0.00 and NeII != 000 and NeIII != 0.00 :

        #II = ((np.log10(NeIII/NeII))/(dellogNe**2) + ( 0.71 + (1.58*(np.log10(SIV/SIII))))/(dellogS**2))/(dellogS**(-2)+ dellogNe**(-2))
        RHI = ((np.log10(NeIII/NeII)) + ( 0.71 + (1.58*(np.log10(SIV/SIII)))))/2
        RHI_unc = delII
   
    return(RHI,RHI_unc)


# SECTION: make plots

def make_fig_10_plot():
    II = getIIforEngel()
    np.savetxt('II_values.txt', II)
    EQW= np.loadtxt("englbrt II") # Here EQW does not have onlt EQWs, it also has atomic lines and uncertainties.
                                  # Order of the columns are the same as table 5 in Engelbracht et al 200
    
    # X is II values and Y is EQWs. These are going to be X and Y axis of the plot.
    X= np.loadtxt("II_values.txt")
    
    Y = list(EQW[0:,0])
    Yerr = list(EQW[0:,1])
    Yerr = (np.array(Yerr)/np.array(Y))*0.434
    Y = [np.log10(a) for a in Y]
    
    IIerr = getIIerrorforEngel(EQW[:,3],EQW[:,5], EQW[:,7],EQW[:,9],EQW[:,2],EQW[:,4], EQW[:,6],EQW[:,8], X)
    Xerr = IIerr
    plt.errorbar(X,Y,Yerr,Xerr,'ks',linewidth=2.0)
    plt.plot(X,Y,'ks',label = 'Engelbracht et al. 2008', markersize=15 , mfc = 'white')
       
    feature = 0
    X1,X1err,Y1,Y1err = getIIvsEQW(feature)
    plt.plot(X1[1], Y1[1], 'b<', markersize=20,linewidth=2.0) #Upper Limit
    plt.errorbar(X1[3], Y1[3] ,Y1err[3],X1err[3], 'bo', markersize=15,linewidth=2.0)
    plt.plot(X1[4], Y1[4] , 'b>', markersize=20,linewidth=2.0)
    plt.errorbar(X1[5], Y1[5] ,Y1err[5],X1err[5], 'bo', markersize=15,linewidth=2.0)
    plt.errorbar(X1[6], Y1[6] ,Y1err[6],X1err[6], 'bo', markersize=15,linewidth=2.0)
    plt.plot(X1[7], Y1[7] , 'b<', markersize=20,linewidth=2.0) # Lowe limit
    plt.errorbar(X1[8], Y1[8] ,Y1err[8],X1err[8], 'bo', markersize=15,linewidth=2.0)
    plt.plot(X1[8], Y1[8] , 'bo', markersize=15,linewidth=2.0, label = 'M31')
    
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

def make_figure_11():    
    fig,axes = plt.subplots(2,1,sharex=True)
    
    plotting(2,0,"Log(EQW_7.7($\mu m$))")
    plotting(6,1,"Log(EQW_11.3($\mu m$))")
    return


def plotting(feature,panel,Ylabel):
    """ Plotting EQWs vs II for M31 regions and also over plotting it with Gordon et al. 2008.
    Uses the fnction fig8Gordon() 

    X values : II values
    Y values : EQWs    """ 

    X1,X1err,Y1,Y1err = getIIvsEQW(feature)  # M31 data
    X2,X2err,Y2,Y2err = fig8Gordon(feature)  # Gordon's data

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


# Feature = which feature u want to plot with II. 0 = 6.2, 2 = 7.7, 4 = 8.6, 6 = 11.3, 8 = 12.7
def fig8Gordon(feature):
    # Returns II and EQWs with their uncertainties. All are in log values. uncertainties are calculated to deal with log scale.
    II,IIerror = getIIforGordon()
    np.savetxt('II_values_gord.txt', II)
    II= np.loadtxt("II_values_gord.txt")
    EQWvals= np.loadtxt("Norm_Gord_EQW.txt")
   
    EQW = (EQWvals[0:,feature])
    EQWerr = (EQWvals[0:,feature+1])
  
    EQWerr = ((EQWerr/EQW)*0.434)
    EQW = [np.log10(a) for a in EQW]

    return II,IIerror,EQW,EQWerr
  
def make_figure_12():
    fig, ax = plt.subplots()
    minorLocator   = MultipleLocator(5)
    ax.xaxis.set_minor_locator(MultipleLocator(4))
    ax.yaxis.set_minor_locator(minorLocator)
    
    plt.tick_params(which='both', width=2)
    plt.tick_params(which='major', length=10)
    plt.tick_params(which='minor', length=7, color='k')
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    englbrt()
    Oxy = np.loadtxt("Oxy_abun.dat", skiprows = 1,dtype="string,f,f")  # Reading the metallicity values
    Oxy = zip(*Oxy)
    
    Combined_EQW = np.loadtxt("EQW.txt") # Reading the normalized EQWs from the file.
    eqw_unc = np.loadtxt("EQW_unc.txt") # Reading uncertainties of the normalized EQWs from the file.
    
    feature = 2  # Which feature you want to plot ? 2 = 7.7 mu feature
    # 0 to 9 represent 5.7, 6.2, 7.7, 8.3, 8.6, 10.7, 11.3, 12.0, 12.7, 17.0 mu dust features.
    
    # Removing IRC3 data.
    X = np.delete(Oxy[1],[2])
    Xerr = np.delete(Oxy[2],[2])
    
    Y = np.delete(Combined_EQW[0:,feature],[2])
    Yerr = np.delete(eqw_unc[0:,feature],[2])
    
    
    getplot_oxy_eqw(Y,X,Xerr,Yerr,'7.7 AF','ro','r') # Change the 5th argument according to the dust feature you plot.
    
    feature = 6  # 11.3 dust feature
    
    Y = np.delete(Combined_EQW[0:,feature],[2])
    Yerr = np.delete(eqw_unc[0:,feature],[2])
    
    getplot_oxy_eqw(Y,X,Xerr,Yerr,'11.3 AF','bs','b')
    return
	

def englbrt():
    """ Plots the Normalized EQWs vs Metallicity for the starburst galaxy sample in Engelbracht et al. 2008.
    This reads a table which has EQW, EQW uncertainty, Oxygen abundance and its uncertainty in four columns."""

    table = np.loadtxt("englbrt_eqw_oxy", skiprows = 2)
    EQW = table[0:,0]  # REading EQWs from the file
    EQWunc = table[0:,1] # REading uncertainties of EQWs from the file
    normalizedEQW = EQW/(np.mean(EQW))   # Normalizing EQWs
    normalizedEQWunc = EQWunc/(np.mean(EQW))  # Normalizing EQW uncertainties 
    normalizedEQWunc = (normalizedEQWunc/normalizedEQW) * 0.434  # Chanding uncertainty to deal with log scale
    Oxy = table[0:,2]  # Reading metallicity values
    Oxyunc = table[0:,3] # Reading uncertainties of metallicity values 
    plt.errorbar(Oxy,normalizedEQW,normalizedEQWunc,Oxyunc,'o',color = '0.75', linewidth=2.0)
    plt.plot(Oxy,normalizedEQW,'o',mfc = 'white', markersize=15)
    plt.yscale('log')
    plt.legend(loc='lower right',prop={'size':20})
    plt.show()
    return
    

#Plots eqw vs oxygen abundance for a given PAH feature
def getplot_oxy_eqw(eqw,oxy,oxy_err,eqw_err,featureL,symbol,col):
    """ Plots the normalized EQWs vs Metallicity.
    Inputs : eqw = 2D array of EQWs (columns) for each region (rows)
             eqw_err = Uncertainty values of the above.
             oxy = An array of oxygen abundance values. (Length should match the number of rows in "eqw")
             oxy_err = Uncertainty values of the above.
             featureL = A string to get the lable of the plot. eg: '11.3 mu' for 11.3 mu feature
             col = A string to get the color of the plot. eg: 'r' for red

    Outputs: Plot of normalized EQWs vs Metallicity.

    Notes : IRC3 is removed from the plot."""
    
    oxy_eqw = []
    for i in range(np.size(eqw)):
        if eqw[i]!=0:
            oxy_eqw.append([oxy[i],eqw[i],oxy_err[i],eqw_err[i]])
    oxy_eqw = np.array(oxy_eqw)
    Yerr = (oxy_eqw[0:,3]/(10**(oxy_eqw[0:,1]))) * 0.434
    Xerr = oxy_eqw[0:,2]
    # We applied a shift of -0.35 for metalicities from M31.
    plt.errorbar((oxy_eqw[0:,0] - 0.35), oxy_eqw[0:,1], Yerr, Xerr,symbol,color = col, linewidth=2.0)
    plt.plot((oxy_eqw[0:,0] - 0.35), oxy_eqw[0:,1],symbol,color = col, markersize=15,label = featureL)
    plt.yscale('log')
    plt.xlabel(" log(O/H) + 12" ,fontsize=28)
    plt.ylabel("EQW ($\mu m$)" ,fontsize=28)
    plt.legend( loc='lower left' ,prop={'size':20} )
    plt.plot(8.058,0.072, 'ws', mec='white') # Creates a white spot for an errorbar from an outlier
    plt.show()
    return
