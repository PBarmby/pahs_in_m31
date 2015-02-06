import numpy as np                     
import matplotlib.pyplot as plt
from pylab import *
from matplotlib.ticker import AutoMinorLocator

"""
This program calculates the ionization index (II) for M31 data and Engelbracht's data for starbusrst galaxies.
It alos plots the EQW of PAH vs II for both data sets.
Documentation is similar to IIeqw.py

"""

def getII(Line, Lineunc):
    """ Calculates RHI values(here I call them II values) for M31 data.
    Inputs : Line :2D matrices with SIV,NeII,NeIII, SIV fluxes in the same order
             Lineunc :2D matrices with SIV,NeII,NeIII, SIV flux uncertainties in the same order

    Output : II values and their errors.

    Notes  : More details about calculating RHI (II) values can be found in the paper.
             We are using different methods to calculate II according to the line they are missing
             SIII is not missing in any region. Generating NAN values will not matter beacause they will not be plotted"""

    SIV,NeII,NeIII,SIII = Line[0], Line[1], Line[2], Line[3]
    SIVunc,NeIIunc,NeIIIunc,SIIIunc = Lineunc[0], Lineunc[1], Lineunc[2], Lineunc[3]
    
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

        II = ( 0.71 + (1.58*(np.log10(SIV/SIII))))
        IIerror = dellogS
        
    if SIV != 0.00 and NeII != 000 and NeIII != 0.00 :

        #II = ((np.log10(NeIII/NeII))/(dellogNe**2) + ( 0.71 + (1.58*(np.log10(SIV/SIII))))/(dellogS**2))/(dellogS**(-2)+ dellogNe**(-2))
        II = ((np.log10(NeIII/NeII)) + ( 0.71 + (1.58*(np.log10(SIV/SIII)))))/2
        IIerror = delII
   
    return II,IIerror



###########################################################################################


def getIIvsEQW(feature) :
    """ Get II values for all the regions and return them with corresponding EQWs.
    Inputs : Feature : Integre that specifies which PAH feature we want to plot with II. (see the note below)

    Output : Arrays of II, II uncertainties, EQW and EQW uncertainties

    Notes  : Feature numbers :0 = combined 7.7, 8.3 and 8.6 , 1 = 6.2, 2 = 7.7 , 6 = 11.3 , 8 = 12.7 PAH features
             This reads 4 files. Myatomicwithselecteduperlimits.txt, AtomicLines_unc2.txt, EQW_combined.dat, eqw_unc.dat """
    

 
    Line= np.loadtxt("Myatomicwithselecteduperlimits.txt", skiprows=1)  # Line intensities of SIV, NeII, NeII, NeIII
    Lineunc= np.loadtxt("AtomicLines_unc2.txt",skiprows=1) # Line uncertainties of SIV, NeII, NeII, NeIII

    II = []
    IIerror = []
    numberofregions = np.shape(Line)[0]
    for i in range(0,numberofregions) :

        IIval,IIerrval = getII(Line[i,:],Lineunc[i,:])
        II.append(IIval)
        IIerror.append(IIerrval)

    #print II, IIerror

    EQW= np.loadtxt("EQW_combined.dat")  # Reading EQW values of PAH liness
    EQWerr= np.loadtxt("eqw_unc.dat")   # Reading Uncertainties of PAH lines
    print np.shape(EQW)
    #feature = 8   # Which PAH feature you want to plot against II ? 0 = combined 7.7, 8.3 and 8.6 ,2 = 7.7, 6 = 11.3 , 8 = 12.7 micron features

    if feature == 0:
        Y = (EQW[0:,2]+EQW[0:,3]+EQW[0:,4])   #/ EQW[0:,6]
        Yerr = list(EQWerr[0:,2]+EQWerr[0:,3]+EQWerr[0:,4])
        Yerr = (np.array(Yerr)/np.array(Y))*0.434
        Y = [np.log10(a) for a in Y]
        ylable = "Log(EQW_8($\mu m$))"

    else :
        Y = list(EQW[0:,feature])
        Yerr = list(EQWerr[0:,feature])
        Yerr = (np.array(Yerr)/np.array(Y))*0.434
        Y = [np.log10(a) for a in Y]
        ylable = "Log(EQW_11.3($\mu m$))"   # Change this according to the feature


    
    II = np.delete(II,[9])  # Removing region 9 from the sample 
    IIerror = np.delete(IIerror,[9])
  
    np.savetxt('II_values_M31.txt', II)

    return II,IIerror,Y,Yerr

def getIIerror(SIV,NeII,NeIII,SIII,Line,II):
#Calculates uncertainties for II
#Inputs: Uncertainties of four lines as SIV,NeII,NeIII,SIII
#      : Table of values with line intensities as Line
#      : Ionisation Index as II
#Output: Uncertainties of II
#Theory1 if f = x/y then del(f) = sqrt( (del(x)/x)^2 + (del(y)/y)^2 ) * y
#Theory2 if f = x+y then del(f) = sqrt( del(x)^2 + del(y)^2 )
#Theory3 if f = log(x) the del(f) = (del(x)/x)*0.434

    IIerr = []
    #Line= np.loadtxt("myatomic II.txt", skiprows=3)
    Line= np.loadtxt("Myatomicwithselecteduperlimits.txt", skiprows=1)
    SIVval,NeIIval,NeIIIval,SIIIval = Line[:,0],Line[:,1],Line[:,2],Line[:,3]
    SIV_SIII = Line[:,0]/Line[:,3]
    NeIII_NeII = Line[:,2]/Line[:,1]
    for i in range(0,(np.size(SIV))):

        delS = math.sqrt(((SIV[i]/SIVval[i])**2) + ((SIII[i]/SIIIval[i])**2))*(SIV_SIII[i])  # Theory1
        delNe = math.sqrt(((NeIII[i]/NeIIIval[i])**2) + ((NeII[i]/NeIIval[i])**2))*(NeIII_NeII[i]) # Theory1
        delII = math.sqrt((delS**2)+(delNe**2)) # Theory2
        IIerr.append((delII/ (10**(II[i])))*0.434 ) # Theory3. Since II values are already in log scale it is taken back too base 10
    return IIerr



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


def getIIforEngel():

    Line= np.loadtxt("englbrt II")

    SIV_SIII = Line[:,8]/Line[:,6]
    NeIII_NeII = Line[:,4]/Line[:,2]
    #print NeIII_NeII,SIV_SIII
    
    #print np.array([lgS])*1.58
    II = (np.log10((np.array([NeIII_NeII]))) + ( 0.71 + (1.58*(np.log10(np.array(SIV_SIII))))))/2
    #II = ((np.log10(NeIII/NeII))/(dellogNe**2) + ( 0.71 + (1.58*(np.log10(SIV/SIII))))/(dellogS**2))/(dellogS**(-2)+ dellogNe**(-2))
    meanII = mean(II)
    #return np.log10(II)
    return (II)
    #II = II/float(meanII)


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
#plt.errorbar(X1,Y1,Y1err,X1err,'bo')

#plt.plot(X1[0], Y1[0],Y1err[0],X1err[0], 'bo', markersize=15,linewidth=2.0,)  
plt.plot(X1[1], Y1[1], 'b<', markersize=20,linewidth=2.0) #Upper Limit
#plt.plot(X1[2], Y1[2],Y1err[2],X1err[2], 'bo', markersize=15,linewidth=2.0) 
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
#ax.yaxis.set_minor_locator(minorLocator)

plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=10)
plt.tick_params(which='minor', length=7, color='k')
plt.show()


"""
Purpose : Plot normalized EQW vs RHI for M31 sample and overplot with that of Gordon et al 2008.

Region  : All regions of M31 except the nucleus. 

Inputs  : 1) "PAHfilenames.dat" file with file names and the corresponding files (Obtained from spec.pro IDL code)
          3) "myatomic new.txt" and "lineerror new.txt" data files

Outputs : Specified at each function below

Notes  : You should have gordonII.py code in the same directory.

"""


import numpy as np                     
import matplotlib.pyplot as plt
from pylab import *
from gordonII import fig8Gordon  # Importing fig8Gordon function from "gordonII.py". 
from matplotlib.ticker import AutoMinorLocator


########################################################################################################################
def getII(Line, Lineunc):
    """ Calculates RHI values(here I call them II values) for M31 data.
    Inputs : Line :2D matrices with SIV,NeII,NeIII, SIV fluxes in the same order
             Lineunc :2D matrices with SIV,NeII,NeIII, SIV flux uncertainties in the same order

    Output : II values and their errors.

    Notes  : More details about calculating RHI (II) values can be found in the paper.
             We are using different methods to calculate II according to the line they are missing
             SIII is not missing in any region. Generating NAN values will not matter beacause they will not be plotted"""
    

    # Calculating uncertainty for log(NeIII/NeII)

    SIV,NeII,NeIII,SIII = Line[0], Line[1], Line[2], Line[3]
    SIVunc,NeIIunc,NeIIIunc,SIIIunc = Lineunc[0], Lineunc[1], Lineunc[2], Lineunc[3]
    
    delNe = math.sqrt(((NeIIIunc/NeIII)**2) + ((NeIIunc/NeII)**2))*(NeIII/NeII)
    dellogNe = (delNe/(NeIII/NeII))*0.434

    #Calculating uncertainty for log(SIV/SIII)
    delS = math.sqrt(((SIVunc/SIV)**2) + ((SIIIunc/SIII)**2))*(SIV/SIII)
    dellogS = (delS/(SIV/SIII))*0.434

    #Calculating uncertainty for II
    delII = (math.sqrt((dellogS**2)+(dellogNe**2)))
    
    if SIV == 0.000 :

        II = np.log10(NeIII/NeII)
        IIerror = dellogNe
        
        
    if NeII == 0.00000 or NeIII == 0.00000 :

        II = ( 0.71 + (1.58*(np.log10(SIV/SIII))))
        IIerror = dellogS
        
    if SIV != 0.00 and NeII != 000 and NeIII != 0.00 :

        II = ((np.log10(NeIII/NeII))/(dellogNe**2) + ( 0.71 + (1.58*(np.log10(SIV/SIII))))/(dellogS**2))/(dellogS**(-2)+ dellogNe**(-2))
        IIerror = delII
        

    return II,IIerror


##########################################################################################
# Normalizing by weighting over uncertainty
# Wi = 1/ eqw_uncertainty^2
# xi = eqw
#  X = sum( Wi * xi) / sum( Wi)  # sum over one region

def normalizedEQW(eqw,eqw_unc):
    #Normalize EQWs using the equations mentioned above. eqw and eqw_unc are the
    # 1D arrays of EQWs and their uncertainties.
    # Returns normalized EQW as an array.

    Num_regions,Num_features = np.shape(eqw_unc)[0],np.shape(eqw_unc)[1]
    #print eqw_unc[0:,1]
    X = []
    for i in range(Num_features) :
        Wi = []
        for j in range(Num_regions) :
    
            if eqw_unc[j,i] == 0:
                Wi.append(0)
            else :
                Wi.append(1./(float(eqw_unc[j,i]**2)))
        x = sum(Wi*eqw[0:,i]) / sum(Wi)
        X.append(x)
        #print len(x)
        eqw[0:,i] = eqw[0:,i]/x  

        
    return eqw

###########################################################################################


def getIIvsEQW(feature) :

    """ Get II values for all the regions and return them with corresponding EQWs.
    Inputs : Feature : Integre that specifies which PAH feature we want to plot with II. (see the note below)

    Output : Arrays of II, II uncertainties, EQW and EQW uncertainties

    Notes  : Feature numbers :0 = combined 7.7, 8.3 and 8.6 , 1 = 6.2, 2 = 7.7 , 6 = 11.3 , 8 = 12.7 PAH features
             This reads 4 files. Myatomicwithselecteduperlimits.txt, AtomicLines_unc2.txt, EQW_combined.dat, eqw_unc.dat """
    

    
    Line= np.loadtxt("Myatomicwithselecteduperlimits.txt", skiprows=1)  # Line intensities of SIV, NeII, NeII, NeIII
    # This file have only selected upper limit values for regions 2, 5 and 8.
    Lineunc= np.loadtxt("AtomicLines_unc2.txt",skiprows=1) # Line uncertainties of SIV, NeII, NeII, NeIII

    II = []
    IIerror = []
    numberofregions = np.shape(Line)[0]

    # This for loop calls the getII function and get II and itz uncertainty for all the regions in M31.
    for i in range(0,numberofregions) :

        IIval,IIerrval = getII(Line[i,:],Lineunc[i,:])
        II.append(IIval)
        IIerror.append(IIerrval)

    #print II, IIerror

    EQW= np.loadtxt("EQW_combined.dat")  # Reading EQW values of PAH liness
    EQWerr= np.loadtxt("eqw_unc.dat")   # Reading Uncertainties of PAH lines
    EQW = normalizedEQW(EQW,EQWerr)
    
    #feature = 8   # Which PAH feature you want to plot against II ? 0 = combined 7.7, 8.3 and 8.6 ,2 = 7.7, 6 = 11.3 , 8 = 12.7 micron features

    if feature == 0:
        Y = (EQW[0:,2]+EQW[0:,3]+EQW[0:,4])   #/ EQW[0:,6]
        Yerr = list(EQWerr[0:,2]+EQWerr[0:,3]+EQWerr[0:,4])
        Yerr = (np.array(Yerr)/np.array(Y))*0.434
        Y = [np.log10(a) for a in Y]
        ylable = "Log(EQW_8($\mu m$))"

    else :
        Y = list(EQW[0:,feature])
        Yerr = list(EQWerr[0:,feature])
        Yerr = (np.array(Yerr)/np.array(Y))*0.434
        Y = [np.log10(a) for a in Y]
        ylable = "Log(EQW_11.3($\mu m$))"   # Change this according to the feature

    return II,IIerror,Y,Yerr


########################################################################################################################
def plotting(feature,panel,Ylabel):
    """ Plotting EQWs vs II for M31 regions and also over plotting it with Gordon et al. 2008.
    This reads the fnction fig8Gordon() in gordonII code.

    X values : II values
    Y values : EQWs    """ 

    X1,X1err,Y1,Y1err = getIIvsEQW(feature)  # M31 data
    X2,X2err,Y2,Y2err = fig8Gordon(feature)  # Gordon's data

    axes[panel].errorbar(X2,Y2,Y2err,X2err,'ks', markersize=15,linewidth=2.0)
    axes[panel].plot(X2,Y2,'gs', markersize=15,linewidth=2.0, mfc = "white", label = 'Gordon et al. 2008')

    #axes[panel].plot(X1[0], Y1[0],Y1err[0],X1err[0], 'bo', markersize=15,linewidth=2.0,)  
    axes[panel].plot(X1[1], Y1[1], 'b<', markersize=20,linewidth=2.0) #Upper Limit
    #axes[panel].plot(X1[2], Y1[2],Y1err[2],X1err[2], 'bo', markersize=15,linewidth=2.0) # Remove Region 3 data point
    axes[panel].errorbar(X1[3], Y1[3] ,Y1err[3],X1err[3], 'bo', markersize=15,linewidth=2.0)
    axes[panel].plot(X1[4], Y1[4], 'b>', markersize=20,linewidth=2.0)
    axes[panel].errorbar(X1[5], Y1[5] ,Y1err[5],X1err[5], 'bo', markersize=15,linewidth=2.0)
    axes[panel].errorbar(X1[6], Y1[6] ,Y1err[6],X1err[6], 'bo', markersize=15,linewidth=2.0)
    axes[panel].plot(X1[7], Y1[7] ,'b<', markersize=20,linewidth=2.0)
    axes[panel].errorbar(X1[8], Y1[8] ,Y1err[8],X1err[8], 'bo', markersize=15,linewidth=2.0)
    axes[panel].plot(X1[8], Y1[8] , 'bo', markersize=15,linewidth=2.0, label = 'M31')


    
    #axes[panel].errorbar(X1,Y1,Y1err,X1err,'o')
    
    axes[panel].set_ylabel(Ylabel,fontsize=24)
    axes[1].set_xlabel("Log(RHI)",fontsize=24)
    axes[panel].tick_params(axis='both', which='major', labelsize=30)
    axes[panel].tick_params(axis='both', which='minor', labelsize=30)
    minorLocator   = AutoMinorLocator(5)
    #axes[0].xaxis.set_minor_locator(minorLocator)
    #axes[0].yaxis.set_minor_locator(minorLocator)
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
    
fig,axes = plt.subplots(2,1,sharex=True)

plotting(2,0,"Log(EQW_7.7($\mu m$))")
plotting(6,1,"Log(EQW_11.3($\mu m$))")

"""
Purpose : Plots normalized EQWs vs metallicity for both M31 data and Engelbracht et al. 2008 data.

Regions  : All regions of M31 except the nucleus, irc3.
          Starburst galaxy sample of engelbracht et al. 2008

Inputs  : 1) Two data files which
          2) Data files of EQWs of 18 dust features for each region.

Outputs : 1) Data file which has the uncertainties of EQWs of 10 dust features from all regions.
          2) Data file which has the EQWS of 11 dust features from all regions.
          3) Two data files of EQWs and their uncertainties.
          4) A table which has the EQW +/- Uncertainty in a readable version for Latex.

Notes  : Variable "feature"

"""



import numpy as np                     
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from pylab import *
import math

fig, ax = plt.subplots()
minorLocator   = MultipleLocator(5)
ax.xaxis.set_minor_locator(MultipleLocator(4))
ax.yaxis.set_minor_locator(minorLocator)

plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=10)
plt.tick_params(which='minor', length=7, color='k')


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

ax.xaxis.set_minor_locator(MultipleLocator(5))
englbrt()
    

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




    
        
      
    




