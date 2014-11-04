"""
Purpose : Plot normalized EQW vs RHI for M31 sample and overplot with that of Gordon et al 2008.

Region  : All regions of M31 except the nucleus. 

Inputs  : 1) "PAHfilenames.dat" file with file names and the corresponding files (Obtained from spec.pro IDL code)
          3) "myatomic new.txt" and "lineerror new.txt" data files

Outputs : Specified at each function below

Notes  : You should have gordonII.py code in the same directory.

"""


import numpy as np                     


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


